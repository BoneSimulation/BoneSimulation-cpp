#include "core/image_loading.hpp"
#include "core/image_processing.hpp"
#include "core/cluster_analysis.hpp"
#include "core/mesh_generation.hpp"
#include "core/block_extraction.hpp"
#include "core/calculix_export.hpp"
// #include "core/reporting.hpp"
#include "utils/logger.hpp"
#include <chrono>
#include <filesystem>
#include <omp.h>
#include <algorithm>

namespace fs = std::filesystem;

struct Config {
    std::string data_dir = "data/bigdataset";
    std::string output_dir = "pictures";
    std::string report_dir = "report";
    std::string calculix_dir = "calculix";
    size_t chunk_size = 50;
    double interpolation_factor = 0.1;
    bool dry_run = false;
    bool test_mode = false;
    bool enable_blocks = true;
    bool enable_calculix = true;
    int num_threads = 0;
};

void print_usage(const char* prog) {
    std::cout << "Usage: " << prog << " [options]\n"
              << "  --test dry             Dry run\n"
              << "  --test small           Small chunk size\n"
              << "  --test complete        Full pipeline\n"
              << "  --no-blocks            Skip block extraction\n"
              << "  --no-calculix          Skip CalculiX export\n"
              << "  --threads N            OpenMP threads\n"
              << "  --help                 Show help\n";
}

int main(int argc, char* argv[]) {
    Config cfg;
    for (int i=1; i<argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--test" && i+1<argc) {
            std::string m = argv[++i];
            if (m=="dry") cfg.dry_run = true;
            else if (m=="small") cfg.interpolation_factor = 0.1;
            else if (m=="complete") cfg.test_mode = true;
        }
        else if (arg == "--no-blocks") cfg.enable_blocks = false;
        else if (arg == "--no-calculix") cfg.enable_calculix = false;
        else if (arg == "--threads" && i+1<argc) cfg.num_threads = std::stoi(argv[++i]);
        else if (arg == "--help") { print_usage(argv[0]); return 0; }
    }
    if (cfg.num_threads>0) omp_set_num_threads(cfg.num_threads);
    auto start = std::chrono::steady_clock::now();

    try {
        bonesim::Logger::info("Loading images from " + cfg.data_dir);
        auto vol = bonesim::load_tiff_stack(cfg.data_dir);
        bonesim::Logger::info("Volume size: " + std::to_string(vol.dim_z) + "x" +
                              std::to_string(vol.dim_y) + "x" + std::to_string(vol.dim_x));

        // Otsu threshold berechnen
        uint16_t th = bonesim::compute_otsu_threshold(vol);
        bonesim::Logger::info("Otsu threshold = " + std::to_string(th));

        // OPTIMIERT: In-place Binärkonvertierung (spart eine komplette Volume-Kopie)
        bonesim::binary_from_volume_inplace(vol, th);
        bonesim::Logger::info("Binary conversion done (in-place)");

        // Konvertiere zu BinaryVolume (uint8_t) für weitere Operationen
        bonesim::BinaryVolume binary;
        binary.dim_z = vol.dim_z;
        binary.dim_y = vol.dim_y;
        binary.dim_x = vol.dim_x;
        binary.data.resize(vol.size());
        for (size_t i = 0; i < vol.size(); ++i) {
            binary.data[i] = static_cast<uint8_t>(vol.data[i]);
        }

        // Original Volume-Daten freigeben (spart ~50% Speicher)
        vol.data.clear();
        vol.data.shrink_to_fit();
        bonesim::Logger::info("Original volume freed");

        // OPTIMIERT: Morphological Closing (Slice-basiert)
        bonesim::Logger::info("Starting morphological closing...");
        auto closed = bonesim::morphological_closing_optimized(binary, 2);
        bonesim::Logger::info("Morphological closing done");

        // Binärvolume freigeben
        binary.data.clear();
        binary.data.shrink_to_fit();

        // OPTIMIERT: Interpolation mit Chunking
        bonesim::Logger::info("Starting interpolation (factor " + std::to_string(cfg.interpolation_factor) + ")...");
        auto interpolated = bonesim::interpolate_volume_chunked(closed, cfg.interpolation_factor, 32);
        // auto interpolated = closed; // only used when you want to use thewhole volume directly
        bonesim::Logger::info("Interpolation done");

        // Closed freigeben
        closed.data.clear();
        closed.data.shrink_to_fit();

        // OPTIMIERT: Cluster-Analyse
        bonesim::Logger::info("Finding largest cluster...");
        auto cluster_res = bonesim::find_largest_cluster_optimized(interpolated);
        if (cluster_res.num_clusters == 0) throw std::runtime_error("No clusters found");
        bonesim::Logger::info("Largest cluster size: " + std::to_string(cluster_res.largest_cluster_size));

        std::string ts = std::to_string(std::time(nullptr));

        // Surface Mesh
        if (!cfg.dry_run) {
            fs::create_directories(cfg.output_dir);
            std::string mesh_path = cfg.output_dir + "/mesh_" + ts + ".vtk";
            bonesim::Logger::info("Generating surface mesh...");
            auto surf_mesh = bonesim::marching_cubes(cluster_res.largest_cluster, 0.05, 0.05, 0.05);
            bonesim::save_surface_as_vtk(surf_mesh, mesh_path);
            bonesim::Logger::info("Surface mesh saved: " + mesh_path);
        }

        // Block-Extraktion (nur 1 Thread, kleine Blöcke)
        if (cfg.enable_blocks && !cfg.dry_run) {
            bonesim::BlockExtractionConfig block_cfg;
            block_cfg.block_size_mm = 8.0;        // 160 Voxel
            block_cfg.step_size_mm = 8.0;         // kein Overlap
            block_cfg.voxel_size_mm = 0.05;
            block_cfg.write_tetra_mesh = true;
            block_cfg.min_voxels_threshold = 50;
            std::string block_dir = cfg.output_dir + "/blocks_" + ts;

            bonesim::Logger::info("Starting block extraction (1 thread, block size 8mm)...");
            auto blocks = bonesim::extract_blocks_parallel(cluster_res.largest_cluster, block_cfg,
                                                           block_dir, ts, 1);
            size_t succ = std::count_if(blocks.begin(), blocks.end(), [](auto& b){ return b.mesh_success; });
            bonesim::Logger::info("Blocks extracted: " + std::to_string(blocks.size()) + ", meshes: " + std::to_string(succ));
        }

        // CalculiX Export
        if (cfg.enable_calculix && !cfg.dry_run && cfg.test_mode) {
            bonesim::Logger::info("Generating tetrahedral mesh for CalculiX...");
            auto inv = bonesim::reverse_binary(cluster_res.largest_cluster);
            auto tetra = bonesim::generate_tetrahedral_mesh(inv, 0.05);
            if (!tetra.tetrahedra.empty()) {
                fs::create_directories(cfg.calculix_dir);
                std::string inp = cfg.calculix_dir + "/bone_" + ts + ".inp";
                bonesim::MaterialProperties mat{18000.0, 0.3};
                bonesim::export_to_calculix(tetra, mat, inp);
                bonesim::generate_run_script(inp, "ccx", cfg.calculix_dir);
                bonesim::Logger::info("CalculiX export done");
            } else {
                bonesim::Logger::warning("Tetrahedral mesh generation failed");
            }
        }

        if (cfg.test_mode && !cfg.dry_run) {
            bonesim::Logger::info("Report generation disabled (libharu missing)");
        }

    } catch (const std::exception& e) {
        bonesim::Logger::error(e.what());
        return 1;
    }

    auto end = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(end-start).count();
    bonesim::Logger::info("Total time: " + std::to_string(elapsed) + " s");
    return 0;
}