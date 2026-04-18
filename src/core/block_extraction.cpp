#include "block_extraction.hpp"
#include <filesystem>
#include <fstream>
#include <omp.h>
#include <iostream>

namespace fs = std::filesystem;
namespace bonesim {

static size_t count_voxels(const BinaryVolume& vol, size_t z0,size_t z1,size_t y0,size_t y1,size_t x0,size_t x1) {
    size_t cnt=0;
    for (size_t z=z0; z<z1; ++z)
        for (size_t y=y0; y<y1; ++y)
            for (size_t x=x0; x<x1; ++x)
                if (vol.at(z,y,x)) ++cnt;
    return cnt;
}

static void compute_blocks(const BinaryVolume& vol, const BlockExtractionConfig& cfg, std::vector<BlockInfo>& blocks) {
    int block_vox = static_cast<int>(cfg.block_size_mm / cfg.voxel_size_mm);
    int step_vox  = static_cast<int>(cfg.step_size_mm / cfg.voxel_size_mm);
    if (block_vox<=0 || step_vox<=0) throw std::runtime_error("Invalid block/step size");
    blocks.clear();
    for (size_t z=0; z+block_vox<=vol.dim_z; z+=step_vox)
        for (size_t y=0; y+block_vox<=vol.dim_y; y+=step_vox)
            for (size_t x=0; x+block_vox<=vol.dim_x; x+=step_vox) {
                BlockInfo info;
                info.index = blocks.size();
                info.z_start=z; info.z_end=z+block_vox;
                info.y_start=y; info.y_end=y+block_vox;
                info.x_start=x; info.x_end=x+block_vox;
                info.voxel_count = count_voxels(vol,z,z+block_vox,y,y+block_vox,x,x+block_vox);
                info.mesh_success = false;
                blocks.push_back(info);
            }
}

static void save_block_npy(const BinaryVolume& parent, const BlockInfo& info, const std::string& path) {
    size_t dz = info.z_end - info.z_start;
    size_t dy = info.y_end - info.y_start;
    size_t dx = info.x_end - info.x_start;
    std::vector<uint8_t> data(dz*dy*dx);
    size_t idx=0;
    for (size_t z=info.z_start; z<info.z_end; ++z)
        for (size_t y=info.y_start; y<info.y_end; ++y)
            for (size_t x=info.x_start; x<info.x_end; ++x)
                data[idx++] = parent.at(z,y,x);
    std::ofstream f(path, std::ios::binary);
    f.write(reinterpret_cast<const char*>(data.data()), data.size());
}

static bool process_one_block(const BinaryVolume& parent, const BlockInfo& info,
                              const BlockExtractionConfig& cfg, const std::string& dir,
                              const std::string& ts) {
    setenv("TETGEN_MEMORY_LIMIT", "2048", 1);  // 2 GB Limit
    size_t dz = info.z_end - info.z_start;
    size_t dy = info.y_end - info.y_start;
    size_t dx = info.x_end - info.x_start;
    BinaryVolume block;
    block.dim_z=dz; block.dim_y=dy; block.dim_x=dx;
    block.data.resize(dz*dy*dx);
    for (size_t z=0; z<dz; ++z)
        for (size_t y=0; y<dy; ++y)
            for (size_t x=0; x<dx; ++x)
                block.at(z,y,x) = parent.at(info.z_start+z, info.y_start+y, info.x_start+x);
    std::string npy = dir + "/block_" + std::to_string(info.index) + "_" + ts + ".npy";
    save_block_npy(parent, info, npy);
    if (cfg.write_tetra_mesh && info.voxel_count >= cfg.min_voxels_threshold) {
        BinaryVolume inv = reverse_binary(block);
        BinaryVolume* mesh_vol = &inv;
        BinaryVolume down;
        if (cfg.use_downsampling && inv.size() > 50'000'000) {
            down = interpolate_volume(inv, 0.5);
            mesh_vol = &down;
        }
        std::string tet_path = dir + "/block_" + std::to_string(info.index) + "_tetra_" + ts + ".vtk";
        Mesh mesh = generate_tetrahedral_mesh(*mesh_vol, cfg.voxel_size_mm);
        bool ok = !mesh.tetrahedra.empty();
        if (ok) save_mesh_as_vtk(mesh, tet_path);
        return ok;
    }
    return false;
}

std::vector<BlockInfo> extract_blocks_parallel(const BinaryVolume& vol,
                                               const BlockExtractionConfig& cfg,
                                               const std::string& out_dir,
                                               const std::string& timestamp,
                                               int max_threads) {
    fs::create_directories(out_dir);
    std::vector<BlockInfo> blocks;
    compute_blocks(vol, cfg, blocks);
    if (max_threads <= 0) max_threads = omp_get_max_threads();
    omp_set_num_threads(max_threads);
    #pragma omp parallel for schedule(dynamic)
    for (size_t i=0; i<blocks.size(); ++i) {
        bool ok = process_one_block(vol, blocks[i], cfg, out_dir, timestamp);
        blocks[i].mesh_success = ok;
    }
    return blocks;
}

} // namespace