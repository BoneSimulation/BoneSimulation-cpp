#pragma once
#include "image_processing.hpp"
#include "mesh_generation.hpp"
#include <string>
#include <vector>

namespace bonesim {
    struct BlockExtractionConfig {
        double block_size_mm = 30.0;
        double step_size_mm = 30.0;
        double voxel_size_mm = 0.05;
        bool write_tetra_mesh = true;
        size_t min_voxels_threshold = 100;
        bool use_downsampling = true;
    };

    struct BlockInfo {
        size_t index;
        size_t z_start, y_start, x_start;
        size_t z_end, y_end, x_end;
        size_t voxel_count;
        bool mesh_success;
    };

    std::vector<BlockInfo> extract_blocks_parallel(const BinaryVolume& vol,
                                                   const BlockExtractionConfig& cfg,
                                                   const std::string& out_dir,
                                                   const std::string& timestamp,
                                                   int max_threads = 0);
}