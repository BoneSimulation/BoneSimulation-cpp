// cluster_analysis.hpp
#pragma once

#include "image_processing.hpp"
#include <cstddef>

namespace bonesim {

    struct ClusterResult {
        BinaryVolume largest_cluster;
        size_t num_clusters;
        size_t largest_cluster_size;
    };

    ClusterResult find_largest_cluster(const BinaryVolume& vol, int connectivity = 1);
    ClusterResult find_largest_cluster_optimized(const BinaryVolume& vol, int connectivity = 1);

} // namespace bonesim