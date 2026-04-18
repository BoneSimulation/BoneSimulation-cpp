#pragma once
#include "cluster_analysis.hpp"
#include <string>
#include <vector>

namespace bonesim {
    struct Statistics {
        size_t num_clusters;
        double porosity;
        size_t largest_cluster_size;
        double mean_cluster_size;
        std::vector<size_t> cluster_sizes;
    };
    Statistics compute_statistics(const BinaryVolume& vol);
    void generate_pdf_report(const Statistics& stats, const std::string& output_path,
                             const std::string& title = "Bone Analysis Report");
}