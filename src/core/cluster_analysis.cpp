// cluster_analysis.cpp
#include "cluster_analysis.hpp"
#include <functional>
#include <numeric>
#include <iostream>

namespace bonesim {

// Einfache Union-Find Struktur
struct DSU {
    std::vector<size_t> parent, rank;
    DSU(size_t n) : parent(n), rank(n, 0) {
        std::iota(parent.begin(), parent.end(), 0);
    }
    size_t find(size_t x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    }
    void unite(size_t x, size_t y) {
        x = find(x); y = find(y);
        if (x == y) return;
        if (rank[x] < rank[y]) parent[x] = y;
        else if (rank[x] > rank[y]) parent[y] = x;
        else { parent[y] = x; rank[x]++; }
    }
};

ClusterResult find_largest_cluster(const BinaryVolume& vol, int connectivity) {
    size_t dim_z = vol.dim_z, dim_y = vol.dim_y, dim_x = vol.dim_x;
    size_t total = dim_z * dim_y * dim_x;
    DSU dsu(total);
    
    auto idx = [&](size_t z, size_t y, size_t x) { return (z * dim_y + y) * dim_x + x; };
    
    // Nachbarschaftsoffsets (26er für 3D)
    int offsets[13][3] = {
        {0,0,1}, {0,1,0}, {0,1,1}, {1,0,0}, {1,0,1}, {1,1,0}, {1,1,1},
        {0,1,-1}, {1,0,-1}, {1,1,-1}, {1,-1,0}, {1,-1,1}, {1,-1,-1}
    };
    
    for (size_t z = 0; z < dim_z; ++z) {
        for (size_t y = 0; y < dim_y; ++y) {
            for (size_t x = 0; x < dim_x; ++x) {
                if (vol.at(z, y, x) == 0) continue;
                size_t cur = idx(z, y, x);
                for (auto& off : offsets) {
                    int nz = z + off[0], ny = y + off[1], nx = x + off[2];
                    if (nz >= 0 && nz < (int)dim_z && ny >= 0 && ny < (int)dim_y && nx >= 0 && nx < (int)dim_x) {
                        if (vol.at(nz, ny, nx)) {
                            dsu.unite(cur, idx(nz, ny, nx));
                        }
                    }
                }
            }
        }
        if (z % 50 == 0) {
            std::cout << "Clustering slice " << z << "/" << dim_z << std::endl;
        }
    }
    
    std::vector<size_t> cluster_sizes(total, 0);
    for (size_t i = 0; i < total; ++i) {
        if (vol.data[i]) {
            cluster_sizes[dsu.find(i)]++;
        }
    }
    
    size_t largest_label = 0;
    size_t largest_size = 0;
    size_t num_clusters = 0;
    for (size_t i = 0; i < total; ++i) {
        if (cluster_sizes[i] > 0) {
            num_clusters++;
            if (cluster_sizes[i] > largest_size) {
                largest_size = cluster_sizes[i];
                largest_label = i;
            }
        }
    }
    
    BinaryVolume largest;
    largest.dim_z = dim_z; largest.dim_y = dim_y; largest.dim_x = dim_x;
    largest.data.resize(total, 0);
    for (size_t i = 0; i < total; ++i) {
        if (vol.data[i] && dsu.find(i) == largest_label) {
            largest.data[i] = 1;
        }
    }
    
    return {std::move(largest), num_clusters, largest_size};
}

ClusterResult find_largest_cluster_optimized(const BinaryVolume& vol, int connectivity) {
    return find_largest_cluster(vol, connectivity);
}

} // namespace bonesim