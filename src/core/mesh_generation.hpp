#pragma once
#include "image_processing.hpp"
#include <vector>
#include <array>
#include <string>

namespace bonesim {
    struct Mesh {
        std::vector<std::array<double,3>> vertices;
        std::vector<std::array<size_t,3>> faces;
        std::vector<std::array<size_t,4>> tetrahedra;
    };

    Mesh marching_cubes(const BinaryVolume& vol, double sx=1.0, double sy=1.0, double sz=1.0);
    Mesh generate_tetrahedral_mesh(const BinaryVolume& vol, double voxel_size,
                                   bool use_downsampling=true, size_t max_volume_size=100'000'000);
    void save_mesh_as_vtk(const Mesh& mesh, const std::string& filename);
    void save_surface_as_vtk(const Mesh& mesh, const std::string& filename);
    BinaryVolume reverse_binary(const BinaryVolume& vol);
}