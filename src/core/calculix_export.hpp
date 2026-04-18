#pragma once
#include "mesh_generation.hpp"
#include <string>

namespace bonesim {
    struct MaterialProperties { double E=18000.0; double nu=0.3; };
    void export_to_calculix(const Mesh& mesh, const MaterialProperties& mat,
                            const std::string& output_path,
                            bool fix_z_plane=true, bool load_z_plane=true,
                            double total_force=100.0);
    void generate_run_script(const std::string& inp_path, const std::string& solver_path,
                             const std::string& output_dir);
}