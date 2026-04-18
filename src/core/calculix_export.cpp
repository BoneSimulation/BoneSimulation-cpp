#include "calculix_export.hpp"
#include <fstream>
#include <algorithm>
#include <filesystem>
#include <omp.h>

namespace fs = std::filesystem;
namespace bonesim {

void export_to_calculix(const Mesh& mesh, const MaterialProperties& mat,
                        const std::string& output_path,
                        bool fix_z_plane, bool load_z_plane,
                        double total_force) {
    if (mesh.tetrahedra.empty()) throw std::runtime_error("No tetrahedra");
    double min_z = mesh.vertices[0][2], max_z = mesh.vertices[0][2];
    for (const auto& v : mesh.vertices) {
        min_z = std::min(min_z, v[2]); max_z = std::max(max_z, v[2]);
    }
    const double eps = 1e-6;
    std::ofstream f(output_path);
    if (!f) throw std::runtime_error("Cannot write " + output_path);
    f << "*HEADING\nBone simulation\n*NODE, NSET=ALLNODES\n";
    for (size_t i=0; i<mesh.vertices.size(); ++i)
        f << i+1 << ", " << mesh.vertices[i][0] << ", " << mesh.vertices[i][1] << ", " << mesh.vertices[i][2] << "\n";
    f << "*ELEMENT, TYPE=C3D4, ELSET=ALLELEMENTS\n";
    for (size_t i=0; i<mesh.tetrahedra.size(); ++i) {
        const auto& tet = mesh.tetrahedra[i];
        f << i+1 << ", " << tet[0]+1 << ", " << tet[1]+1 << ", " << tet[2]+1 << ", " << tet[3]+1 << "\n";
    }
    f << "*MATERIAL, NAME=BONE\n*ELASTIC\n" << mat.E << ", " << mat.nu << "\n";
    f << "*SOLID SECTION, ELSET=ALLELEMENTS, MATERIAL=BONE\n";
    if (fix_z_plane) {
        f << "*BOUNDARY\n";
        for (size_t i=0; i<mesh.vertices.size(); ++i)
            if (std::abs(mesh.vertices[i][2]-min_z) < eps)
                f << i+1 << ", 1, 3, 0.0\n";
    }
    if (load_z_plane && total_force>0) {
        std::vector<size_t> top;
        for (size_t i=0; i<mesh.vertices.size(); ++i)
            if (std::abs(mesh.vertices[i][2]-max_z) < eps) top.push_back(i);
        if (!top.empty()) {
            double per_node = total_force / top.size();
            f << "*CLOAD\n";
            for (size_t n : top) f << n+1 << ", 3, " << per_node << "\n";
        }
    }
    f << "*STEP\n*STATIC\n1.0, 1.0\n*END STEP\n";
    f.close();
}

void generate_run_script(const std::string& inp_path, const std::string& solver_path,
                         const std::string& output_dir) {
    fs::create_directories(output_dir);
    std::string script = output_dir + "/run_simulation.sh";
    std::ofstream f(script);
    if (!f) throw std::runtime_error("Cannot create " + script);
    fs::path inp(inp_path);
    std::string base = inp.stem().string();
    std::string dir = inp.parent_path().string();
    f << "#!/bin/bash\nexport OMP_NUM_THREADS=" << omp_get_max_threads() << "\n";
    f << "cd \"" << dir << "\"\n" << solver_path << " -i " << base << "\n";
    f.close();
    fs::permissions(script, fs::perms::owner_exec | fs::perms::group_exec | fs::perms::others_exec, fs::perm_options::add);
}

} // namespace