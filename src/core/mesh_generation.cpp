#include "mesh_generation.hpp"
#include <vtkSmartPointer.h>
#include <vtkMarchingCubes.h>
#include <vtkImageData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkTetra.h>
#include <vtkCleanPolyData.h>
#include <vtkSmoothPolyDataFilter.h>
#include <cstdlib>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <cstring>
#include <unistd.h>
#include <algorithm>
#include <iostream>

namespace fs = std::filesystem;

namespace bonesim {

static vtkSmartPointer<vtkImageData> binary_to_vtk(const BinaryVolume& vol) {
    auto img = vtkSmartPointer<vtkImageData>::New();
    img->SetDimensions(vol.dim_x, vol.dim_y, vol.dim_z);
    img->SetSpacing(1.0, 1.0, 1.0);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 1);
    auto* ptr = static_cast<unsigned char*>(img->GetScalarPointer());
    std::memcpy(ptr, vol.data.data(), vol.data.size());
    return img;
}

Mesh marching_cubes(const BinaryVolume& vol, double sx, double sy, double sz) {
    auto vtk_img = binary_to_vtk(vol);
    auto mc = vtkSmartPointer<vtkMarchingCubes>::New();
    mc->SetInputData(vtk_img);
    mc->SetValue(0, 0.5);
    mc->Update();

    auto poly = mc->GetOutput();

    // Optional: Mesh glätten und reinigen (verbessert tetgen-Qualität)
    auto cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
    cleaner->SetInputData(poly);
    cleaner->SetTolerance(1e-6);
    cleaner->Update();

    auto smoother = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
    smoother->SetInputConnection(cleaner->GetOutputPort());
    smoother->SetNumberOfIterations(20);
    smoother->SetRelaxationFactor(0.1);
    smoother->SetFeatureAngle(60.0);
    smoother->Update();

    auto smoothPoly = smoother->GetOutput();

    Mesh mesh;
    vtkPoints* pts = smoothPoly->GetPoints();
    mesh.vertices.resize(pts->GetNumberOfPoints());
    for (vtkIdType i = 0; i < pts->GetNumberOfPoints(); ++i) {
        double p[3];
        pts->GetPoint(i, p);
        mesh.vertices[i] = {p[0] * sx, p[1] * sy, p[2] * sz};
    }

    vtkCellArray* cells = smoothPoly->GetPolys();
    cells->InitTraversal();
    vtkIdType npts;
    const vtkIdType* ids;
    while (cells->GetNextCell(npts, ids)) {
        if (npts == 3) {
            mesh.faces.push_back({static_cast<size_t>(ids[0]),
                                  static_cast<size_t>(ids[1]),
                                  static_cast<size_t>(ids[2])});
        }
    }
    return mesh;
}

void save_mesh_as_vtk(const Mesh& mesh, const std::string& filename) {
    if (!mesh.tetrahedra.empty()) {
        auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        auto pts = vtkSmartPointer<vtkPoints>::New();
        for (const auto& v : mesh.vertices) pts->InsertNextPoint(v[0], v[1], v[2]);
        grid->SetPoints(pts);
        for (const auto& tet : mesh.tetrahedra) {
            auto tetra = vtkSmartPointer<vtkTetra>::New();
            tetra->GetPointIds()->SetId(0, tet[0]);
            tetra->GetPointIds()->SetId(1, tet[1]);
            tetra->GetPointIds()->SetId(2, tet[2]);
            tetra->GetPointIds()->SetId(3, tet[3]);
            grid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
        }
        auto writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
        writer->SetFileName(filename.c_str());
        writer->SetInputData(grid);
        writer->Write();
    } else if (!mesh.faces.empty()) {
        auto poly = vtkSmartPointer<vtkPolyData>::New();
        auto pts = vtkSmartPointer<vtkPoints>::New();
        for (const auto& v : mesh.vertices) pts->InsertNextPoint(v[0], v[1], v[2]);
        poly->SetPoints(pts);
        auto tris = vtkSmartPointer<vtkCellArray>::New();
        for (const auto& f : mesh.faces) {
            vtkIdType ids[3] = {static_cast<vtkIdType>(f[0]),
                                static_cast<vtkIdType>(f[1]),
                                static_cast<vtkIdType>(f[2])};
            tris->InsertNextCell(3, ids);
        }
        poly->SetPolys(tris);
        auto writer = vtkSmartPointer<vtkPolyDataWriter>::New();
        writer->SetFileName(filename.c_str());
        writer->SetInputData(poly);
        writer->Write();
    }
}

void save_surface_as_vtk(const Mesh& mesh, const std::string& filename) {
    save_mesh_as_vtk(mesh, filename);
}

BinaryVolume reverse_binary(const BinaryVolume& vol) {
    BinaryVolume res = vol;
    for (auto& v : res.data) v = 1 - v;
    return res;
}

// Tetraeder-Meshing mit tetgen (verbesserte Parameter + Mesh-Glättung)
Mesh generate_tetrahedral_mesh(const BinaryVolume& vol, double voxel_size,
                               bool use_downsampling, size_t max_volume_size) {
    Mesh result;
    // 1. Oberflächenmesh mit Glättung erzeugen
    BinaryVolume working = vol;
    double eff_voxel = voxel_size;
    if (use_downsampling && vol.size() > max_volume_size) {
        double factor = std::pow(static_cast<double>(max_volume_size) / vol.size(), 1.0 / 3.0);
        factor = std::min(factor, 0.5);
        working = interpolate_volume(vol, factor);
        eff_voxel = voxel_size / factor;
    }
    Mesh surface = marching_cubes(working, eff_voxel, eff_voxel, eff_voxel);
    if (surface.faces.empty()) return result;

    // Temporäre Dateien
    char temp_name[] = "/tmp/tetgen_XXXXXX";
    int fd = mkstemp(temp_name);
    if (fd == -1) return result;
    close(fd);
    std::string base(temp_name);
    std::string poly_file = base + ".poly";
    std::string node_file = base + ".node";
    std::string ele_file = base + ".ele";

    // 2. .poly Datei schreiben (verbesserte Qualität: Keine Löcher, Regionen)
    std::ofstream poly(poly_file);
    poly << "# Part 1 - node list\n";
    poly << surface.vertices.size() << " 3 0 0\n";
    for (size_t i = 0; i < surface.vertices.size(); ++i) {
        poly << i+1 << " " << surface.vertices[i][0] << " " << surface.vertices[i][1] << " " << surface.vertices[i][2] << "\n";
    }
    poly << "# Part 2 - facet list\n";
    poly << surface.faces.size() << " 1\n";
    for (size_t i = 0; i < surface.faces.size(); ++i) {
        poly << "1\n";
        poly << "3 " << surface.faces[i][0]+1 << " " << surface.faces[i][1]+1 << " " << surface.faces[i][2]+1 << "\n";
    }
    poly << "# Part 3 - hole list (none)\n0\n";
    poly << "# Part 4 - region list (none)\n0\n";
    poly.close();

    // 3. tetgen mit robusten Parametern: -pq1.2a0.1
    // p: tetrahedralisiert, q: Qualität (min. Winkel 1.2 Grad), a: maximale Volumenbeschränkung
    std::string cmd = "tetgen -pq1.2a0.1 " + poly_file;
    int ret = system(cmd.c_str());
    if (ret != 0) {
        std::cerr << "tetgen failed with code " << ret << std::endl;
        // Aufräumen
        std::remove(poly_file.c_str());
        return result;
    }

    std::string node_file_final = base + ".1.node";
    std::string ele_file_final = base + ".1.ele";

    std::ifstream node(node_file_final);
    if (!node) {
        std::cerr << "tetgen output file missing: " << node_file_final << std::endl;
        std::remove(poly_file.c_str());
        return result;
    }
    size_t num_nodes, dim, attr, boundary;
    node >> num_nodes >> dim >> attr >> boundary;
    result.vertices.resize(num_nodes);
    for (size_t i = 0; i < num_nodes; ++i) {
        size_t idx; double x, y, z;
        node >> idx >> x >> y >> z;
        result.vertices[idx-1] = {x, y, z};
    }
    node.close();

    std::ifstream ele(ele_file_final);
    if (!ele) {
        std::cerr << "tetgen output file missing: " << ele_file_final << std::endl;
        std::remove(poly_file.c_str());
        std::remove(node_file_final.c_str());
        return result;
    }
    size_t num_tets, nodes_per_tet, attrs;
    ele >> num_tets >> nodes_per_tet >> attrs;
    result.tetrahedra.resize(num_tets);
    for (size_t i = 0; i < num_tets; ++i) {
        size_t idx, v0, v1, v2, v3;
        ele >> idx >> v0 >> v1 >> v2 >> v3;
        result.tetrahedra[i] = {v0-1, v1-1, v2-1, v3-1};
    }
    ele.close();

    // 5. Aufräumen
    std::remove(poly_file.c_str());
    std::remove(node_file.c_str());
    std::remove(ele_file.c_str());
    std::remove((base + ".1.node").c_str());
    std::remove((base + ".1.ele").c_str());
    std::remove((base + ".1.face").c_str());

    return result;
}

} // namespace bonesim