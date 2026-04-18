#include "image_loading.hpp"
#include <tiffio.h>
#include <filesystem>
#include <algorithm>
#include <cstring>
#include <iostream>

namespace fs = std::filesystem;
namespace bonesim {

static std::vector<std::string> get_tiff_files(const std::string& dir) {
    std::vector<std::string> files;
    for (const auto& entry : fs::directory_iterator(dir)) {
        auto ext = entry.path().extension();
        if (ext == ".tif" || ext == ".tiff")
            files.push_back(entry.path().string());
    }
    std::sort(files.begin(), files.end());
    return files;
}

static std::vector<uint16_t> load_single_tiff(const std::string& path, size_t& width, size_t& height) {
    TIFF* tif = TIFFOpen(path.c_str(), "r");
    if (!tif) throw std::runtime_error("Cannot open " + path);
    uint32_t w, h;
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
    width = w; height = h;
    std::vector<uint16_t> img(w * h);
    for (uint32_t row = 0; row < h; ++row)
        TIFFReadScanline(tif, &img[row * w], row, 0);
    TIFFClose(tif);
    return img;
}

Volume load_tiff_stack(const std::string& directory) {
    auto files = get_tiff_files(directory);
    if (files.empty()) throw std::runtime_error("No TIFF files in " + directory);
    size_t width, height;
    auto first = load_single_tiff(files[0], width, height);
    Volume vol;
    vol.dim_z = files.size();
    vol.dim_y = height;
    vol.dim_x = width;
    vol.data.resize(vol.dim_z * vol.dim_y * vol.dim_x);
    for (size_t z = 0; z < vol.dim_z; ++z) {
        size_t w, h;
        auto img = load_single_tiff(files[z], w, h);
        if (w != width || h != height) throw std::runtime_error("Inconsistent image dimensions");
        std::memcpy(&vol.data[z * width * height], img.data(), width * height * sizeof(uint16_t));
    }
    return vol;
}

Volume load_tiff_volume(const std::string& filepath) {
    TIFF* tif = TIFFOpen(filepath.c_str(), "r");
    if (!tif) throw std::runtime_error("Cannot open " + filepath);
    uint32_t width, height;
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
    size_t depth = 0;
    do { depth++; } while (TIFFReadDirectory(tif));
    TIFFClose(tif);

    Volume vol;
    vol.dim_x = width; vol.dim_y = height; vol.dim_z = depth;
    vol.data.resize(depth * height * width);
    tif = TIFFOpen(filepath.c_str(), "r");
    for (size_t z = 0; z < depth; ++z) {
        TIFFSetDirectory(tif, z);
        for (uint32_t y = 0; y < height; ++y) {
            TIFFReadScanline(tif, &vol.data[vol.index(z, y, 0)], y, 0);
        }
    }
    TIFFClose(tif);
    return vol;
}

} // namespace