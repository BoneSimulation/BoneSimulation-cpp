// image_processing.hpp
#pragma once

#include <cstdint>
#include <cstddef>
#include <vector>
#include <string>

// Vorwärtsdeklarationen
namespace bonesim {
    struct Volume;
    struct BinaryVolume;
    struct ClusterResult;
}

#include "image_loading.hpp"

namespace bonesim {

    // Struktur für binäres Volume
    struct BinaryVolume {
        std::vector<uint8_t> data;
        size_t dim_z, dim_y, dim_x;

        size_t index(size_t z, size_t y, size_t x) const {
            return (z * dim_y + y) * dim_x + x;
        }

        uint8_t& at(size_t z, size_t y, size_t x) {
            return data[index(z, y, x)];
        }

        const uint8_t& at(size_t z, size_t y, size_t x) const {
            return data[index(z, y, x)];
        }

        size_t size() const { return data.size(); }

        static BinaryVolume from_volume(const Volume& vol, uint16_t threshold);
        void save_as_tiff(const std::string& filepath) const;
    };

    // Funktionen
    uint16_t compute_otsu_threshold(const Volume& vol, size_t max_samples = 1000000);
    void binary_from_volume_inplace(Volume& vol, uint16_t threshold);
    BinaryVolume morphological_closing(const BinaryVolume& vol, int radius = 2);
    BinaryVolume morphological_closing_optimized(const BinaryVolume& vol, int radius);
    BinaryVolume interpolate_volume(const BinaryVolume& vol, double scaling_factor);
    BinaryVolume interpolate_volume_chunked(const BinaryVolume& vol, double scaling_factor, size_t chunk_size = 32);
    BinaryVolume crop_volume(const BinaryVolume& vol,
                             size_t crop_z_top, size_t crop_z_bottom,
                             size_t crop_y_front, size_t crop_y_back,
                             size_t crop_x_left, size_t crop_x_right);

} // namespace bonesim