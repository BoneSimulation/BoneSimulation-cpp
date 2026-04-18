#pragma once
#include <string>
#include <vector>
#include <cstdint>

namespace bonesim {
    struct Volume {
        std::vector<uint16_t> data;
        size_t dim_z, dim_y, dim_x;
        size_t index(size_t z, size_t y, size_t x) const { return (z * dim_y + y) * dim_x + x; }
        uint16_t& at(size_t z, size_t y, size_t x) { return data[index(z,y,x)]; }
        const uint16_t& at(size_t z, size_t y, size_t x) const { return data[index(z,y,x)]; }
        size_t size() const { return data.size(); }
    };

    Volume load_tiff_stack(const std::string& directory);
    Volume load_tiff_volume(const std::string& filepath);
}