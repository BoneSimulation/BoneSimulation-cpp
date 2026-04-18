#pragma once
#include <string>
#include <vector>
#include <cstdint>

namespace bonesim {
    bool write_binary_to_tiff(const std::vector<uint8_t>& data,
                              size_t dim_z, size_t dim_y, size_t dim_x,
                              const std::string& filepath);
}