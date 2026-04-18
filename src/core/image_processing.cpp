// image_processing.cpp
#include "image_processing.hpp"
#include "image_loading.hpp"
#include <algorithm>
#include <random>
#include <cstring>
#include <tiffio.h>
#include <queue>
#include <cmath>
#include <stdexcept>
#include <numeric>
#include <iostream>

namespace bonesim {

// Otsu threshold
uint16_t compute_otsu_threshold(const Volume& vol, size_t max_samples) {
    std::vector<uint16_t> samples;
    samples.reserve(std::min(max_samples, vol.size()));
    if (vol.size() <= max_samples) samples = vol.data;
    else {
        std::random_device rd; std::mt19937 gen(rd());
        std::uniform_int_distribution<size_t> dist(0, vol.size()-1);
        for (size_t i = 0; i < max_samples; ++i) samples.push_back(vol.data[dist(gen)]);
    }
    const int max_val = 65535;
    std::vector<size_t> hist(max_val+1,0);
    for (auto v : samples) hist[v]++;
    size_t total = samples.size();
    double sum = 0;
    for (int i=0; i<=max_val; ++i) sum += i * hist[i];
    double sum_b = 0; size_t w_b = 0; double var_max = 0; uint16_t thresh = 0;
    for (int i=0; i<=max_val; ++i) {
        w_b += hist[i];
        if (w_b == 0) continue;
        size_t w_f = total - w_b;
        if (w_f == 0) break;
        sum_b += i * hist[i];
        double m_b = sum_b / w_b;
        double m_f = (sum - sum_b) / w_f;
        double var = w_b * w_f * (m_b - m_f) * (m_b - m_f);
        if (var > var_max) { var_max = var; thresh = i; }
    }
    return thresh;
}

// In-place Binärkonvertierung
void binary_from_volume_inplace(Volume& vol, uint16_t threshold) {
    for (auto& v : vol.data) {
        v = (v > threshold) ? 1 : 0;
    }
}

// BinaryVolume aus Volume
BinaryVolume BinaryVolume::from_volume(const Volume& vol, uint16_t threshold) {
    BinaryVolume res;
    res.dim_z = vol.dim_z; res.dim_y = vol.dim_y; res.dim_x = vol.dim_x;
    res.data.resize(vol.size());
    for (size_t i=0; i<vol.size(); ++i) res.data[i] = (vol.data[i] > threshold) ? 1 : 0;
    return res;
}

// Morphologisches Closing (einfach, 3D)
BinaryVolume morphological_closing(const BinaryVolume& vol, int radius) {
    // Dilatation
    BinaryVolume dilated = vol;
    for (size_t z = 0; z < vol.dim_z; ++z) {
        for (size_t y = 0; y < vol.dim_y; ++y) {
            for (size_t x = 0; x < vol.dim_x; ++x) {
                if (vol.at(z,y,x) == 0) continue;
                for (int dz = -radius; dz <= radius; ++dz) {
                    for (int dy = -radius; dy <= radius; ++dy) {
                        for (int dx = -radius; dx <= radius; ++dx) {
                            if (dz*dz + dy*dy + dx*dx > radius*radius) continue;
                            int nz = z + dz, ny = y + dy, nx = x + dx;
                            if (nz >= 0 && nz < (int)vol.dim_z &&
                                ny >= 0 && ny < (int)vol.dim_y &&
                                nx >= 0 && nx < (int)vol.dim_x) {
                                dilated.at(nz,ny,nx) = 1;
                            }
                        }
                    }
                }
            }
        }
    }
    // Erosion
    BinaryVolume closed;
    closed.dim_z = vol.dim_z; closed.dim_y = vol.dim_y; closed.dim_x = vol.dim_x;
    closed.data.resize(vol.size(), 0);
    for (size_t z = 0; z < vol.dim_z; ++z) {
        for (size_t y = 0; y < vol.dim_y; ++y) {
            for (size_t x = 0; x < vol.dim_x; ++x) {
                bool inside = true;
                for (int dz = -radius; dz <= radius && inside; ++dz) {
                    for (int dy = -radius; dy <= radius && inside; ++dy) {
                        for (int dx = -radius; dx <= radius && inside; ++dx) {
                            if (dz*dz + dy*dy + dx*dx > radius*radius) continue;
                            int nz = z + dz, ny = y + dy, nx = x + dx;
                            if (nz < 0 || nz >= (int)vol.dim_z ||
                                ny < 0 || ny >= (int)vol.dim_y ||
                                nx < 0 || nx >= (int)vol.dim_x ||
                                dilated.at(nz,ny,nx) == 0) {
                                inside = false;
                                break;
                            }
                        }
                    }
                }
                if (inside) closed.at(z,y,x) = 1;
            }
        }
    }
    return closed;
}

// Optimiertes Closing (Slice-basiert, weniger Speicher)
BinaryVolume morphological_closing_optimized(const BinaryVolume& vol, int radius) {
    BinaryVolume result;
    result.dim_z = vol.dim_z;
    result.dim_y = vol.dim_y;
    result.dim_x = vol.dim_x;
    result.data.resize(vol.size(), 0);
    
    std::vector<uint8_t> slice(vol.dim_y * vol.dim_x);
    std::vector<uint8_t> temp(vol.dim_y * vol.dim_x, 0);
    
    for (size_t z = 0; z < vol.dim_z; ++z) {
        // Kopiere Slice
        for (size_t y = 0; y < vol.dim_y; ++y) {
            std::memcpy(&slice[y * vol.dim_x], &vol.data[vol.index(z, y, 0)], vol.dim_x);
        }
        
        // Dilatation
        for (size_t y = 0; y < vol.dim_y; ++y) {
            for (size_t x = 0; x < vol.dim_x; ++x) {
                if (slice[y * vol.dim_x + x]) {
                    for (int dy = -radius; dy <= radius; ++dy) {
                        for (int dx = -radius; dx <= radius; ++dx) {
                            if (dx*dx + dy*dy <= radius*radius) {
                                int ny = y + dy, nx = x + dx;
                                if (ny >= 0 && ny < (int)vol.dim_y && nx >= 0 && nx < (int)vol.dim_x) {
                                    temp[ny * vol.dim_x + nx] = 1;
                                }
                            }
                        }
                    }
                }
            }
        }
        
        // Erosion
        for (size_t y = 0; y < vol.dim_y; ++y) {
            for (size_t x = 0; x < vol.dim_x; ++x) {
                bool inside = true;
                for (int dy = -radius; dy <= radius && inside; ++dy) {
                    for (int dx = -radius; dx <= radius && inside; ++dx) {
                        if (dx*dx + dy*dy <= radius*radius) {
                            int ny = y + dy, nx = x + dx;
                            if (ny < 0 || ny >= (int)vol.dim_y || nx < 0 || nx >= (int)vol.dim_x ||
                                !temp[ny * vol.dim_x + nx]) {
                                inside = false;
                                break;
                            }
                        }
                    }
                }
                if (inside) {
                    result.data[result.index(z, y, x)] = 1;
                }
            }
        }
        
        std::fill(temp.begin(), temp.end(), 0);
        
        if (z % 50 == 0) {
            std::cout << "Closing slice " << z << "/" << vol.dim_z << std::endl;
        }
    }
    
    return result;
}

// Einfache Interpolation
BinaryVolume interpolate_volume(const BinaryVolume& vol, double scaling_factor) {
    size_t new_z = static_cast<size_t>(vol.dim_z * scaling_factor);
    size_t new_y = static_cast<size_t>(vol.dim_y * scaling_factor);
    size_t new_x = static_cast<size_t>(vol.dim_x * scaling_factor);
    BinaryVolume res;
    res.dim_z = new_z; res.dim_y = new_y; res.dim_x = new_x;
    res.data.resize(new_z * new_y * new_x, 0);
    for (size_t z = 0; z < new_z; ++z) {
        size_t src_z = std::min(static_cast<size_t>(z / scaling_factor), vol.dim_z - 1);
        for (size_t y = 0; y < new_y; ++y) {
            size_t src_y = std::min(static_cast<size_t>(y / scaling_factor), vol.dim_y - 1);
            size_t src_idx = (src_z * vol.dim_y + src_y) * vol.dim_x;
            size_t dst_idx = (z * new_y + y) * new_x;
            for (size_t x = 0; x < new_x; ++x) {
                size_t src_x = std::min(static_cast<size_t>(x / scaling_factor), vol.dim_x - 1);
                res.data[dst_idx + x] = vol.data[src_idx + src_x];
            }
        }
    }
    return res;
}

// Chunk-basierte Interpolation
BinaryVolume interpolate_volume_chunked(const BinaryVolume& vol, double scaling_factor, size_t chunk_size) {
    size_t new_z = static_cast<size_t>(vol.dim_z * scaling_factor);
    size_t new_y = static_cast<size_t>(vol.dim_y * scaling_factor);
    size_t new_x = static_cast<size_t>(vol.dim_x * scaling_factor);
    
    BinaryVolume res;
    res.dim_z = new_z; res.dim_y = new_y; res.dim_x = new_x;
    res.data.resize(new_z * new_y * new_x, 0);
    
    for (size_t z = 0; z < new_z; z += chunk_size) {
        size_t z_end = std::min(z + chunk_size, new_z);
        for (size_t zz = z; zz < z_end; ++zz) {
            size_t src_z = std::min(static_cast<size_t>(zz / scaling_factor), vol.dim_z - 1);
            for (size_t y = 0; y < new_y; ++y) {
                size_t src_y = std::min(static_cast<size_t>(y / scaling_factor), vol.dim_y - 1);
                size_t src_idx = (src_z * vol.dim_y + src_y) * vol.dim_x;
                size_t dst_idx = (zz * new_y + y) * new_x;
                for (size_t x = 0; x < new_x; ++x) {
                    size_t src_x = std::min(static_cast<size_t>(x / scaling_factor), vol.dim_x - 1);
                    res.data[dst_idx + x] = vol.data[src_idx + src_x];
                }
            }
        }
        std::cout << "Interpolated chunk " << z / chunk_size + 1 << "/" << (new_z + chunk_size - 1) / chunk_size << std::endl;
    }
    return res;
}

// Crop Volume
BinaryVolume crop_volume(const BinaryVolume& vol, 
                         size_t crop_z_top, size_t crop_z_bottom,
                         size_t crop_y_front, size_t crop_y_back,
                         size_t crop_x_left, size_t crop_x_right) {
    size_t new_z = vol.dim_z - crop_z_top - crop_z_bottom;
    size_t new_y = vol.dim_y - crop_y_front - crop_y_back;
    size_t new_x = vol.dim_x - crop_x_left - crop_x_right;
    
    BinaryVolume res;
    res.dim_z = new_z; res.dim_y = new_y; res.dim_x = new_x;
    res.data.resize(new_z * new_y * new_x);
    
    for (size_t z = 0; z < new_z; ++z) {
        for (size_t y = 0; y < new_y; ++y) {
            std::memcpy(&res.data[res.index(z, y, 0)],
                       &vol.data[vol.index(z + crop_z_top, y + crop_y_front, crop_x_left)],
                       new_x);
        }
    }
    return res;
}

void BinaryVolume::save_as_tiff(const std::string& filepath) const {
    TIFF* tif = TIFFOpen(filepath.c_str(), "w");
    if (!tif) throw std::runtime_error("Cannot create " + filepath);
    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, dim_x);
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH, dim_y);
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_LZW);
    for (size_t z = 0; z < dim_z; ++z) {
        if (z > 0) TIFFWriteDirectory(tif);
        for (size_t y = 0; y < dim_y; ++y) {
            TIFFWriteScanline(tif, const_cast<uint8_t*>(&data[index(z,y,0)]), y, 0);
        }
    }
    TIFFClose(tif);
}

} // namespace bonesim