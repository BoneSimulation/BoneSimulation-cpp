#include "file_io.hpp"
#include <tiffio.h>
#include <iostream>

namespace bonesim {
    bool write_binary_to_tiff(const std::vector<uint8_t>& data,
                              size_t dim_z, size_t dim_y, size_t dim_x,
                              const std::string& filepath) {
        TIFF* tif = TIFFOpen(filepath.c_str(), "w");
        if (!tif) return false;
        TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, dim_x);
        TIFFSetField(tif, TIFFTAG_IMAGELENGTH, dim_y);
        TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
        TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
        TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
        TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_LZW);
        for (size_t z = 0; z < dim_z; ++z) {
            if (z > 0) TIFFWriteDirectory(tif);
            for (size_t y = 0; y < dim_y; ++y) {
                size_t idx = (z * dim_y + y) * dim_x;
                TIFFWriteScanline(tif, const_cast<uint8_t*>(&data[idx]), y, 0);
            }
        }
        TIFFClose(tif);
        return true;
    }
}