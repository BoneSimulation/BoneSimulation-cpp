#include "reporting.hpp"
#include <hpdf.h>
#include <numeric>
#include <iostream>

namespace bonesim {

    Statistics compute_statistics(const BinaryVolume& vol) {
        auto res = find_largest_cluster(vol);
        Statistics stats;
        stats.num_clusters = res.num_clusters;
        stats.largest_cluster_size = res.largest_cluster_size;
        stats.mean_cluster_size = static_cast<double>(res.largest_cluster_size) / std::max(1ul, res.num_clusters);
        size_t total = vol.size();
        size_t bone = std::accumulate(vol.data.begin(), vol.data.end(), 0ul);
        stats.porosity = 1.0 - static_cast<double>(bone)/total;
        return stats;
    }

    void generate_pdf_report(const Statistics& stats, const std::string& output_path, const std::string& title) {
        HPDF_Doc pdf = HPDF_New(NULL, NULL);
        if (!pdf) { std::cerr << "PDF init failed\n"; return; }
        HPDF_Page page = HPDF_AddPage(pdf);
        HPDF_Page_SetSize(page, HPDF_PAGE_SIZE_A4, HPDF_PAGE_PORTRAIT);
        HPDF_Page_BeginText(page);
        HPDF_Page_SetFontAndSize(page, HPDF_GetFont(pdf, "Helvetica-Bold", NULL), 20);
        HPDF_Page_TextOut(page, 50, 800, title.c_str());
        HPDF_Page_EndText(page);
        HPDF_Page_BeginText(page);
        HPDF_Page_SetFontAndSize(page, HPDF_GetFont(pdf, "Helvetica", NULL), 12);
        double y = 750;
        auto line = [&](const std::string& s) {
            HPDF_Page_TextOut(page, 50, y, s.c_str()); y -= 20;
        };
        line("Number of clusters: " + std::to_string(stats.num_clusters));
        line("Porosity: " + std::to_string(stats.porosity*100) + " %");
        line("Largest cluster size: " + std::to_string(stats.largest_cluster_size) + " voxels");
        line("Mean cluster size: " + std::to_string(stats.mean_cluster_size) + " voxels");
        HPDF_Page_EndText(page);
        HPDF_SaveToFile(pdf, output_path.c_str());
        HPDF_Free(pdf);
    }

} // namespace