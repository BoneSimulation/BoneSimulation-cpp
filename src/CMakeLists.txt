cmake_minimum_required(VERSION 3.20)

find_package(OpenMP REQUIRED)
find_package(VTK REQUIRED COMPONENTS CommonCore CommonDataModel FiltersCore IOGeometry IOLegacy)
find_package(Eigen3 REQUIRED)
find_package(TIFF REQUIRED)
find_package(Boost REQUIRED COMPONENTS system thread)

include_directories(${VTK_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${CGAL_INCLUDE_DIRS})

add_executable(bonesim
        main.cpp
        core/image_loading.cpp
        core/image_processing.cpp
        core/cluster_analysis.cpp
        core/mesh_generation.cpp
        core/block_extraction.cpp
        core/calculix_export.cpp
        utils/file_io.cpp
        utils/logger.cpp
)

target_link_libraries(bonesim
        ${VTK_LIBRARIES}
        ${EIGEN3_LIBRARIES}
        ${TIFF_LIBRARIES}
        ${CGAL_LIBRARIES}
        ${Boost_LIBRARIES}
        OpenMP::OpenMP_CXX
)

set_target_properties(bonesim PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
)