get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)

project(${PROJECT_NAME})
cmake_minimum_required(VERSION 3.19)

set(CMAKE_CXX_STANDARD 17)

include_directories(${PROJECT_SOURCE_DIR}/external/eigen)
include_directories(${PROJECT_SOURCE_DIR}/external/voronoi/src)
include_directories(${PROJECT_SOURCE_DIR}/external/stb)


add_executable(${PROJECT_NAME}  src/VoronoiApp.cpp src/ImageSampler.cpp)
configure_file(src/run.py ${CMAKE_BINARY_DIR}/run.py COPYONLY)