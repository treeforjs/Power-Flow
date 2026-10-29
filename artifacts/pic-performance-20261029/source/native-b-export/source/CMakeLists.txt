cmake_minimum_required(VERSION 3.24)
project(NativeMagneticFieldExport LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
find_package(MFEM 4.10 CONFIG REQUIRED)
add_executable(native_b_export native_b_export.cpp)
target_include_directories(native_b_export PRIVATE
    /home/treefor/projects/runtimes/mfem-gpu-pic-20261029/source/include
    /home/treefor/projects/runtimes/mfem-gpu-pic-20261029/source/src)
target_compile_options(native_b_export PRIVATE -ffunction-sections -fdata-sections)
target_link_options(native_b_export PRIVATE -Wl,--gc-sections)
target_link_libraries(native_b_export PRIVATE mfem)
