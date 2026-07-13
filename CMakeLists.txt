cmake_minimum_required(VERSION 3.26)
project(mhdlab_core LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(_mhd_core src/mhd_core/bindings.cpp src/mhd_core/reduced.cpp)
target_include_directories(_mhd_core PRIVATE src/mhd_core)
target_compile_definitions(_mhd_core PRIVATE VERSION_INFO="${PROJECT_VERSION}")
install(TARGETS _mhd_core LIBRARY DESTINATION mhdlab)
