cmake_minimum_required(VERSION 3.24)
project(pic_electrode_heat_source_v5 LANGUAGES CXX)
# Owner v5 component, its focused tests and the audit reproducer port.
set(ELECTRODE_HEAT_SOURCE_V5_BUILD_TESTS ON CACHE BOOL "" FORCE)
add_subdirectory(native/mfem_electrode_heat_source_v5)
add_subdirectory(audit-port)
