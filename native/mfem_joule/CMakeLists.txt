cmake_minimum_required(VERSION 3.20)
project(power_flow_mfem_joule LANGUAGES CXX)

find_package(MFEM 4.10 CONFIG REQUIRED)
find_package(MPI REQUIRED COMPONENTS CXX)
if(NOT MFEM_USE_MPI)
  message(FATAL_ERROR "The installed MFEM must support MPI")
endif()
set(MFEM_SOURCE_ROOT "" CACHE PATH "Matching external MFEM source tree")
foreach(path miniapps/electromagnetics/joule.cpp miniapps/electromagnetics/joule_solver.cpp
             miniapps/common/pfem_extras.cpp)
  if(NOT EXISTS "${MFEM_SOURCE_ROOT}/${path}")
    message(FATAL_ERROR "Missing matching external source: ${path}")
  endif()
endforeach()

# Compile the unmodified upstream miniapp and helpers against installed MFEM.
# No MFEM sources are copied into this project or modified by the build.
set(COMMON_ROOT "${MFEM_SOURCE_ROOT}/miniapps/common")
add_library(power_flow_mfem_common STATIC
  "${COMMON_ROOT}/mesh_extras.cpp" "${COMMON_ROOT}/fem_extras.cpp"
  "${COMMON_ROOT}/particles_extras.cpp" "${COMMON_ROOT}/pfem_extras.cpp"
  "${COMMON_ROOT}/dist_solver.cpp")
target_link_libraries(power_flow_mfem_common PUBLIC mfem MPI::MPI_CXX)
target_compile_features(power_flow_mfem_common PUBLIC cxx_std_17)
add_executable(power-flow-joule
  "${MFEM_SOURCE_ROOT}/miniapps/electromagnetics/joule.cpp"
  "${MFEM_SOURCE_ROOT}/miniapps/electromagnetics/joule_solver.cpp")
target_link_libraries(power-flow-joule PRIVATE power_flow_mfem_common)
add_executable(power-flow-joule-inspect inspect_snapshot.cpp)
target_link_libraries(power-flow-joule-inspect PRIVATE mfem MPI::MPI_CXX)
target_compile_features(power-flow-joule-inspect PRIVATE cxx_std_17)
