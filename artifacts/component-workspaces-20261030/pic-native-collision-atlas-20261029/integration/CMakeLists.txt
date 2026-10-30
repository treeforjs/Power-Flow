cmake_minimum_required(VERSION 3.24)
project(IndexedCensusHelper LANGUAGES CXX)
find_package(MFEM 4.10 CONFIG REQUIRED)
find_package(OpenSSL REQUIRED)
add_executable(test_complete_indexed_census test_complete_indexed_census.cpp)
target_include_directories(test_complete_indexed_census PRIVATE . ../include
 /home/treefor/projects/runtimes/pic-water-native-subflow-20261029/source
 /home/treefor/projects/runtimes/pic-water-native-subflow-20261029/inputs
 /home/treefor/projects/runtimes/pic-water-native-subflow-20261029/inputs/include
 ../native-v6_1/source/include ../native-v6_1/inputs)
target_compile_options(test_complete_indexed_census PRIVATE -std=c++20 -O2 -ffp-contract=off -Wall -Wextra -Wpedantic)
target_link_libraries(test_complete_indexed_census PRIVATE
 /home/treefor/projects/runtimes/pic-native-collision-atlas-20261029/build/libphysical_collision_atlas.a mfem OpenSSL::Crypto)
add_library(physical_collision_atlas_fulltime STATIC ../source/physical_collision_atlas.cpp)
target_include_directories(physical_collision_atlas_fulltime PUBLIC ../include PRIVATE fulltime-inputs ../inputs)
target_compile_options(physical_collision_atlas_fulltime PRIVATE -std=c++20 -O2 -ffp-contract=off -Wall -Wextra -Wpedantic)
target_link_libraries(physical_collision_atlas_fulltime PUBLIC OpenSSL::Crypto)
add_executable(test_full_time_helper_api test_full_time_helper_api.cpp)
target_include_directories(test_full_time_helper_api PRIVATE . ../include
 /home/treefor/projects/runtimes/pic-water-full-time-subflow-20261029/source
 /home/treefor/projects/runtimes/pic-water-full-time-subflow-20261029/inputs
 /home/treefor/projects/runtimes/pic-water-full-time-subflow-20261029/inputs/include
 /home/treefor/projects/runtimes/pic-water-full-time-subflow-20261029/native-draft/include
 /home/treefor/projects/runtimes/pic-water-full-time-subflow-20261029/native-draft/inputs)
target_compile_options(test_full_time_helper_api PRIVATE -std=c++20 -O2 -ffp-contract=off -Wall -Wextra -Wpedantic)
target_link_libraries(test_full_time_helper_api PRIVATE physical_collision_atlas_fulltime mfem OpenSSL::Crypto)
