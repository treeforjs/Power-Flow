cmake_minimum_required(VERSION 3.24)
project(NativeCollisionAtlas LANGUAGES CXX)
find_package(MFEM 4.10 CONFIG REQUIRED)
find_package(OpenSSL REQUIRED)
add_library(physical_collision_atlas STATIC source/physical_collision_atlas.cpp)
target_include_directories(physical_collision_atlas PUBLIC include PRIVATE inputs)
target_compile_options(physical_collision_atlas PRIVATE -std=c++20 -O2 -ffp-contract=off -Wall -Wextra -Wpedantic)
target_link_libraries(physical_collision_atlas PUBLIC OpenSSL::Crypto)
add_executable(test_native_collision_atlas tests/test_native_collision_atlas.cpp)
target_include_directories(test_native_collision_atlas PRIVATE inputs native-v6_1/source/include native-v6_1/inputs native-v6_1/vendor)
target_compile_options(test_native_collision_atlas PRIVATE -std=c++20 -O2 -ffp-contract=off -Wall -Wextra -Wpedantic)
target_link_libraries(test_native_collision_atlas PRIVATE physical_collision_atlas mfem OpenSSL::Crypto)
add_executable(test_root_header_coexist tests/test_root_header_coexist.cpp)
target_include_directories(test_root_header_coexist PRIVATE
 /home/treefor/projects/runtimes/pic-water-native-subflow-20261029/source
 /home/treefor/projects/runtimes/pic-water-native-subflow-20261029/inputs
 /home/treefor/projects/runtimes/pic-water-native-subflow-20261029/inputs/include
 native-v6_1/source/include native-v6_1/inputs)
target_compile_options(test_root_header_coexist PRIVATE -std=c++20 -O2 -ffp-contract=off -Wall -Wextra -Wpedantic)
target_link_libraries(test_root_header_coexist PRIVATE physical_collision_atlas mfem OpenSSL::Crypto)
