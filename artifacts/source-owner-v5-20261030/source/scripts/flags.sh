# Shared include/define flags for direct compiles (same as the CMake targets).
V3=/home/treefor/projects/runtimes/pic-electrode-heat-feedback-20261029/native/mfem_electrode_heat_source_v3
V2=/home/treefor/projects/runtimes/pic-electrode-heat-feedback-20261029/native/mfem_electrode_heat_feedback_v2
M=/home/treefor/projects/runtimes/pic-electrode-heat-source-v5-20261030/native/mfem_electrode_heat_source_v5
MFEM=/home/treefor/projects/runtimes/mfem-4.10-pic-20261029
COUNT=/home/treefor/projects/runtimes/pic-water-exact-count-subflow-20261029/native/source/include/raikou/exact_particle_count.h
FLAGS="-std=gnu++20 -DHEAT_SOURCE_EXACT_COUNT_HEADER=\"$COUNT\" -I$V3/hall-owner/include -I$V3/hall-owner/src -I$M -I$V3 -I$V2 -I$V2/source-kernel/include -isystem $MFEM/install-mpi/include -isystem $MFEM/tpl/hypre/include -isystem $MFEM/tpl/metis/include -O2 -ffp-contract=off -fno-fast-math"
