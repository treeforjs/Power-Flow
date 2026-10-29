Native magnetic field export, 2026-10-29

Ready data:
  /home/treefor/Documents/Codex/2026-10-28/we/outputs/electron-mapping-paused-26p48ns/native-b/manifest.json

The manifest indexes 28 existing frame times, including the paused frame at
exact source time 26.478684721998123 ns. Each CSV contains 129600 unchanged
source snapshot cell centres, in y-major/x-fast order, and columns:
  time_s,time_ns,x_m,y_m,Bx_T,By_T,Bz_T,Bmag_T

Bx, By, Bz, and their Euclidean magnitude are absolute Tesla values.
No frame normalization or magnitude rescaling was performed. A renderer can
use the same fixed Tesla contour levels across all frames. The global |B|
range is 0.0193218443238617 to 58.419107320956606 T.

Sampling uses the exact sealed production source:
  /home/treefor/projects/runtimes/mfem-gpu-pic-20261029/source/examples/mfem_heating_history_pic.cpp
  /home/treefor/projects/runtimes/mfem-gpu-pic-20261029/source/examples/native_heating_fields.hpp

The utility renames and includes the original driver, then calls its original
Fields::value(x,y,t,true). The original driver entry is discarded by linker
garbage collection and is never called. Only frozen MFEM is linked; CUDA,
cuDSS, PIC driver, and PIC field-solve entries are absent.

The original sampler transforms PIC coordinates to native heating coordinates
with x_native=x_PIC-0.0045 m and y_native=y_PIC-0.0045 m. It evaluates the native
Az grid-function gradient through the original MFEM mesh transformations:
  B = (dAz/dy, -dAz/dx, 0).
The exact original Fields::set interval selection and clamped linear temporal
interpolation combine the original native endpoint vectors. No analytic or
inferred magnetic field is used.

The original Fields result was compared against independent direct original
NativeField endpoint evaluations and their midpoint interpolation at seven
grid points in the final frame's native interval. All 63 component comparisons
passed with maximum absolute difference 0 T. Every output component/magnitude
was checked finite while exporting. Per-frame and global min/max are in the
manifest. All sealed source, native field input, source frame, utility, and
output hashes are recorded; input hashes were checked unchanged after export.

Exact time_s source text is retained in time_s_decimal. time_ns_exact and the
by_time_ns keys are its exact decimal conversion. Numeric CSV/JSON time values
use binary double; the sampler used the unchanged source time_s double.

Build:
  taskset -c 6 nice -n 15 \
    flock /home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-profile/testing-core.lock \
    /bin/bash /home/treefor/projects/runtimes/pic-native-b-export-20261029/build.sh

Run:
  taskset -c 6 nice -n 15 \
    flock /home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-profile/testing-core.lock \
    /home/treefor/projects/runtimes/arch-base-20261028/run-frozen \
    /usr/bin/env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 \
    /home/treefor/projects/runtimes/pic-native-b-export-20261029/build/native_b_export \
    /home/treefor/projects/Power-Flow/runs/restore_field_mpi8_domain40mm \
    /home/treefor/projects/runtimes/pic-native-b-export-20261029/jobs.tsv \
    /home/treefor/Documents/Codex/2026-10-28/we/outputs/electron-mapping-paused-26p48ns/native-b

No GPU workload, simulation advance, production signal, cgroup modification,
renderer edit, or shared-source edit occurred.
