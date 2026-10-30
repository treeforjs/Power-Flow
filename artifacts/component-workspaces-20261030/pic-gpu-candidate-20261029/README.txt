cuDSS SOLVE CUDA Graph candidate, 2026-10-29

State: CPU build complete. GPU capture qualification, checkpoint comparison,
and performance qualification remain pending. Default execution is the original
direct CUDA solve path. No active deployment or sealed source was changed.

Source was copied from:
  /home/treefor/projects/runtimes/mfem-gpu-pic-20261029/source
The initial hashes are in source-original-sha256.json.
Only src/cudss_poisson.hpp, CMakeLists.txt, and the new qualification program
differ. The Poisson residual guard, particle stepping, timestep rules,
checkpoint format, field assembly, and diagnostics source remain unchanged.

Build:
  taskset -c 0 nice -n 15 \
    /home/treefor/projects/runtimes/arch-base-20261028/run-frozen \
    /usr/bin/cmake -S PREFIX/source -B PREFIX/build \
    -DCMAKE_BUILD_TYPE=Release \
    -DMFEM_DIR=/home/treefor/projects/runtimes/mfem-4.10-pic-20261029/install-mpi/lib/cmake/mfem
  taskset -c 0 nice -n 15 \
    /home/treefor/projects/runtimes/arch-base-20261028/run-frozen \
    /usr/bin/cmake --build PREFIX/build --parallel 1
PREFIX is /home/treefor/projects/runtimes/pic-gpu-candidate-20261029.
All compilation uses the same frozen GNU compiler, Release flags, MFEM,
CUDA12.8 runtime, cuBLAS, and cuDSS0.6.0.5 as the sealed executable.

Documentation audit
The archived CUDA12.9.1 cuDSS release notes identify cuDSS0.6.0. Its general
description says SOLVE is asynchronous when hybrid memory and hybrid execution
are disabled. The source already disables both. SOLVE is also documented as
the only phase that allocates no device memory. The 0.6 description does not
contain an explicit CUDA Graph support promise; the CUDA13.1 cuDSS description
does, with restrictions concerning allocation and synchronous modes. Graph
capture on the pinned 0.6 library is therefore an experimental capability to
qualify, not an assumed supported production feature.

The 0.6 configuration documentation says IR_TOL is ignored and the configured
IR_N_STEPS number of steps is executed exactly. IR_N_STEPS remains 2, and no
tolerance, reordering, factorization, pivot, or numerical precision setting
changes. The unchanged outer guard requires:
  ||A*x-b||_2 <= max(1e-30, 2e-13*||b||_2)

The 0.6 documentation explicitly warns that atomic operations prevent a general
guarantee of bitwise reproducibility. No deterministic mode exists in the
pinned headers. This candidate enforces a stricter first-RHS bitwise comparison,
and requires external checkpoint comparisons before any adoption.

Implementation
RAIKOU_PIC_CUDSS_GRAPH selects off (default), required, or auto.

off preserves synchronous transfers and the original direct SOLVE call.
required creates a nonblocking stream and installs stream-ordered
cudaMallocAsync/cudaFreeAsync callbacks before analysis. Stable pinned host
RHS and answer buffers live as long as the solver. Analysis and factorization
remain outside capture. No later analysis or factorization occurs.

On the first actual RHS, an uncaptured unchanged SOLVE warms up the library
and produces a reference answer. Capture records pinned H2D, the full SOLVE
including both refinement iterations, and pinned D2H. Nested graphs are
inspected to require solve kernels and transfer nodes. A graph replay of that
same RHS must match the direct reference byte for byte before graph use is
marked qualified. Each later solve copies RHS into its stable host buffer,
launches the graph, synchronizes that stream, and copies the host answer into
the original MFEM answer vector. The outer residual guard still runs.

The first solve has one extra reference solve and capture/qualification cost.
Steady graph replay preserves the original solve and refinement work.
Graph and graph executable are destroyed before their cuDSS data, handle,
device buffers, host buffers, and stream.

required fails on unsupported capture or any qualification/runtime/numerical
error. auto permits fallback only for explicit cudaErrorStreamCaptureUnsupported,
cudaErrorNotSupported, CUDSS_STATUS_NOT_SUPPORTED, or the documented lack of
stream-ordered memory-pool support. cudaErrorStreamCaptureInvalidated alone is
fatal. Any unrelated error remains fatal even if a secondary unsupported
capture error is also returned. A first-RHS mismatch is always fatal.

Set RAIKOU_PIC_CUDSS_GRAPH_RECORD to an absolute JSON path to record qualified,
unsupported, failed, or fallback status. Successful capture records kernel,
copy, and graph node counts plus first_rhs_bitwise_qualified=true.

Qualification commands, for an isolated GPU window only
Set the same sealed CUDA library path as the baseline:
  CUDA=/home/treefor/projects/runtimes/mfem-cuda-benchmark-20261029/cuda-libs/nvidia
  LIBS="$CUDA/cu12/lib:$CUDA/cublas/lib:$CUDA/cuda_runtime/lib"

Small-system test:
  /home/treefor/projects/runtimes/arch-base-20261028/run-frozen \
    /usr/bin/env LD_LIBRARY_PATH="$LIBS" OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 \
    RAIKOU_PIC_CUDSS_GRAPH_RECORD=/absolute/output/small-graph.json \
    PREFIX/build/cudss_graph_qualification

This test requires graph mode, compares 32 varied RHSs with independently
factorized direct solves bitwise, and applies the unchanged residual bound.
It includes a zero RHS and varied scales. It has been compiled, not run.

Full checkpoint replay:
  PREFIX/bin/run-candidate required /absolute/output/full-graph.json PIC_ARGUMENT...
Use a copied checkpoint directory, the same arguments as the sealed baseline,
and a bounded physical end time. Compare final checkpoint particle/field
data and histories against the sealed baseline, then benchmark without
profiling. required is the acceptance mode; an auto fallback run cannot
demonstrate a graph speedup.

For Nsight tracing, place the explicit env and build/heating_history_pic_gpu
argv directly after the profile options in capture-cuda-replay. Avoid nesting
run-candidate or another frozen launcher as the target of a profiler session.

No GPU workload, production signal, driver change, package installation,
kernel setting change, or active source edit occurred during construction.

Official sources, with downloaded content hashes in docs/sources.json:
  https://docs.nvidia.com/cuda/archive/12.9.1/cudss/release_notes.html
  https://docs.nvidia.com/cuda/archive/12.9.1/cudss/general.html
  https://docs.nvidia.com/cuda/archive/12.9.1/cudss/functions.html
  https://docs.nvidia.com/cuda/archive/12.9.1/cudss/types.html
  https://docs.nvidia.com/cuda/archive/13.1.0/cudss/general.html
  https://docs.nvidia.com/cuda/archive/12.8.1/cuda-runtime-api/group__CUDART__STREAM.html
  https://docs.nvidia.com/cuda/archive/12.8.1/cuda-runtime-api/group__CUDART__GRAPH.html
