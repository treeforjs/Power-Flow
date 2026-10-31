Resolution-series preparation -- 2026-10-31 (Power-Flow coupled PIC water breakdown)

Purpose: everything needed to run the mesh series (200/100/50/25 um) and the macroparticle series (1e4/1e5/1e6 carriers per
species) on the GPU-resident path as soon as the coupled code is integrated (pf-pilot-1 LIMITATIONS L6, pre-integration audit
A14), and early detection of what will not scale. Results of record: RESULTS.txt. Nothing sealed was modified: Hall corrective
v3.1, Hall 8dddec0, hall-gpu-resident-v1, native candidate003 and the C19 trace are read in place and hashed before and after
every build and run (build_run.py); native-full-time-candidate007 and the other listed worker directories were not read or touched.

Layout
  build_run.py              receipted build/run inside arch-base-20261028/run-frozen (sealed-input hashes, leases or SimQueue
                            context, binary hashes); outputs /home/treefor/projects/.workspace/outputs/resolution-series-prep-20261031/
  lent.sh                   one CPU-only job on a lent core (0, 1, 3) with its lease, nice 15; refuses after RETURN-CORES-0-3
  submit_gpu_series.sh      GPU runs as SimQueue jobs (simq submit -p powerflow -c 1 --gpus 1 --gpu-mem X), one job per run
  run_gpu_series.sh, run_gpu_dtscan.sh, run_gpu_series2.sh
                            the first GPU runners (any-core.sh + in-process GPU locks taken after the core lease: wrong lock
                            order, caused a deadlock at 10:37 and were withdrawn; kept as a record, do not use)
  graded/graded_mesh.hpp    graded transverse Yee mesh (uniform core over the 9 mm kinetic box, growth to a 35-45 mm box) with
                            the Hall operators; graded/run_graded_study.sh (sizes, DC port inductance versus box and grading)
  meshes/                   the series meshes in the form Hall admits (MeshSpec + identities, sizes, admission record):
                            mesh-200um.json, mesh-100um.json, mesh-50um.json, mesh-25um.json, mesh-100um-axial8.json,
                            mesh-50um-axial16.json, INDEX.json; SIZES.json (also variants too large to build: axial 8/16/32 at
                            25 um, isotropic z); make_mesh_files.py, size_table.py; graded-<h>um-box<35|45>mm-q<1p1|1p2>.json
                            (node coordinates, cell sizes, counts, memory, DC port inductance) and GRADED-INDEX.json;
                            make_graded_files.py
  hall-scaling/field_scaling.cpp
                            driver over the sealed Hall v3.1 public types (MeshSpec other than the default; no Hall source
                            changed): mesh-report, admit (instrumented copy of the DC seed solves, compared bitwise with the
                            Model), dc-floor (attainable accuracy of the DC seed solves), cpu-vacuum, gpu-vacuum, gpu-dtscan,
                            gpu-coupled (closed Hall step with the GPU push and 2..1e5 carriers), graded-dc, graded-sizes;
                            summarize_gpu.py
  precond/                  deterministic GPU geometric multigrid preconditioner for the Hall edge operator (standalone):
                            gmg.hpp (hierarchy, Galerkin operators, CPU reference), gmg_core.h (row arithmetic shared by host and
                            NVRTC), gmg_kernels.cu, gmg_gpu.hpp (hall_gpu::Session), test_gmg.cpp (gmg_solve)
  timestep/dt_limits.py     CFL/displacement, gyro, plasma, Debye and field-only limits per mesh (coordinator v7 constants)
  merge-split/              conservative macroparticle merge/split: CPU reference, NVRTC kernels, tests, INTEGRATION.txt (what
                            native v8 and Hall must add), README.txt, RESULTS.txt
  runmatrix/run_matrix.py   run matrix (mesh x carriers x dt) with memory and wall-time estimates; RUN-MATRIX.txt
  tests/check_results.py    acceptance checks over the recorded outputs (tests/TESTS-RESULT.json); tests/run_sanitizers.sh
                            (ASan + UBSan, build s1)
  finalize.sh               regenerates dt limits (timestep/DT-LIMITS.txt), run matrix (runmatrix/RUN-MATRIX.txt) and the checks

How to run (cwd under /home/treefor/projects; every compile and run goes through run-frozen)
  build:        ./lent.sh 3 /usr/bin/python -I build_run.py build LABEL [--sanitize]
                (or simq submit -p powerflow -c 1 --mem 6G -t 30 --borrow -C <this dir> -- /usr/bin/python -I build_run.py build LABEL)
  CPU run:      ./lent.sh K /usr/bin/python -I build_run.py run LABEL NAME field_scaling admit {OUT} {TRACE} 90 4 1
  GPU run:      simq submit -n NAME -p powerflow -c 1 --mem 8G --gpus 1 --gpu-mem 4G -t 12 -C <this dir> -- \
                    /usr/bin/python -I build_run.py run LABEL NAME field_scaling gpu-vacuum {OUT} {TRACE} {NVRTC} {UUID} 180 4 1 400 1e-14 1e-14 100
                or flock testing-gpu.lock flock testing-gpu-<UUID>.lock core-lease-20261030/any-core.sh --wait /usr/bin/python -I build_run.py run ...
                (GPU locks first, then the core lease; build_run.py refuses a GPU run unless the locks are held or SimQueue
                allocated the GPU). Every GPU job here is well under 10 minutes.
  dt limits:    timestep/dt_limits.py OUTDIR MESH_REPORT.json GPU_VACUUM.json... --dtscan GPU_DTSCAN.json...
  run matrix:   runmatrix/run_matrix.py OUTDIR runmatrix/inputs.json
  checks:       tests/check_results.py OUTPUT_ROOT RESULT.json ; tests/run_sanitizers.sh K
  merge/split:  merge-split/README.txt

Build labels of record: b1 (mesh report, first admissions and field-only runs), b3 (DC accuracy diagnostic, axial-variant
admissions), b5/b6 (GMG), b9 (proxy, coupled and smooth-GMG runs via SimQueue), b10/b11 (graded meshes, box study, uniform
large boxes), s1 (ASan+UBSan). b1-b8 take the GPU locks inside the program after the core lease (the withdrawn order); from b8 on
the programs take no lock and run under SimQueue or the GPU-locks-first wrappers. Every receipt records the binary sha256; the
flag own_sources_unchanged_since_build is false for runs whose sources were edited later (new modes added), which does not
affect the binary. Receipts: OUTPUT/LABEL/build/build-receipt.json and run-receipt.json next to every run.

Running the series once the coupled code is integrated (what the coordinator needs)
  0. Box: use the graded 35-45 mm meshes (meshes/graded-*.json, RESULTS F9 / 1b) once Hall accepts per-index coordinates and a
     kinetic box separate from the field box; until then the uniform 9 mm meshes are the reference (port inductance 11.0 nH
     instead of ~16.2-16.4 nH, gap voltage ~32 % low) and Hall v3.1 admits uniform 35/45 mm boxes only at 200 and 100 um.
  1. One mesh option that sets three things together (coordinator v7 driver_v7.cpp hard-codes the 200 um mesh in three places):
     the Hall MeshSpec (load_inputs: make_hall_model(..., hc::MeshSpec{}, ...)), the source-setup mesh (load_inputs:
     axial::Mesh g(45, 4, .009, 1) -> make_source_setup: one source slice per foil per axial cell, streams per body x slice x face,
     thermal slice cells 0.005/h x 0.0002/h) and the chemistry collision grid (request(): rq.grid_counts = {45, 45, 5}). Use
     meshes/mesh-*.json for the MeshSpec and check the run against its configuration_sha256 / motion_mesh_sha256 (admitted meshes)
     and kinetic_geometry_sha256. Whether the collision grid refines with the field mesh is a decision to record per member (it
     sets the MCC statistics per bin; a fixed grid keeps the chemistry statistics identical across the mesh series).
  2. Hall: v3.2 (absolute Gauss accumulated allowance) is required for any window past ~1 ns; 25 um additionally needs a seed
     solve that passes its own check (RESULTS section 3.2). Configure CarrierLimits for >= 3 x carriers per species (default cap
     1e5, default host budget 1 GiB = ~196k carriers) and the push engine capacity (<= 25e6).
  3. Weights: the streams' w0_electron / w0_neutral (identity_v7.hpp make_source_setup) set the birth weights; merge/split
     (merge-split/INTEGRATION.txt, native v8 + Hall) holds the census at the member's carrier count. A member's target count per
     species and the merge/split policy digest belong in its configuration record.
  4. dt: coordinator v7 controller unchanged (c_plasma 0.03, c_gyro 0.025, c_disp 0.1, particle cap 10 fs); the dt/2 member halves
     every factor and the cap. On 25 um the displacement criterion (8.3 fs for relativistic electrons) binds before the cap.
  5. Pre-run: Hall-only, per mesh, with the dt cap raised (field accuracy and CG iteration budget: RESULTS 3.3, 4, 5), then a fresh
     native genesis bound to that Hall state.
  6. Record the metrics of RESULTS section 8 and the execution envelope; run members via SimQueue (one core, --gpus 1, --gpu-mem from
     runmatrix device bytes + 20 %).
  7. Note for review on 25 um: the source settings use local_linear_vacuum_extent_m = 25e-6 (identity_v7.hpp), equal to one cell of
     the 25 um mesh.
