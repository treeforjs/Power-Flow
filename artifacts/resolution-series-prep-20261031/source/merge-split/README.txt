merge-split -- conservative macroparticle merge/split transaction (2026-10-31)
Schema merge-split-transaction-v1. Pilot plan item H7 (weights controllable without breaking the exact
ledgers, for the 1e4/1e5/1e6 macroparticle series) on the CPU reference and the GPU-resident path (5c).
Integration into native v8 and Hall v4: INTEGRATION.txt. Measured numbers: RESULTS.txt.

1. Layout
  include/ms_core.h        shared host/NVRTC core: classification (Hall cell, native kinetic energy formula,
                           velocity bin), TwoSum/TwoProduct, exact superaccumulator (units 2^-1074), group
                           formation per (cell, species) segment, merge kinematics, split check, Philox
                           stream layout (core: pic-arrival-owner-v1 philox_core.h, included unchanged)
  include/merge_split.hpp  CPU reference API: Carriers, Exact, Result (groups, outputs, removed IDs, relocation
  src/merge_split.cpp      records, per-species exact ledger, counters), transact(), canonical encoding and
                           SHA-256 digest, partitioned execution (rank-independence tests), policy identity
  include/relocation.hpp   relocation current (port of the Hall axial::Path deposit), Words free charge and
  src/relocation.cpp       the Hall endpoint continuity defect/scale, relocation work, test fields
  gpu/ms_kernels.cu        NVRTC kernels (no atomics; fixed-order sums; exact integer ledger reduction)
  gpu/ms_gpu.hpp/.cpp      host launcher on the hall-gpu-resident-v1 common session and deterministic
                           primitives (stable LSD radix sort, exclusive scan, segment bounds)
  tests/fixtures.hpp       electron / H2O+ carrier sets (Maxwellian + keV tail, 9 x 9 x 25 mm box), policy
  tests/test_cpu.cpp       unit, kinematics, split-only, ledger closure (1e6 per species, 45/90/180/360
                           meshes), relocation continuity/work, partitions 1/2/5/8 contiguous and shuffled
  tests/test_gpu.cpp       CPU/GPU bit identity (every array, every ledger term, digest), repeat runs, timing
  ms_build.py              receipted build/run (run-frozen; sealed-input hashes before/after)

2. Build and run (core leases and GPU locks are taken by the caller)
  build:  flock .../pic-profile/testing-core3.lock taskset -c 3 nice -n 15 /usr/bin/python -I ms_build.py build LABEL
  CPU:    flock .../testing-core3.lock taskset -c 3 nice -n 15 /usr/bin/python -I ms_build.py run LABEL cpu-full \
              test_cpu --n 1000000 --out '{OUT}/results-cpu.json'
  ASan:   ... ms_build.py run LABEL cpu-asan test_cpu_asan --n 1000000 --out '{OUT}/results-cpu-asan.json'
  GPU:    flock .../testing-gpu.lock flock .../testing-gpu-<UUID>.lock core-lease-20261030/any-core.sh --wait \
              timeout 600 /usr/bin/python -I ms_build.py run LABEL gpu-full test_gpu --nvrtc '{NVRTC}' \
              --uuid '{UUID}' --n 1000000 --repeat 3 --out '{OUT}/results-gpu.json'
  Compiler: /usr/bin/c++ (GCC 16.2.1) -std=c++20 -O2 -ffp-contract=off -fno-fast-math; ASan+UBSan target
  -O1 -g -fsanitize=address,undefined -fno-sanitize-recover=all. NVRTC 12.8 (dlopen'ed, pinned path) with
  --fmad=false --prec-div=true --prec-sqrt=true --ftz=false, target taken from the device (sm_86 here).
  Outputs: /home/treefor/projects/.workspace/outputs/resolution-series-prep-20261031/merge-split/LABEL/.

3. Design
  Grouping. Census carriers are classified (Hall mesh cell of x + shift; K = pic_kinetic_energy_J; bin =
   optional spatial sub-bin of the cell (sub[3] divisions per axis, 1 = none) x energy bin (per-species
   edges) x cube-sphere direction cell) and put in canonical order (cell, species, bin, native ID). Each (cell, species) segment is processed sequentially:
   merge when N > merge_above (reduce toward merge_target): in each bin, a group starts at a merge-eligible
   carrier and takes the next eligible carriers of the bin while k < min(k_max, excess + 2) and the running
   weight sum stays exactly representable (TwoSum error zero); k >= k_min forms a group, otherwise the start
   carrier is left. Split when 0 < N < split_below (grow toward split_target): eligible carriers with
   w >= split_min_weight become 2^m clones with weight w 2^-m (formed only when w 2^-m is normal).
  Merge kinematics (k -> 2, same species): c = sum w x / W (double-word sums); if c is not strictly inside
   the group's cell, the first input position is used. U = P/W, kappa = E/W from double-word sums of w u and
   w K; |u_out|^2 = (kappa/m)(2 + kappa/(m c^2)) (so K(u_out) = kappa in real arithmetic); q^2 = |u_out|^2 -
   |U|^2 (negative: skipped, no real solution); u_a,b = U +/- q n with n perpendicular to U, azimuth from
   Philox via the rational map (cos, sin) = ((x^2 - y^2)/s, 2xy/s) of a rejection-sampled disk point;
   U == 0: isotropic n (Marsaglia). Only + - * / sqrt: host and device bits agree.
  Split: clones at the parent position with the parent velocity: exact in count, charge, momentum, energy.
   (Velocity-perturbed clones cannot conserve both momentum and energy: K is strictly convex, so at fixed
   sum w_i u_i the energy sum w_i K(u_i) is minimal, equal to w K(u), only when all u_i = u.)
  Exactness. Count: every merge group's weight sum is exactly representable by construction, outputs W/2;
   split weights are exact powers-of-two fractions. Ledger: every term (w, w m u_d as 4 binary64 words,
   w K as 2 words; Dekker TwoProduct with an exponent-range guard) goes into a 68-digit fixed-point
   superaccumulator (units 2^-1074); the totals are exact integers, independent of order and partitioning.
   Booked: merge_momentum_residual[3] and merge_energy_residual per species (exact signed integer, its
   round-to-nearest binary64 and canonical binary64 expansion); in + booked == out exactly.
  Relocation and Gauss. Merge inputs move (zero time, at the step begin) from x_i to c inside their cell;
   each move is a relocation record (start x + shift, D = (c + shift) - (x + shift), charge w q) deposited with
   the Hall path-deposit arithmetic into J_reloc, which joins the motion current in the same field step,
   so the endpoint continuity guard and absolute Gauss close; the field work on J_reloc is relocation_work_J.
   With c the weighted centroid the group dipole is unchanged: the uniform-field work cancels to rounding;
   field gradients inside a cell leave a second-order work that grows with the relocation length, which the
   spatial sub-bins bound (the 4 axial cells are 6.25 mm long; RESULTS.txt section 4).
  IDs. Outputs get M_s + rank (rank = exclusive scan of per-group output counts in canonical group order),
   like births (B_s of the step becomes M_s + n_out).
  Determinism. The result is a pure function of the carrier set: canonical order is global (partitions
   are sorted locally and k-way merged), groups are pure functions of their inputs, IDs come from a scan
   in canonical order, ledgers are exact. The GPU path computes the same functions (one thread per segment
   and per group; stable radix sorts; integer-exact ledger reduction) and is compared word for word.

4. Read-only inputs
  pic-arrival-owner-v1-20261030/include/pf_arrival/philox_core.h (Philox4x32-10 core);
  native-full-time-candidate006/source/include/raikou/exact_particle_count.h (tests only);
  hall-gpu-resident-v1-20261030/common (session, NVRTC, primitives); NVRTC 12.8 and CUDA headers at the
  pinned paths. Arithmetic ported (not linked): native pic_kinetic_energy_J; Hall axial::Path deposit,
  point/charge deposit, Words, incidence, continuity scale. Nothing outside merge-split/ was modified.
