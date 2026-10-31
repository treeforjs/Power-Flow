Hall corrective owner v3.1 -- 2026-10-31 (device-resident field step)

Entry points: source/root_api.hpp (namespace hall_corrective_v3; the v3 Owner, unchanged behaviour) and
source/resident_owner.hpp (ResidentOwner: the same protocol with the field state resident on the GPU).
Derived from Hall corrective v3 (runtimes/hall-corrective-v3-20261030, read in place, unchanged). Sealed inputs used
unchanged: Hall 8dddec0 axial/foil field headers, kinetic_boundary.{hpp,cpp}, heat-v2 bindings, hall-gpu-resident-v1
common/ and push/, native candidate003 (coordinator v5 link).

1. What moved to the device (source/gpu/resident_kernels.cu, resident_field.{hpp,cpp}, resident_abi.h)
  Two device slots (accepted, candidate) hold E, B words and the complete charge ledger (free net words, gross,
  conductor and exterior words, bank residue) as structure-of-arrays. One prepare runs on the device:
    k_operator    conduction S from sigma (cells ascending), a = (.5 dt)(S+Y), Jacobi diagonal (faces ascending,
                  pow(c,2) precomputed on the host)                                          [only when sigma or dt change]
    k_face_rhs, k_rhs   curl words of E, estimated B, H B words, ampere words (faces ascending), compensated rhs
    gpu1 CG       the hall-gpu-resident-v1 deterministic Jacobi CG kernels and control flow (cg_kernels.hpp: the
                  same kernel bytes as the gpu1 backend) on the device rhs; a = .5 dt (S+Y), diag from k_operator
    k_update_E/B  E words + increment, compensated midpoint E, curl(mid E), B words, compensated midpoint B
    k_currents    ampere(mid B), displacement/solid/exterior/metal/total/short currents, guard terms
    k_ledger      metal/exterior incidence words (edges ascending per node), conductor/exterior ledger words, bank
                  transfer (incl. the x87-exact residue step), KCL node sums, endpoint free charge copy
    k_cells       cell Joule power (12 edges in cell order), |mid B| at cell centres
    k_continuity, k_admit_nodes/edges/cells, k_gauss_increment, k_omega_p     guards and controller quantities
    k_reduce1/2   fixed two-stage reduction trees (double-word sums, exact maxima, first failing index)
    k_leaves, k_gather_blocks   SHA-256 of 2048-byte blocks of every state array; sampled blocks for the host
  Every kernel is compiled by NVRTC through the common session with --fmad=false --prec-div=true --prec-sqrt=true
  --ftz=false; no atomics; every sum has a fixed order.
  Host per step: sigma(T) (8 slice owners), the thermal owner (MFEM heat-v2 Couplers, unchanged), the independent
  residual (rhs and increment downloaded), the Merkle roots and a random leaf sample, the controller receipt.
  Transfers per vacuum step: sigma (97 KB) up; rhs + increment (2 x 333 KB), 200 metal-cell Joule powers, ~1400
  leaf digests (45 KB), 16 sampled blocks (32 KB) and the reduction results (1.5 KB) down.

2. Parity (what is bitwise, what is guard-verified)
  Bitwise equal to the Hall v3 host arithmetic for the same increment (verified every step, tests/test_resident.cpp):
    E, B, E_low, B_low, clock, revision, free net/gross, conductor and exterior words, bank residue; authoritative
    midpoint E and B; cell Joule power, therefore every heat slice (heat digest = the unchanged v3 record);
    rhs (so rhs norm, independent residual and CG iteration count are equal to the gpu1 backend); KCL and Ampere
    residual maxima (order-independent maxima of bitwise per-element values); the controller quantities omega_p dt,
    |B| max, omega_c dt (equal to their host form). Coupled steps: motion records (v3 digest, v2 JSON seal, census
    digest, trial-field digest), fixed-point iteration counts.
  Increment: device CG (gpu1 kernels) = the v3 gpu1 backend bit for bit; with ResidentOptions::host_cpu_pcg the
    sealed CPU PCG runs on the downloaded rhs and the step equals Hall v3 cpu-pcg bit for bit (parity mode).
  Guard-verified (computed on the device in double-word arithmetic with fixed trees; Hall v3 uses x87 long double
    and sequential sums): absolute Gauss defect/bound, divB/bound, endpoint continuity ratio, Gauss increment and its
    activity, energy sums (field change, plasma work, Joule, circuit, outer loss/supply), terminal voltage, region
    inflows, the increment-norm indicator. The energy sums differ from the host sequential sums by <= 0.0013 of the
    host's own rounding bound (n-1) eps sum|t| over 200 steps; every guard decision is the same on every step run.
  Independent FP64 residual: stays on the host with the sealed apply arithmetic (HostApply, preallocated, bitwise
    equal to the sealed lambda) and the sealed long-double norm, on the downloaded rhs and increment and a host
    conduction from the metal cells (bitwise equal to Mesh::conduction). A device-side bitwise residual is not
    possible at this cost: the sealed norm is a sequential x87 sum (41,676 dependent 80-bit steps).
  Refusal texts and order: the sealed texts; admission/continuity refusals are evaluated before the solve result is
    used; with several simultaneous failures of host-side preconditions the first reported refusal can differ.

3. Identities
  Field state: hall-v3.1-field-state-merkle-v1 = SHA-256 over api, configuration, block size and, per array (E, B,
  E_low, B_low, free_net_hi, free_net_lo, free_gross, conductor_hi/lo, exterior_hi/lo, residue when present): name,
  element count, leaf count, leaf digests; then time, revision, declarations, residue flag. Leaf = SHA-256(u64 array
  tag, u64 block index, u64 byte count, block bytes), 2048-byte blocks, hashed on the device. The host recomputes
  16 random leaves per digest (with the seed, 1,096 state leaves), every leaf with ResidentOptions::verify_all.
  Detection: a single corrupted block is found with probability 16/1096 = 1.5 % per digest; a persistent corruption
  is also caught by the preflight recomputation of all leaves of both slots on the device (compared with the digested
  leaves) and by the next admissions. Midpoint role seal: hall-v3.1-authoritative-midpoint-merkle-v1 (same scheme).
  Configuration: backend gpu1-resident-v1 with a descriptor binding the rule text, kernel/ABI/CG source digests,
  NVRTC options, target architecture and launch geometry. Heat, kinetic-ledger, conductivity, motion mesh, material
  and every motion/native-facing digest are the v3 records unchanged.

4. Other changes (CPU path unchanged in behaviour and identities)
  - Field-step upper bound as configuration: Model(..., ThermalDeclaration, double dt_max = 1e-14) is used by the
    field interval guard, the closed proposal clock guard, the thermal guard and the motion admission; it enters the
    configuration identity only when it differs from 1e-14, so default identities equal Hall v3.
  - Controller receipt (control_receipt.hpp): omega_p dt (free gross charge, electron mass), |B| max and omega_c dt
    (authoritative midpoint), drive current change and C19 knot spacing, increment norm, third-difference
    (trapezoidal LTE) estimate; filled on the device every resident step, host form host_control_receipt(); the
    particle quantities are the motion budget's maximum_gyro and maximum_displacement_cells.
  - Motion: MotionContext::defer_record_v3 (resident owner) computes the binary v3 record digest only for the accepted
    iterate (MotionAssembly::seal_v3); census ordering and the deposit-support/charge sorts skip the sort when the
    census is already in ascending ID order (same result).
  - Owner preflight: the accepted state's admission and field/heat digest recomputation run once per owner
    generation; the published candidate (admitted and re-hashed by its own final preflight) starts the new generation
    verified, so a step's preflight no longer re-hashes the unchanged accepted state (same values; a mismatch still
    refuses). heat_for no longer repeats the candidate admission that Model::prepare performs as its last statement.
  - Total-field (midpoint) digest computed once per Midpoint object (Midpoint::total_field_sha256, same record);
    prepare_closed_proposal passes the digest of the exact midpoint objects to the provider for one motion() call
    (VerifiedField, one-shot), and motion admission skips its recomputation only for those objects and that digest.
  - step_loads is generic over the owner type; Digest::bytes (raw leaf bytes); CurrentDrive::knot_interval.

5. Layout
  source/                      v3 sources (copied) with the changes above; resident_owner.hpp; control_receipt.hpp
  source/gpu/resident_*.{h,cu,hpp,cpp}, cg_kernels.hpp   device-resident step
  tests/test_resident.cpp      units, parity-vacuum, parity-coupled, negative, vacuum-run, coupled-run, gauss-trace,
                               dt-study (the last two are diagnostics)
  tests/profile_v3_baseline.cpp, profile_prepare.hpp   phase profile of Hall v3 (sealed sources) incl. an instrumented
                               copy of Model::prepare checked word for word against the candidate of record
  tests/compare_v3_evidence.py comparison of the v3 component outputs with the Hall v3 evidence of record
  run_hall_v31.py              build | run | gpu (receipted, frozen launcher, leases, sealed-input hashes)
  run_coordinator_parity.py    coordinator v5 T1-T7 on v3.1 vs the Hall v3 run (run-p3-221714)
  run_acceptance.sh            acceptance batch

6. Results of record: /home/treefor/projects/.workspace/outputs/hall-corrective-v3.1-20261031/RESULTS.txt (final build a5)

7. Open
  - The absolute Gauss allowance is fixed at construction while the per-step admitted Gauss increments accumulate:
    the vacuum step is refused at 1.236 ns (step 48,650 at 10 fs; 1.233 ns at 100 fs, step 4,833), independent of
    dt. Analysis and a proposed per-node accumulated allowance: RESULTS.txt section 6. No guard was changed.
  - Per-carrier host work is unchanged in kind: the v2 JSON seal, the census digest and the binary v3 record are
    verified by the sealed native (x87 kinetic sum compared for equality, streaming digests), so exact sums and
    Merkle motion digests need a native version that accepts them (RESULTS.txt section 5).
  - DepositSupportBuilder runs 4-6x slower in the resident test process than in the baseline profiler process with
    the same code and data (2.8 vs 0.64 ms at 1e4 carriers; repeat calls equal, so not a warm-up effect); not
    explained yet.
  - Checkpoint save/restore of the resident owner (accepted_state() downloads the state; no restore path yet).
  - Cross-architecture bit identity (sm_89, sm_120) not shown.
