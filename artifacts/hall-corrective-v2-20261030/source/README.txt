Hall corrective owner v2 — 2026-10-30

Entry point: source/root_api.hpp, namespace hall_corrective_v2 (GPU backends: source/gpu/gpu_backends.hpp,
namespace hall_corrective_v2::gpu). Derived from the Hall corrective owner at commit 8dddec0
(development/hall-field-correction-20261030, development/hall-canonical-driver-20261029); nothing in the Hall
repository, in hall-gpu-resident-v1-20261030 or in native candidate003 was modified or branched. Sealed inputs
used unchanged: 8dddec0 axial/foil field headers, the v2 heat owner (axial_thermal_v2.hpp and heat-v2 bindings),
kinetic_boundary.{hpp,cpp}; hall-gpu-resident-v1 common/ and push/ (compiled as is).

Changes relative to 8dddec0

1. Binary canonical identities (encoding.hpp). Same convention as native v5: u64 little-endian, f64 as the
   binary64 bit pattern, bool one byte, text = u64 length + bytes, sequences = u64 count + elements in index
   order, every record starts with a version tag. SHA-256 is streamed through a fixed 16 KiB buffer (no JSON value,
   no growing byte string). Records: hall-v2-field-configuration-v1 (api, complete mesh incl. all operator arrays
   and masks, unchanged material binding digest, C19 digest and roles, initial sigma/E/B/potential/feed/
   background/short current, constraint parameters, the endpoint-continuity rule text, the independent-residual
   rule text, field-solve backend name + descriptor digest), hall-v2-field-state-v1 (api, configuration digest,
   E, B, E_low, B_low, time, revision, charge ledger: free net words, gross, conductor and exterior words and
   declarations), authoritative-compensated-midpoint-v2, hall-v2-conductivity-v1, hall-v2-heat-state-v1 (unchanged
   material digest + per slice T, inventory words, watermarks, cumulative, time, revision).
   Hashed once per immutable object: the configuration at Model construction; a field candidate only when it is
   accepted by the fixed point (Model::seal, discarded iterates are not hashed); the owner caches its accepted
   field and heat digests and replaces them at publication by nothrow string swaps. Preflight still performs one
   independent recomputation of the accepted and candidate field/heat digests. Unchanged values are kept where the
   encoding did not change: motion mesh identity (now computed from the mesh JSON text cached at Model
   construction), material binding digest, census/total-field/kinetic-policy/motion record digests.
   Measured on the default mesh (41,676 edges, 39,015 faces, 14,812 nodes): field-state digest 0.95 ms binary vs
   8.0 ms for the 8dddec0 JSON text of the same state (8.4x); per-step identity work in the Stage 1 step fell from
   ~117 ms to a few ms (see the v5 receipts).
   JSON destructor terminations: no nonempty nlohmann value is destroyed on a v2 normal path. The only JSON built
   per step in v2 is the sealed motion record required by native admission; it is built by assignments and moves
   (no initializer-list temporaries) and drained before destruction.

2. Node-plane defect (Stage 1 defect 1). The endpoint charge/current guard is now
     |after - before - div(J) dt|_n <= 1e-12 (gross_before + gross_after + |div(J) dt|)_n + A_n,
     A_n = eps [ sum_{p in S(n)} Q_p (2C + 14 + (16C + 67) R_p) + 24 M_n sum_{p in S(n)} Q_p R_p ],
   eps = 2^-52, C = nx+ny+nz, Q_p = |w q| (double-word), R_p = sum_d |dx_d|/h_d, S(n) the carriers whose start/
   end/path cell range (widened by eps (3 c_d + 8) cells) contains node n, M_n = |S(n)|. In exact arithmetic the
   CIC point deposits and the cut-at-planes, two-point-Gauss path current conserve charge per node exactly; A_n
   bounds the rounding (cell-coordinate rounding of the start and end deposits, macro-charge and product
   roundings, quadrature coordinates, cut parameters, endpoint representation, recursive edge summation); the
   derivation is in census_charge.hpp. A_n is built only from the actual post-first-exposure census and the actual
   motion endpoints (DepositSupportBuilder; no caller-supplied tolerance) and is part of the configuration rule.
   Evidence (runs b4, b5): the Stage 1 reproducer (tests/repro_hall_node_plane_v2.cpp, cases unchanged) now admits every
   case incl. z = 0 at weights 1 and 1e9 (worst defect/tolerance 7.8e-4); random pairs on x, y and z node planes
   (offsets 0..1e-7 m) are admitted; with the absolute term removed the 8dddec0 refusal reproduces; 1 e of invented
   charge, 4x the absolute term invented, and a 1 nm endpoint shift without current are refused. The audit A09/A10
   counterexamples (artifacts/audit-reproducer-ports-20261030/hall-field) are still rejected: -e x 1e9 inconsistent
   charge at admit, prepare and bootstrap; 0.01 T divB face at admit, prepare, bootstrap and resealed restore, down
   to 1e-13 T; the compensated midpoint stays authoritative and closed proposals bind motion to its digest.

3. Field-solve backend interface (solver_backend.hpp). FieldSolveBackend {name, descriptor_sha256, solve}; the
   Model is constructed with a backend (default cpu-pcg = the sealed foil::solve) and binds its descriptor into the
   configuration identity, so a GPU run has a different configuration and field-state identity. After every solve
   Model::prepare performs the mandatory independent FP64 residual check with the sealed matrix-free apply and the
   sealed long-double norm, ||b - Ax|| <= max(1e-30, 2e-13 ||b||), for every backend including the CPU one; a 0.1 %
   increment error and a 1e-6 one-component error are refused by it (tests). Model::prepare is private (Owner and
   the closed proposal only); the owner supplies its cached accepted digest.
   gpu1-factored-cg (gpu/gpu_backends.cpp): the deterministic FP64 Jacobi CG of hall-gpu-resident-v1/field
   (factored operator, row-per-thread, fixed tree 512 x 256, no atomics, same stopping rule) running on the
   common hall_gpu::Session; its descriptor binds the algorithm text, kernel source digest, NVRTC options, target
   architecture and NVRTC version.

4. Motion providers (owner.hpp MotionProvider, motion_receipt.hpp, motion_provider.cpp). motion_receipt.hpp is
   the v2 equivalent of the sealed header: identical hall_driver types, layout and accessors, with
   hall_corrective_v2::MotionAssembly as an additional friend; motion_provider.cpp replaces the sealed
   native_arena_motion.cpp in the link (one definition per program) with the identical CPU arithmetic, refusal
   order/texts, digests and record seal, plus MotionAssembly::assemble for device arrays. Native candidate003 and
   native v5 adapters only forward-declare hall_driver::MotionReceipt and call hall_driver::motion_record, so both
   admit the v2 receipt unchanged. CpuMotionProvider = sealed arithmetic; GpuMotionProvider = hall-gpu-resident-v1
   push MotionEngine (carriers uploaded once per proposal, midpoint per iterate) with the host Count2 packet check
   placed in the sealed refusal order, endpoint free charge taken from the device (bitwise equal to
   charge_of_endpoints). The field CG uses the same session: one device session, one ordered non-blocking stream.
   Evidence (runs b4, b5): v2 CPU records equal the sealed 8dddec0 provider on 12 fixtures (record bytes, seal, endpoint
   charge words, refusal texts; 2 to 1e4 carriers, node planes, neutral carrier, carrier in metal, divergent packet);
   GPU records equal both on every fixture; GPU push + CPU field reproduces the CPU candidate field/heat digests,
   midpoint and iteration count; GPU push + GPU field at z = 0 is accepted (independent residual 4.6e-4 of target,
   Gauss, divB, energy) and bit-identical over three in-process runs and two fresh processes; teardown drained,
   not quarantined.
   Device-side right-hand-side assembly was not done: the RHS is one part of the compensated-word host arithmetic
   of Model::prepare (curl/ampere words, update, midpoint, Gauss increment, admission), which together take ~15 ms
   per iterate on the default mesh and stay on the host; moving the RHS alone would not make the step resident.

5. Motion record v3 (binary) and configurable Hall carrier cap (requested after native v5 was released).
   Every fixed-point iterate now produces only the binary "motion v3" digest of its record
   (record_v3_sha256(), tag canonical-full-time-total-field-motion-v3-binary). The v2 JSON-text seal is computed
   once per step, for the accepted motion only (MotionAssembly::seal_v2, called by prepare_closed_proposal at the
   fixed point), because native candidate003 and native v5 both still recognize only the v2 record and verify its
   JSON-text seal; motion_record() refuses an unsealed receipt. Owner::preflight recomputes the v3 digest.
   v2 record -> v3 record mapping (same content, canonical binary encoding, fields in the v2 key order):
     schema text                 -> tag "canonical-full-time-total-field-motion-v3-binary" + text of the v2 schema
     binding (14 keys)           -> text "binding", phase text, accepted_step, staged_arena_revision, trial_stage (u64),
                                    begin_s, end_s, dt_s (f64), staged_arena, trial_field, mesh, material, species,
                                    common_clock, kinetic_policy digests (text)
     budget (7 numbers)          -> text "budget", 7 f64 in the order kinetic_change, field_work, work_defect,
                                    continuity, plasma_dt, maximum_gyro, maximum_displacement_cells
     census_sha256, coverage, current_A_m_sha256, endpoint_charge_C_sha256 -> text key + text value
     endpoints[]                 -> text "endpoints", u64 count, per carrier u64 id, 3 f64 position, 3 f64 proper
                                    velocity, f64 charge, mass, weight (ascending native ID)
     initial_charge_C_sha256     -> text key + text value
     kinetic_policy              -> text "kinetic_policy", policy schema text, 12 f64 (width, gap, short_z, feed_z,
                                    lower foil l/r/b/t, upper foil l/r/b/t), 4 policy texts
     native_DKD_called, native_Poisson_called -> text key + u64 0
     native_to_field_m           -> text key + 3 f64
     neutral_motion_owner        -> text key + text value
     paths[]                     -> text "paths", u64 count, per carrier u64 id, 3 f64 start, 3 f64 displacement,
                                    f64 flight; then the motion owner text once
     seal_sha256 (v2)            -> not part of v3 (v3 is the record digest itself)
   The array digests (motion-array-binary64-v1) are unchanged. Measured at 1e6 carriers (b5): v3 digest ~0.2 s,
   v2 JSON seal 2.3 s, GPU motion 0.52 s (device 0.10 s), CPU motion 1.42 s; GPU and CPU v3 digests equal.
   Hall carrier cap: CarrierLimits {maximum_carriers (default 100000, <= 2^27), memory_budget_bytes (default 1 GiB)}
   given to the Model and bound into its configuration identity; the cap is admitted against the budget with a
   per-carrier host model of 4096 B (census, endpoints, paths, deposit support, staging and the JSON record built for
   native admission) + 256 MiB, else refused ("Hall carrier limits exceed the memory budget"). The cap is enforced
   in census order/identity, census digest, charge census and both providers with the sealed refusal text for the
   default cap. Evidence (b5): 1e5 default refuses 100002 carriers, a 1.5e5 cap admits and preflights them; GPU and
   CPU providers agree at 1.5e5 and 1e6 carriers (cap 2e6, 16 GiB budget).

Layout
  source/encoding.hpp            canonical binary encoder, streaming SHA-256
  source/field_contract.hpp      words, absolute Gauss/divB admission, C19 drive, binary mesh/state/ledger records
  source/solver_backend.hpp      backend interface, cpu-pcg backend
  source/field_model.hpp         Model (private prepare, backend + independent residual, endpoint tolerance, seal)
  source/census_charge.hpp       canonical carrier charge, DepositSupportBuilder (tolerance derivation)
  source/motion_receipt.hpp      v2 hall_driver motion types, MotionContext, MotionAssembly
  source/motion_provider.cpp     v2 motion provider TU (replaces sealed native_arena_motion.cpp in a link)
  source/owner.hpp, root_api.hpp, checkpoint.hpp   owner protocol, closed proposal, JSON export/restore
  source/gpu/gpu_backends.{hpp,cpp}                Device (one session), gpu1 field solver, GPU motion provider
  tests/test_hall_v2.cpp         CPU component qualification + motion dump
  tests/test_hall_v2_gpu.cpp     GPU qualification
  tests/sealed_motion_reference.cpp   same dump with the sealed 8dddec0 provider
  tests/repro_hall_node_plane_v2.cpp  Stage 1 node-plane reproducer on v2
  run_hall_v2.py                 build | cpu | gpu, receipted (frozen launcher, leases, sealed-input hashes)

Results of record: /home/treefor/projects/.workspace/outputs/hall-corrective-v2-20261030/b5/
  build/build-receipt.json, cpu-182242/cpu-receipt.json, gpu-182553/gpu-receipt.json
  (b4 = the same component evidence before item 5).
Coordinator integration and all-owner evidence: runtimes/pic-water-common-commit-v5-20261030.

Open
  - The sealed axial::push work-path fixed point (12 iterations, |dD| <= 1e-18 m) does not converge for one
    carrier of a random node-plane fixture (seed 22: electron 1e-9 m above a z node plane, 8.4e6 m/s) already in
    the DC seed field; the CPU and GPU providers refuse identically ("axial native work path did not converge").
    Sealed particle code, not the endpoint guard; 1 of 20 seeds x 500 carriers.
  - The v2 JSON-text record seal is still computed once per accepted step because native003 and native v5
    recognize only the v2 record (~2 us per carrier; 0.2 s at 1e5, 2.0-2.3 s at 1e6) and native v5 re-verifies it;
    removing it needs a native recognizer for the v3 record above.
  - Per iterate the census digest and the v3 digest are recomputed (identical inputs within one proposal); caching
    the census digest per proposal and hashing v3 only for the accepted iterate would remove ~0.4 s/step at 1e6.
  - Cross-architecture bit identity (sm_89, sm_120) not shown; the GPU descriptor binds the target architecture.
  - Device-side RHS and a device-resident compensated update/admission are not implemented.
