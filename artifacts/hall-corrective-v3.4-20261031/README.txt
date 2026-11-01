Hall corrective owner v3.4 -- 2026-10-31

v3.4 = v3.3 (sealed: runtimes/hall-corrective-v3.3-20261031, build s1) with the coupled step's per-carrier and per-mesh host
work moved to the device for the resident owner, and the divB check in coefficient form with a rounding floor. Entry points
unchanged: source/root_api.hpp (Owner, host arithmetic) and source/resident_owner.hpp (ResidentOwner, field state on the GPU).
Results of record: /home/treefor/projects/.workspace/outputs/hall-corrective-v3.4-20261031/RESULTS.txt.

0. Host-work modes (gpu::HostWork, ResidentOptions::host_work)
 device (default): the resident owner's coupled proposal runs on the device (section 1) when the Model uses the native v7
   motion forms (MotionForms::native_v7; the v3 forms keep the v3.3 path: their census v1 and record v3 are streaming hashes
   over the carrier texts).
 v33: the v3.3 host path of every moved item (GpuMotionProvider host arrays, host census admission, host DepositSupportBuilder,
   host bank_transfer, host census v2 and record v7, host preflight recomputation). Kept as the verification mode: both modes
   give the same decisions, refusal texts (same text, same refusal point), ledgers, physical state and digests (same Model ->
   same configuration identity -> equal block-Merkle state digests). The independent residual rule (section 2.2) is the v3.4
   rule in both modes (host evaluation in v33, device in device).
 HALL_HOST_WORK=v33|device in the environment sets the default (qualification switch for the unchanged v3.3 suites).
 The host Owner (Model::prepare) is unchanged except section 3 (divB) and stays the coordinator's host path.

1. The device coupled step (gpu/device_step_v34.{hpp,cpp}, gpu/step_kernels_v34.cu, gpu/support_v34.hpp, gpu/x87_v34.hpp)
 One DeviceStep per resident owner on the owner's session and stream; NVRTC programs compiled once per device and shared
 (gpu::Device::cached: push program, step program, native v7 forms, resident field and CG programs).
 1.1 Census admission (first statement of the proposal; v3.3: residents() copies, charge_of_census on the host)
   Host census (native v7: trial.native_census): one gather pass copies the 80-byte GapParticle of each row into pinned
   staging and checks the row identity (ascending IDs, nonempty packet/state texts) for all rows and for the residents; an
   unsorted census is sorted as v3.3 sorts it. Native's census v2 tables (StepInputs::census_packet_keys/state_indices =
   trial.census_v7(), in ascending native ID, the order of the sorted rows) are copied; without them the keys are computed
   from the texts on the host (v3.3 rule).
   Device census (native v7.1 census_v7_device: StepInputs::device_rows/keys/states/count, packet_of): no host pass.
   On the device: rows -> push carriers (SoA), newborn flags/flights/face start rule (k34_newborn_place, sink_core.hpp),
   charge_of_particles of the residents (k34_records_charge: the sealed check order and the Words CIC weights; records
   sorted stably by node; k34_charge_nodes: net words in carrier order, x87 gross step emulated), free_charge_guard and
   same_net_charge against the accepted slot's free charge (x87 emulated, k34_census_nodes), then the accepted slot
   admitted with that census charge (ResidentField::admit_device). Census v2: Forms::census_pack/census_check/block_list on
   the device; the digest header on the host (O(1)). Refusals: the v3.3 texts at the v3.3 points (identity and census-charge
   refusals at once; census v2 codes, the newborn list and the resident-neutral bound deferred to the motion admission).
 1.2 Motion (per fixed-point iterate)
   The v3 push kernels (MotionEngineV3::run_device) on the device trial midpoint (ResidentField prev buffers): no carrier
   upload per iterate, no output download. The fold of push_fold_v32.hpp on the device: exact K, W (scale dt), G_K, C_K, Q,
   G_W (exact_sum_v34.cuh: the native v7 superaccumulator, order- and partition-independent), continuity (div in ascending
   edge order, binary64 as gradient_transpose), E^ (node marks, exact maximum), the endpoint free-charge guard, the
   MotionAssembly::assemble operand checks; refusals resolved on the host from O(1) status words in the order of
   MotionAssembly::admit / GpuMotionProvider::motion / fold_record. Sinks compacted on the device; their canonical records
   downloaded (O(sinks)); SinkRecord/NewbornRecord Hall fields on the host (O(sinks + newborns)).
   Host waits: the device path waits with Session::synchronize (cuStreamSynchronize; the context is created with
   CU_CTX_SCHED_BLOCKING_SYNC, so the host thread sleeps). The v3.3 event waits (cuEventSynchronize on events created
   with CU_EVENT_DEFAULT) busy-wait on the host for the whole device time; the v3.3 path keeps them, the device path
   does not use them (run_device's edge total, the support record total, the residual's rhs/x marker).
 1.3 Deposit support of the endpoint continuity guard (v3.3: DepositSupportBuilder, x87, 1.6 us/carrier per iterate)
   Bit-preserving device form (no identity change): every long-double operation of the sealed builder is emulated exactly
   (x87_v34.hpp: 64-bit significand, nearest even; add/sub through an exact 192-bit window, mul through the exact 128-bit
   product, floor, conversion with the sealed nextafter step). Per carrier: Q, R, the node box, t1, t2 (k34_support_carriers);
   (node, carrier) records in carrier order, stable radix sort by node (k34_support_emit); per node the sequential x87 sums
   s1, s2 in ascending carrier order -- the sealed order -- and the finish (k34_support_nodes). The array stays on the
   device and the field step reads it in place. Qualification: emulator vs the host long double on 4e6 random operations and
   chains; emulated builder vs DepositSupportBuilder on 40-60 random sets (node planes, 1e-15 offsets, weights 1..1e12, both
   meshes): every node bitwise; device vs host builder on every parity step (test_v34).
   The 4-6x resident-process slowdown of HANDOVER-V3.4 item 2: moot for the support (no host loop). Measured in-process
   with test_v34 hostloop (the host DepositSupportBuilder at 1e4 carriers, an arithmetic loop and a 32 MB stream, before
   any CUDA context, with the blocking-sync context, with the resident owner and push engine allocated, after three
   coupled device steps): the same times in every phase, so the CUDA context does not slow host code. The cost that
   scaled with the device work in the resident process was the busy-wait of the event waits above (host CPU time equal to
   the device time), not a slowdown of host arithmetic (RESULTS.txt section 4).
 1.4 Bank transfer (v3.3: a fresh 2 x nn FreeCharge and a full-node scan per sink/newborn, O((S + B) nn))
   bank_transfer_sparse: deposit_charge restricted to the eight CIC nodes of each deposit, the same Words / x87 / binary64
   operations per node in the same deposit order (the v3.3 node scan order only orders different nodes); uploaded as a sparse
   list and scattered on the device (k_bank_scatter). O(sinks + newborns) host. dense_of() equals bank_transfer bitwise
   (every parity step).
 1.5 Record v7 (accepted iterate only) and preflight
   At the fixed point the canonical section arrays (endpoints 80 B, paths 64 B, sinks 136 B, births 24 B) are built on the
   device (k34_record_build / k34_birth_build) and hashed by Forms::block_list; the header (O(1)) and the array digests of
   the current, initial and endpoint charge (motion-array-binary64-v1, streaming: O(mesh) host) on the host.
   Trial-field digest of the binding (total-midpoint-E-B-binary64-v1, streaming, O(mesh) host): computed once per step,
   for the accepted iterate (MotionAssembly::set_trial_field_sha256 before the seal). The iterates before it are
   discarded; the device motion reads the trial buffers, not the digest; the fixed point leaves the trial buffers
   unchanged (ResidentField::midpoint_unchanged copies only when they differ). v3.3 computed it per iterate: same value,
   same decisions (the digest enters no check before the seal except its own format).
   ResidentOptions::host_record_arrays (default true) downloads endpoints/paths and the trial midpoint into the proposal for
   native v7's host admission (memcpy-only, O(N)); false for the native v7.1 device lane.
   Preflight (ResidentOwner::preflight_device): the v3.3 field/heat recheck (device leaves, both admissions, midpoint) and,
   recomputed on the device, the census v2 root, the residents' census charge against the accepted free charge (and the
   admission with it), the four section roots and the current root, the endpoint free charge from the endpoint records
   compared bitwise with the candidate slot; O(1) binding checks; the trial digest compared with the one computed for the
   binding and the stored midpoint compared bitwise with the trial buffers. The census span must be the one of the proposal
   (a different span is refused; in-place edits of the host rows between prepare and preflight are not re-read: native's
   census is immutable for the trial).
 1.6 Native v7.1 device lane (native candidate007.1, HALL-LINK-7.1.txt, API-071.txt; folded in)
   Forms: Hall v3.4 builds against the native v7.1 forms library (candidate007.1/forms: the v7 forms unchanged plus the
   block forms and Forms::block_lists; HALL_FORMS_V71). The NVRTC kernels embed the v7.1 forms_core.h (the v7 core with
   the SHA-256 rounds refactored; same digests).
   Item 1 (one session, one stream): DeviceStep::session() and DeviceStep::forms() are the owner's hall_gpu::Session and
   pf_v7::gpu::Forms (the coordinator constructs FullTimeDevice(session, forms) and attaches it to the arena). Every Hall
   write of the record arrays and E/B/J is complete when prepare returns (the device step synchronises at the seal).
   Item 2 (census in place): StepInputs::device_rows/device_keys/device_states/device_count = trial.census_v7_device()
   (rows, packet_keys, state_indices, carriers); packet_of = trial.carrier_packet (O(sinks + newborns) calls);
   census_states = the registered states, ascending (their digest is checked against the device census's states_sha256).
   Hall reads the rows in place (k34_rows_to_soa, census charge, census v2 recomputed on the device: Hall's census digest
   is its own, equal to native's census_sha256).
   Item 3 (record arrays): DeviceStep::record_view(): endpoints (80 B), paths (64 B), sinks (136 B), births (24 B),
   canonical layouts, ascending native ID, with their counts; unchanged until the next prepare of the owner
   (ResidentJointProposal::device_seal_generation() names the seal generation of the view that belongs to it).
   Item 4 (A1 inputs): record_view().E/B (the trial midpoint: the resident field's trial buffers, unchanged by the fixed
   point) and J (the push current), n_E = edges, n_B = faces.
   Item 5 (block forms): ResidentOptions::field_digest_form = blocks_v71: at the seal the roots of field_E, field_B
   (trial midpoint) and current_J are computed by one fused Forms::block_lists call; the binding's trial_field_sha256 :=
   pf_v7::trial_field_blocks_sha256(E root, B root) (MotionAssembly::set_trial_field_sha256), the header's
   current_A_m_sha256 := pf_v7::current_blocks_sha256(J root); the field owner's midpoint digest
   (MidpointView::total_field_sha256) is the same value. No trial-field download, no host hashing of E, B or J. The
   endpoint and initial charge digests of the header stay motion-array-binary64-v1 (API-007 defines no block form for
   them; O(nodes) host hashing at the seal). The sealed header digests are kept in the receipt
   (MotionReceipt::record_v7_array_digests) for the coordinator's record adapter (no host recomputation).
   Preflight in the block form: the J root and the trial E/B roots recomputed on the device and compared with the sealed
   ones and with the binding; a host copy of the current, when present, is hashed in the block form.
   Lane configuration: host_work device, host_record_arrays = false, field_digest_form = blocks_v71 (requires the device
   path; refused otherwise). The streaming forms (field_digest_form = streaming_v7, the default) remain for native v7.
   Item 6: steps with births or sinks and refusals run native's host path, which reads Hall's device arrays back.
   Item 7 (tests): test_v34 parity runs a third owner in the lane configuration in lockstep with the v33 and device owners
   (same decisions, refusal texts, accepted states and digests; its trial-field and current digests equal the v7.1 CPU
   forms of the device owner's host arrays; its record digest equals a host recomputation with those digests; the roots
   of Forms::block_lists, the step's own kernels and pf_v7::block_list agree); test_v34 tamper refuses one-ulp changes of
   the device current J, an endpoint and a path word of a lane proposal; test_v34 blocks compares the three root paths.
   The coordinator copy (coordinator-v7.1-copy, pair hall-v3.4 + native v7.1) runs native's admission in the device
   lane (FullTimeExecution::device, the default; run option --native-execution host|device|device_verified). Its
   qualification compares native device mode with native host mode on the same Hall records (identical step digests);
   native v7.1's device_verified mode cannot be used with device-only records (native defect, section 4a).

2. Resident field step (gpu/resident_field.cpp, gpu/resident_kernels.cu)
 2.1 Device inputs: plasma current, endpoint free charge and deposit support read in place (ResidentInputs::*_dev), sparse
   bank (ResidentInputs::bank_sparse). No upload per iterate.
 2.2 Independent residual, rule hall-v3.4-independent-residual-v1 (identity change; resident backend only):
   ||b - A x||_2 = sqrt(RN(exact sum of r_e^2)), r = b - A x with the sealed apply arithmetic, ||b|| likewise, admitted iff
   ||b - A x|| <= max(1e-30, 2e-13 ||b||); a nonfinite factor or an overflowing sum fails the check. v3.3 used the x87
   sequential sum (foil::norm); both are rounding-level evaluations of the same norm. Device: gather kernels k_res_scale/
   k_res_face/k_res_edge (no code shared with the CG kernels) and k_res_sum/k_res_finish; host (v33 mode): HostApply + the
   pf_v7::ExactSum CPU reference. Same bits in both modes (parity steps). The rule text is bound through the backend
   (FieldSolveBackend::residual_rule(); host backends keep the v3.3 text, so the Owner's identities are unchanged by it).
 2.3 Descriptor hall-v3.4-resident-field-descriptor-v1 (rules, kernel sources, forms core, exact-sum header).

3. divB check, coefficient form with the rounding floor (divb-coefficient-gross-v1; approved addition)
 3.1 Rule (host admit_constraints and device k_admit_cells):
     divB_c = sum_axes (B_a - B_b) fl(1/h)           (compensated words; fl(1/h) = 1/h, the curl coefficient magnitude)
     bound_c = 128 eps sum_axes (|B_a| + |B_b|) fl(1/h) + 8 eps^2 G_c
     G_c(new) = G_c(old) + ((|dB_xa| + |dB_xb|) fl(1/dx) + (|dB_ya| + |dB_yb|) fl(1/dy)) + (|dB_za| + |dB_zb|) fl(1/dz)
     for every step, dB_f = RN(B_f(new) - B_f(old)) of the compensated face words (binary64, this order; Model::prepare
     divb_gross_step and k_divb_gross give the same bits).
 3.2 Derivation. The update is B_f -= dt S_f, S_f = sum_k c_k midE_k with c_k = +-fl(1/h_t) (axial_mesh.hpp). Per cell the
   discrete divergence of the increment is sum_faces +-fl(1/h_n) S_f; an edge e of the cell lies on two of its faces with
   normals n, n' and tangential spacings t, t', and enters with fl(1/h_n) fl(1/h_t) and -fl(1/h_n') fl(1/h_t'); for the edge
   along the third axis {n, t} = {t', n'}, so the two products are the same pair of binary64 factors and cancel exactly in
   real arithmetic: div(curl) = 0 identically with the stored coefficients. The v3.3 check divided by h; then the pairs are
   fl(1/h_t)/h_n and fl(1/h_t')/h_n', which differ by fl(1/dz)/dy - fl(1/dy)/dz = 1.43e-16/(dy dz) on x and y edges (z
   edges cancel because dx = dy): the seed's 6.6e-15 T/m at the box-edge cells (0,22,5)/(44,22,5) (ratio 2.5e-3) came from
   this mismatch, not from the field. In coefficient form the check measures only rounding: the B words are updated in
   double-word arithmetic (v3.3), each step adds at most a few eps^2 |dB_f| fl(1/h) per face (compensated two-product and
   four-term sums, error <= ~2 eps^2 of the face increment), and the check's own evaluation adds O(eps^2) of
   sum |B| fl(1/h) (covered by the 128 eps term). Summed over steps, the accumulated divergence is bounded by
   ~4 eps^2 G_c; the floor 8 eps^2 G_c covers it with a factor 2 (and the O(n eps) relative rounding of G_c's own binary64
   accumulation, n the step count, ~1e-10 at 1e6 steps). Without the floor a cell whose |B| passes through zero after a long
   history would be refused for divergence that is pure accumulated O(eps^2) rounding (DIVB-ANALYSIS.txt section 1: the
   bound falls with the instantaneous |B|, the history does not). The floor is a rounding allowance, not a physical tolerance:
   8 eps^2 = 3.9e-31 of the cell's total variation.
 3.3 State: G_c per cell in ChargeLedger::divb_gross_T_m (empty = zero), in the host state digest (text
   hall-v3.4-divb-gross-ledger-v1 + one binary64 per cell) and the device block-Merkle digest (array divb_gross), in both
   checkpoints (key divb_gross_T_m; resident checkpoint hall-v3.4-resident-checkpoint-v1); zero at bootstrap ("bootstrap
   divB gross ledger must be zero"); restore takes the sealed value (the state digest covers it). The rule text is in every
   Model's configuration identity (identity change for the Owner too; physics unchanged).
 3.4 A09: 50 T/m single-face perturbation stays refused (bootstrap G = 0: margin as v3.3; mid-run the floor is ~1e-30 T/m).
   One-ulp G_c tamper refused (device slot: leaves; unsealed host and resident checkpoints: seal). RESULTS section 3.

4. Identities (v3.3 -> v3.4)
 Every Model: configuration identity binds divb-coefficient-gross-v1; the host state digest adds the divB gross ledger.
 Resident Models: backend descriptor (new kernels, residual rule), block-Merkle state digest hall-v3.4-field-state-merkle-v1
 with the divb_gross array.
 Native v7.1 lane configuration (field_digest_form = blocks_v71): the record's trial_field_sha256 and current_A_m_sha256
 in the block forms of API-071.txt 1.1, hence a different record v7 digest for the same motion (accepted by native v7.1;
 the record digest definition is unchanged). State digests and physical state are those of the streaming form. Unchanged: the deposit support values (bit-preserving), census v2, record v7, K, W, scales,
 sinks/births records, motion arithmetic, field arithmetic, refusal texts, heat and kinetic-ledger digests.
 Physical state: identical to v3.3 (the divB rule and the residual rule are checks; no state value depends on them).

4a. Approved decisions (2026-10-31)
 Approved with the v3.4 request: the divB check-scale fix (section 3: coefficient form, rounding floor 8 eps^2 G_c, the
 G_c state ledger); exact-sum deposit support with a documented identity change (not used: the bit-preserving x87
 emulation of section 1.3 keeps the sealed values, no identity change).
 D1 (a) Independent residual of the resident field step: rule hall-v3.4-independent-residual-v1 (section 2.2: exact sum
   of the squared residual, same tolerance; device in device mode, host in v33 mode, same bits). Alternative not taken:
   the v3.3 x87 sequential norm (O(mesh) rhs/x download and a host loop per solve).
 D2 (a) Block-form trial-field and current digests (section 1.6 item 5) for the Hall v3.4 + native v7.1 pairing
   (ResidentOptions::field_digest_form = blocks_v71); the streaming forms stay the default elsewhere. Alternative not
   taken: streaming everywhere (native hashes E/B/J on the host per admission).
 D3 (a) Coordinator copy, native v7.1 lane: the coordinator takes Hall's midpoint digest; the two independent checks are
   native's device recomputation of the block-form digest from Hall's E/B and Hall's bitwise preflight comparison of
   the stored midpoint with the trial buffers. Alternative not taken: the coordinator downloads E/B and recomputes
   (O(mesh) host per step).
 Native v7.1 (candidate007.1, sealed) device_verified mode reads a dangling record when the record carries device arrays
 (RESULTS.txt 10.1; reported for native v8); the coordinator-copy qualification uses native device mode against native
 host mode on the same records instead.

5. Remaining host work of the Hall coupled step (device mode) and why
 O(1) per step plus O(sinks + newborns), except:
  - host census input (native v7): one gather pass (80-byte copy + two text-length reads per row); none with the native v7.1
    device census;
  - host record arrays for native v7 host admission (host_record_arrays): endpoint/path download, memcpy only; none with
    the device lane;
  - O(mesh) once per step: the streaming trial-field digest of the binding (download of the trial E/B + SHA-256; native
    v7 requires the streaming form; the block form is computed on the device for native v7.1), the motion-array digests of
    the initial and endpoint charge at the seal (and of the current in the streaming form; API-007 defines no block form
    for the charge arrays); with host record arrays, the preflight's digest of the stored host midpoint (it guards the
    copy native v7 admits);
  - O(cells): the thermal owner (conductivity per proposal, Joule power of the metal cells per step: ThermalV3 is host-resident);
  - O(mesh/256): block-Merkle roots over the device leaves (32 bytes per 2 KB) and a 16-leaf host sample.
 Measured (RESULTS.txt section 1, host CPU per coupled step at 1e6 carriers, 200 / 50 um): v3.3 path 3536 / 3807 ms;
 device path, native v7 configuration 28.9 / 67.2 ms (0.029 / 0.067 us per carrier); native v7.1 lane configuration
 3.8 / 18.4 ms (0.004 / 0.018 us per carrier; flat in the carrier count from 1e4 to 1e6: what remains is per mesh).
 The step is then device-bound; device time per step, the timing table against v3.3 and the per-kernel split (latency-
 bound node walks at 200 um, FP64-side kernels at 50 um): RESULTS.txt sections 12 and 13.

6. Layout
  source/gpu/device_step_v34.{hpp,cpp}, step_kernels_v34.cu, step_abi_v34.h, support_v34.hpp, support_host_v34.hpp,
    x87_v34.hpp, exact_sum_v34.cuh                      the device coupled step
  source/gpu/resident_field.*, resident_kernels.cu, resident_abi.h   device inputs, device residual, divB ledger
  source/gpu/hall_gpu_push_v3.*                       run_device (outputs left on the device), shared program
  source/resident_owner.hpp                           device path (prepare_closed_proposal_device, preflight_device)
  source/field_contract.hpp, field_model.hpp, checkpoint.hpp, owner.hpp   divB rule and ledger, StepInputs tables
  source/motion_receipt.hpp, motion_provider.cpp      device receipts (MotionAssembly::device, header_v7, set_*)
  tests/test_x87.cpp       x87 emulation and the emulated builder vs the host (CPU)
  tests/test_v34.cpp       parity (both modes in lockstep), tamper, divb, blocks, hostloop, bench
  run_hall_v34.py          build | run | gpu (receipted; --env HALL_HOST_WORK=...); PF_BUILD_ISA_FLAGS (per-site ISA
                           flags, whitespace-separated, default empty) is appended to every C++ compile and recorded as
                           isa_flags in the build receipt (also run_coordinator_parity.py, coordinator-v7.1-copy/run_v7_1.py)
  tools/isa_hook_rebuild.sh  rebuild with the hook empty and compare the binaries with the builds of record
  run_acceptance.sh, run_bench.sh, run_coordinator_parity.py (coordinator v5 T1-T7 vs the v3.3 run), run_v33_rerun.sh,
    run_prerun.sh (SimQueue: 20 ns vacuum pre-run 200 um, windows 100/50 um)
  tools/bench_table.py, tools/compare_prerun_v33.py   bench tables; v3.4 pre-run vs the sealed v3.3 pre-run (state)
  tools/ptx_count.py, tools/kernel_timer/             measurement only, outside the Hall build and its qualification:
    static PTX instruction counts of the device programs (NVRTC, session options); per-kernel device time from a
    forwarding libcuda.so.1 (kernel_timer/README.txt; measure.sh runs the test_v34 bench lane at 1e5 and 1e6 carriers)
  coordinator-v7.1-copy/   copy of coordinator v7.1 (runtimes/pic-water-common-commit-v7.1-20261031, unchanged there)
    with the Hall v3.4 resident pair and the native v7 entry: bindings_v7_1.hpp / bindings_adapters_v7_1.hpp
    (PF_BINDING_PAIR_HALL_V34_NATIVE_V7: StepInputs with native's census v7 tables, resident Model), genesis/joint step
    templated on the field owner, run_v7_1.py (--hall v3.4), run_profile.sh (prof6 settings)
