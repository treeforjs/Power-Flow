Canonical motion and v2 axial thermal integration ports

This module has 123 passing component checks on sealed 1f9cf7e. It preserves
every earlier Hall/Maxwell/thermal-v1 source and evidence artifact. It is ready
for the assigned native/root owners to bind to their private full-time outer
transaction. It does not implement a competing arena extension or claim that
the physical coupled application is complete.

Files for the owners
  outer-transaction-contract.txt: full-time first chemistry -> transport/source/
    heat -> second chemistry design and the required staged two-RNG root API.
  motion_receipt.hpp, native_arena_motion.cpp, MOTION-API.txt: sealed immutable
    motion receipt, native owner admission contract, census and receipt hashes.
  axial_thermal_v2.hpp: actual electrode_heat_v2::Coupler for both foil bodies in
    each axial slice, using the sealed v2 material callbacks and source identity.
  axial_joint_v2.hpp: C19 current-control Maxwell/motion/v2-heat preparation and
    common preflight; no arena or RNG publication inside this component.
  qualification.json, motion-api-seal.json: exact hashes, counts, tested scope
    and remaining owner integration. evidence/ contains command receipts and
    validation summaries, including exact restart file comparison.

Joint port call sequence
  The native outer owner supplies its PRIVATE post-first-chemistry census and
  MotionBinding. accepted_step equals the field/heat accepted revision; begin
  equals their accepted physical time. staged_arena_revision and arena hash
  describe the private inventory after first exposure. Do not use rounded view.
  Select the material with material_binding_sha256() and mesh with
  motion_mesh_sha256(mesh, native_options). Both are checked before work.
  Call prepare_axial_joint(mesh, fields, accepted_field, native_options,
      private_census, common_binding, thermal, C19_trace, active_layer_ids).
  The function prepares every particle iterate from that same census and every
  heat trial from accepted thermal state. It converges on IDENTICAL rounded
  midpoint E/B vectors, within the existing 16-iteration cap. The trial field
  seal in the returned MotionReceipt is exactly the midpoint of accepted and
  proposed Maxwell fields. No nearby-stage source or particle work is relabeled.
  The result owns the field proposal, field receipt, MotionReceipt and v2 heat
  trial. Native owner admits the motion through its typed private operation,
  then root prepares second chemistry with the first stage's two staged RNGs.
  Root checks preflight_axial_joint AND every other owner before a final commit.
  publish_field_heat_prevalidated performs only the field/heat nonthrowing swaps;
  it does not commit a native population, source or RNG. No callbacks, I/O or
  fallible computation may be inserted between the root's common swaps.

Thermal choices and ownership
  Namespace electrode_heat_v2; explicit evaluated_solid_1500 scenario, range
  298.2..1500 K, rho_ref=8000 kg/m3, fixed reference volume Cp-as-Cv approximation.
  This is a conditional thermal material, not mechanics or a vapor EOS. Mechanical
  qualification and higher-temperature scenarios are not extended here.
  The mesh slice height is dz exactly once. Each cell's admitted Joule power is
  the same Maxwell edge quadratic assembled into the actual v2 finite-element
  owner. Accepted temperature changes the next conductivity. No axial thermal
  diffusion is claimed. This version has zero source patches and no desorption.
  Source Wfield and Qwall must not be folded into Joule. Actual source-v3 coupling
  and population/ownership handoffs belong to the separate source/native owners.
  Restoring all slices builds private candidate owners and swaps only after all
  files and common clocks validate. Root still owns atomic checkpoint directory
  publication, cross-owner restart and crash consistency.

Evidence
  thermal-v2-001: 23 checks. Actual C19/Joule transfer, material/namespace binding,
  changed conductivity and next Maxwell response, failed and stale trials,
  aggregate restore rollback and exact fresh-process restart of 9 field/heat
  fixture files. Sixteen variable steps; total Joule 1.2994128239385296e-7 J,
  Tmax 300.0000006747159 K. This is not a canonical arena checkpoint.
  motion-v2-002: 52 checks (supersedes the 30-check motion-002 port). Native Boris momentum exactly matches in the constant
  total-field fixture; exact carrier identities/weights and same packet
  kinematics survive; neutral drift included. Path continuity 4.56705e-29 C,
  work defect -2.21727e-27 J. Manufactured variable-step C/2-transport-C/2
  refinement ratio 4.0000924752. Actual source/SSA are not in that manufactured
  fixture. Probed linker traps observe zero native Poisson/D-K-D transport calls.
  joint-motion-v2-001: 48 checks (supersedes joint-v2-001). Twelve variable C19 steps, exact shared midpoint field
  seal, at most 3 nonlinear iterations; maximum particle/field work defect
  5.99820e-23 J. Maximum combined balance 4.95789e-15 J includes explicit thermal
  roundoff already covered by the unchanged v2 budget. Tmax 300.00000050603694 K.
  Test-only transient carrier recurrence exercises the real motion proposal
  implementation; it is not a substitute persistent arena or actual SSA run.
  motion-001 failed to compile because a test helper named apply collided with
  std::apply through argument-dependent lookup. Renaming that test helper fixed
  compilation. Its full source snapshot and log remain in the output directory.
  motion-v2-001 also had a test variable name collision; the failed snapshot is
  retained. The neutral boundary tests passed after renaming that test variable.
  No physical, numerical or resource bound was loosened.

Reproduce each family on CPU6 under the exclusive existing lease
  env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    CUDA_VISIBLE_DEVICES= taskset -c 6 nice -n 15 \
    flock -n /home/treefor/sim-workspace/powerflow/outputs/pic-profile/testing-core.lock \
    python3 development/hall-canonical-driver-20261029/run_thermal_v2.py UNIQUE_RUN
  Use run_motion.py or run_joint_v2.py with a different unique label for the other
  families. Runners snapshot inputs, retain failures, enforce 2 GiB/180 CPU s/
  240 wall s per command, and use the frozen runtime/MFEM dependencies. They
  compile serially with -ffp-contract=off -fno-fast-math. No GPU or production run.

Immutable bindings
  bindings.json records original archive/source hashes. The 17 canonical v6.2
  identity files and 25 heat-v2 files match their immutable copies and originals.
  native-base-verification.json is the arena owner's 359-check v6.2 qualification,
  including sanitizer results; it is not counted among these 123 component checks.
  Never link the copied baseline arena alongside the native owner's extension.

Remaining physical application gates
  Native owner's full-time private outer state, bootstrap, mode-aware reaction
  staging and semantic checkpoint replay; root's actual parent SSA on both sides
  with both RNG streams; source-v3 live population/ownership handoffs; root common
  all-owner rollback/crash consistency; wide exterior and axial convergence.
  The present finite field box is an operator fixture. It has not replaced the
  root's wide exterior or established a physical axial instability result.

Current native-owner helper
  Include native_joint_binding.hpp and use prepare_native_joint on the actual
  PrivateTransportTrial after first SSA. Root need not duplicate Hall binding
  setup. See root_joint_prepare_example.hpp for both exposures and common commit;
  root still owns that coordinator. The helper uses staged_arena_revision(),
  exact private census/prefix and the actual proposed v2 heat-state seal. The
  helper/example syntax-instantiated against the current native and ParentSubflow
  draft APIs; runtime qualification awaits their sealed implementation.
  Motion receipt v2 adds a mandatory independent kinetic_policy_sha256. The new
  kinetic_boundary helper retains native finite-target xy/short/feed policy. It
  provides neutral first-hit/remaining-time/kinetic proposals without removing
  inventory or assigning wall heat. Tests widen the field box from 9 to 17 mm
  and increase axial air padding while native particle bounds remain identical.
  Source population and charged partial-flight admission are still explicit gates.

Caller correction after sealed 8ea1798
  native_joint_binding.hpp now wraps the actual motion audit in the native owner's
  exact arena_total_field_motion_receipt_v1 envelope and uses the actual provider
  proposal seal as transport_trial_sha256. The private arena prefix stays bound
  separately in the provider's motion input. Native acceptance does not replace
  actual field/heat proposal preflight. The numerical component files and prior
  evidence are byte-identical. Eight focused envelope/identity checks and root
  example syntax pass in native-envelope-001; the 123 component matrices were
  not repeated. Full native linking and two-exposure runtime qualification are
  still pending. The native API and its owner implementation were not modified.
