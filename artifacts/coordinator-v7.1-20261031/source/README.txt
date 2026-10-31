Coordinator v7.1 -- Hall v3.2 + native v7 binding pair, mesh selected at launch, less coordinator host work -- 2026-10-31

Coordinator v7 (pic-water-common-commit-v7-20261031, build b5) is read in place. v7.1 is a versioned copy (namespace
raikou::native_water::v7_1, files *_v7_1). The step controller, predicate v7 (common_commit_v7.hpp, unchanged), every owner
operation, identity text and refusal point are those of v7. Changes:
  1. owner pair (Hall corrective v3.2, native v7) as one binding table entry plus one adapter section;
  2. the mesh is selected at launch (Hall MeshSpec, source-setup mesh and collision grid from one option);
  3. the coordinator's own per-step host work is reduced with identical outputs (section 4);
  4. a carrier-load fixture and per-step timings of the coordinator and the device for host profiles (section 5).

Layout
  source/bindings_v7_1.hpp           Hall/native owner names, owner step bounds (coupled_dt_max)
  source/bindings_adapters_v7_1.hpp  one section per admitted pair: native admission predicate, motion record digest and
                                     admission, census states, Hall StepInputs and Model construction
  source/pair_hall_v32_native_v7.hpp record adapter of the (Hall v3.2, native v7) pair (section 1)
  source/exact_weight_sum_v7_1.hpp   exact sum of nonnegative binary64 values (ledger weights, section 4)
  source/identity_v7_1.hpp, joint_step_v7_1.hpp, genesis_v7_1.hpp, step_controller_v7_1.hpp, common_commit_v7.hpp
  tests/driver_v7_1.cpp              modes run | checks | vacuum; fixtures pilot, qual, dense, load
  run_v7_1.py                        build LABEL [--hall v3|v3.1|v3.2] [--native v6|v7] | run | baseline (coordinator v7 b5)
  tools/profile_report.py            host profile tables (section 5); tools of v7 (compare, calibration, receipts, ...)

1. Owner pair (Hall v3.2, native v7)
  run_v7_1.py: HALL_BINDINGS['v3.2'] (hall-corrective-v3.2-20261031/source, checked against its sealed artifact copy; compiles
  against the native v7 forms headers), NATIVE_BINDINGS['v7'] (native-full-time-candidate007, source identity 3f88fab5...,
  core archive libpic_native_v7.a sha256 af2bd9e7..., forms/include), PAIRS[('v3.2', 'v7')] (define
  PF_BINDING_PAIR_HALL_V32_NATIVE_V7). Default pair of v7.1: (v3.2, v7); (v3, v6) and (v3.1, v6) stay admitted.
  bindings_adapters_v7_1.hpp section of the pair:
    Model(..., dt_max, MotionForms{native_v7 = true, neutral_stream_dt_max_s = 1e-14}) (the native default bound);
    Hall StepInputs carry census_states = trial.census_v7().states (native's census v2 state table);
    motion record digest = MotionReceipt::record_v7_sha256(); exposures check compares it;
    admission = pair_hall_v32_native_v7.hpp: FullTimeMotionRecordV7 from the accepted Hall proposal: header fields of the
      binding, budget (Hall's exact kinetic change and field work), census v2 digest (Hall's), array digests
      (pf_v7::array_sha256 of current, initial and endpoint charge: O(mesh)), kinetic policy, frame, neutral bound; spans:
      endpoints (receipt GapParticle array), paths (receipt CarrierPath array, same 64-byte layout as FullTimeMotionPath,
      checked by static_assert of size and field offsets), sinks and births (converted, O(sinks + newborns)); A1 inputs
      trial_E, trial_B (the proposal's midpoint) and current_A_m; record_sha256 = Hall's digest. Native recomputes the
      exact sums, the motion-work-gross-v1 tolerance, the four section roots, the record and census digests and refuses
      any difference. No per-carrier conversion or hashing on the coordinator side.
  Hall v3.2 + native v7 as built: 47/47 checks (c1-checks-cpu: the 31 v6.1 cases, 14 V7, 2 V7.1); 48/48 from build c4 on
  (one more V7.1 case).

2. Parity with coordinator v7 (CPU backends; parity/*.json, compare_v61_v7.py)
  v7.1 on the v7 pair (Hall v3.1 + native v6), build c2: 2000-step pilot at 10 fs (with frames) and qual 20 steps are
  bit-identical to v7 b5 (genesis, ledger schema, deterministic ledger, step digests, frames).
  v7.1 on the new pair (Hall v3.2 + native v7), build c1, same windows, and the adaptive dense run (blob 8e10, 0.5 ps):
    identical physics: every deterministic ledger column except gauss_max_abs_C, frames (densities, E, potential), probes
    (E, B at the probe entries and the field energies), the control ledger, and the controller's decisions (dt bits,
    criterion, every decision input);
    gauss_max_abs_C differs by design: Hall v3.2 reports the checked value max |D - L| (absolute Gauss defect minus the
    admitted Gauss ledger; equal to v3.1's max |D| while the ledger is zero, i.e. at step 1 only);
    identity changes (all by design): Hall field configuration digest (v3.2 rules: admitted Gauss ledger,
    motion-work-gross-v1, native v7 motion forms) and the field-state digest (it adds the admitted Gauss ledger); native
    source identity (3f88fab5...) and the native accepted state prefix (motion record v7, census v2 and exact K in the
    motion rows, API-007.txt 5); hence the coordinator configuration digest, the genesis common clock and genesis
    source seal, every step's common clock (binds the configuration) and the accepted source seal. Unchanged: heat,
    kinetic-ledger, chemistry RNG and arrival digests.
  Hall v3.2 stops at about 4.15 ns on the absolute divB guard (under analysis by the Hall owner); every window used here
  ends before 0.8 ns.

3. Mesh selected at launch
  --mesh FILE (resolution-series-prep-20261031/meshes/mesh-*.json, schema resolution-series-mesh-input-v1; default: the
  200 um mesh) sets the Hall MeshSpec, the source-setup mesh (one source slice per foil per axial cell, patches and streams)
  and, with --collision-grid mesh (default), the chemistry collision grid (transverse mesh cells x 5 axial bins);
  --collision-grid fixed200 keeps 45 x 45 x 5 (identical chemistry binning across a mesh series). The built Hall mesh is
  checked against the file (node/edge/face/cell counts and the dx, dy, dz, zmin bits); the file digest, MeshSpec and
  encode_mesh_sha256 are recorded in run-summary.json (mesh).
  Foil planes (identity_v7_1.hpp hall_geometry): the 200 um special case of v6/v7 is replaced by one rule for every mesh:
  the binary64 value nearest to the plain plane value whose Hall cell coordinate ((p - w/2) + side/2)/h is an exact node
  plane or lies on the metal side (searching toward the metal), so a conductor-bank deposit at a foil face puts no
  rounding weight on a vacuum node. It gives the v6 planes on 200 um bit for bit (V7.1 check) and the same values on
  100 and 50 um (left x faces 2.0e-3, lower-foil outer face 1.8e-3; other planes plain).
  Admitted by Hall: the uniform 200, 100 and 50 um meshes (resolution-series INDEX.json).
  Runs (build c6, CPU backends): qual fixture on 100 um (c6-qual-100um, 20 steps: 948 electron and 874 H2O births, 12
  ionization events, the 10 foil-face sinks) and on 50 um (c6-qual-50um, 10 steps: 468 / 433 births, 10 events, 10 sinks),
  mesh checked against the file, water inventory and H/O nuclei exact every step, Gauss defect <= 0.0063 of its bound; the
  load-fixture profiles on 50 um (section 5).

4. Coordinator host work removed (identical outputs)
  Per step v7 traversed the accepted census four times (native rate_population in both observe() calls, in the
  controller and in the ledger totals), walked the packet map twice, credited exact weights with native ArenaCount
  (BIGNUM) per carrier (up to four credits per carrier: species weight and H/O nuclei), built the controller's field
  bounds with a 27-cell loop per cell, sorted all 8N plasma node terms globally and hashed the chemistry RNG state for
  both the digest line and the receipt. v7.1:
    one traversal per accepted state (joint_step accepted_population): species check of every packet and carrier,
      ledger totals, the next decision's controller inputs and frames use it; observe() repeats the per-carrier species
      check only for an accepted state other than the verified one (cache keyed by the native accepted state prefix;
      the channel-binding and catalog checks run on every call);
    exact weights and nuclei by an exact superaccumulator (exact_weight_sum_v7_1.hpp: units of 2^-1074, 35 limbs,
      round to nearest even once): the same exact values and the same rounded ledger values as the ArenaCount credits
      (V7.1 check against ArenaCount on 2000 random sums; every run below bit-identical);
    controller: field bounds only when charged carriers exist, as three 1-D neighbourhood-max passes (same values);
      plasma node terms bucketed by node with only each node's terms sorted (same order within a node, same sums; build c6);
    one chemistry RNG digest per step for the digest line and the receipt;
    probe field energies (O(mesh)) every --probe-energy-every steps (default 1, as v7).
  --host-work v7 runs the v7 form of this work in the same binary (profile baseline). Evidence (parity/*.json): builds c4,
  c5 and c6 (both modes) against c1 for qual 20 steps (frames), the adaptive dense run and the 2000-step pilot: genesis,
  deterministic ledger, step digests, frames, probes, control ledger and step-receipts.jsonl byte-identical.

5. Host profile
  Load fixture (--fixture load --load-carriers N: 49 % electrons and 49 % H2O+ co-located with weight 1, 2 % H2O of weight
  1e6, in the gap; chemistry limits sized to N), GPU field solver and GPU push, dt 10 fs, one core, GPU locks taken before the
  core lease; runs prof-* (build c5, --host-work v7 = the v7 host work) and prof6-* (build c6, v7.1 host work), 8 / 4 / 2 steps
  at 1e4 / 1e5 / 1e6 carriers on the 200 and 50 um meshes. Tables: outputs .../profile/ (tools/profile_report.py: median of
  the steps after the first, every ledger, sub-phase and coordinator timing column grouped by owner). Peak RSS 1.9 GB (1e6,
  200 um) and 2.4-2.5 GB (1e6, 50 um); genesis 1-6 s.
  Host step per owner (ms, v7.1):  1e4/200  1e4/50  1e5/200  1e5/50  1e6/200  1e6/50
    total                             150    1194      964    1724    22305    9781
    chemistry subflow                  19      36      297     177    15895    2639
    Hall                               94    1037      423    1278     4038    4889
    native                             11      25      116     112     1219    1169
    coordinator                        12      42       71      79      717     670
    source                             12      42       12      40       15      44
    arrival                           0.2     2.0      0.3     0.4      1.4     1.2
    device (field solve + push run)   6.6      44       31      64      262     299
  Host time per carrier at 1e6: 22.3 us (200 um) and 9.8 us (50 um); slopes (1e5 -> 1e6) and the O(mesh) items:
  HOST-WORK-REMAINING.txt.
  Coordinator own host work per step, v7 -> v7.1 (same binary pair of modes, identical outputs: step digests, receipts and
  deterministic ledgers compared at 1e4/200 um and 1e6/50 um):
                                   1e4/200  1e4/50  1e5/200  1e5/50  1e6/200  1e6/50
    v7                                24.0    55.9    233.7   237.8   2150.9  2236.5
    v7.1                              11.6    42.1     71.0    78.5    717.3   670.2
    of which observe() x2 (v7: 4.6+3.6 .. 386+387; v7.1: 3.7+1.9 at 1e6), controller inputs (851 -> 299 at 1e6/200 um),
    totals (517 -> 405), receipt (1.9 -> 0.05). Remaining in the coordinator: one population traversal (0.72 us/carrier) and
    O(mesh) digests (about 20 ms on 50 um).

6. Builds and results of record (outputs/pic-water-common-commit-v7.1-20261031/)
  c1 first build of the pair (Hall v3.2 + native v7); c2 the v7 pair (Hall v3.1 + native v6) for parity with v7; c3 timing
  columns; c4 host-work removal; c5 --host-work switch and load chemistry limits; c6 bucketed plasma sums (build of record,
  Hall v3.2 + native v7). Checks: 47/47 (c1), 48/48 (c4, c5, c6; c2 on the old pair). Parity files: parity/*.json; profile
  tables: profile/*; HOST-WORK-REMAINING.txt (copy).
  Cores: lent cores 0 and 1 only (cores 2 and 3 were returned during the work); GPU runs took both GPU locks before the
  core lease.

7. Open
  - The coordinator's remaining O(N) traversal needs a typed native population view (native v8) or a device reduction to
    become O(1) on the host; the ledger energy sums need an exact or device form to stay bit-stable when made incremental.
  - Hall v3.2 stops at about 4.15 ns on the absolute divB guard (Hall owner analysis); the heat/source owner and the sealed
    neutral stream keep the 10 fs step bound (coordinator v7 README section 7).
  - Mesh series beyond 50 um (25 um) is not admitted by Hall; collision-grid refinement per member is a configuration choice
    recorded in run-summary.json (--collision-grid).
