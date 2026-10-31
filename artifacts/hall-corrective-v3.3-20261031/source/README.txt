Hall corrective owner v3.3 -- 2026-10-31

v3.3 = v3.2 (sealed: Power-Flow/artifacts/hall-corrective-v3.2-20261031) with the magnetic update and the DC seed kept
divergence-free to O(eps^2) (decision (b1) of the v3.2 divB analysis, outputs hall-corrective-v3.2-20261031/DIVB-ANALYSIS.txt):
 1. B update: B_f -= dt S_f with S_f = sum_k c_k midE_k as a compensated exact-product sum (Words: exact two-products,
    compensated four-term sum), B words updated in Words; v3.2 rounded the four-term sum to binary64 first, which left an
    O(eps) divergence per step that accumulates per unit simulated time. Host Model::prepare and device k_update_B use the
    same Words operations in the same order (bitwise host/device parity).
 2. DC seed (source/seed_v33.hpp): the 8dddec0 dc_initial with B0 = double-word curl of the seed vector potential (B, B_low)
    and the Coulomb vector solve's independent residual evaluated in double-word arithmetic (refinement cycles <= 8, target
    unchanged 2e-13 ||b||; the binary64 residual evaluation error exceeds the target on the 25 um mesh). Background force
    from the rounded B0.
 2a. 25 um mesh (360x4+1, 2,601,366 edges): the seed is now admitted (test_v32 seed-check: 5 refinement cycles, double-word
    residual 4.98e-13 within the 2e-13 ||b|| target; 130 s on one core); 200 um: 1 cycle (same CG result as before).
 3. Unchanged: every guard and tolerance, the divB check (division form, bound 128 eps sum|B|/h), the Gauss ledger,
    motion-work-gross-v1, the native v7 forms. The configuration identity binds the seed rule and the update rule text.
 4. Diagnostics: ResidentField::gauss_profile can return the per-cell |divB| and bound; the pre-run series records the
    per-cell maximum ratio, the smallest bound (its floor at zero crossings) and the ratio there; test_v32 compare reports
    ulp differences; test_v32 seed-check constructs a Model on a given mesh.
 Host-work items (device DepositSupportBuilder, device independent residual, device digests/preflight, the v7.1 documents)
 are deferred to v3.4.
 Remaining divergence: the division-form check differs from the curl's coefficient form by fl(1/dz)/dy - fl(1/dy)/dz =
 1.43e-16/(dy dz) on x and y edges; acting on the seed's vector potential it leaves ~0.3 eps relative divergence in some cells
 (largest per-cell ratio 0.0025 over 0.75 -> 2.75 ns). Proposal and measurements: outputs RESULTS.txt.

=== v3.2 README (unchanged below) ===

Hall corrective owner v3.2 -- 2026-10-31

Entry points: source/root_api.hpp (namespace hall_corrective_v3; Owner, host arithmetic) and source/resident_owner.hpp
(ResidentOwner, field state resident on the GPU). Derived from Hall corrective v3.1 (runtimes/hall-corrective-v3.1-20261031,
read in place, unchanged). Sealed inputs used unchanged: Hall 8dddec0 axial/foil field headers, kinetic_boundary.{hpp,cpp},
heat-v2 bindings, hall-gpu-resident-v1 common/ and push/, native candidate003 (coordinator v5 link), and the native v7
producer library native-full-time-candidate007/forms (compiled in place, API-007.txt sections 1-3 with amendment A1).
v3.2 changes:
  1. admitted Gauss ledger: the absolute Gauss admission is measured relative to a per-node ledger of the mismatch
     admitted by accepted steps (section 1);
  2. motion-work-gross-v1: the work tolerance of the motion on the gross exchange (MOTION-WORK-GROSS-V1.txt; native v7
     amendment A1) (section 2);
  3. native v7 motion forms selectable per Model (census v2, exact kinetic/work sums, motion record v7, configurable
     resident-neutral interval bound) (section 3);
  4. checkpoint save/restore of the resident owner (section 4).
Everything else (field arithmetic, every other guard and tolerance, refusal texts, the motion arithmetic) is v3.1.

1. Admitted Gauss ledger
 1.1 Quantities (per node n of a state; Words = compensated binary64 pairs, field_contract.hpp)
   Q_n   field charge: incidence of the E words times -M (field_charge)
   D_n   absolute Gauss defect: Q_n - free_n - conductor_n - exterior_n
   B_n   unchanged bound: eps0 max(1e-30, 2e-13 ||b_DC||_2) + 256 eps (activity_n + gross_n + |conductor_n| + |exterior_n|)
   L_n   admitted Gauss ledger (new state array, ChargeLedger::gauss_admitted_C, empty = zero)
 1.2 Accumulation (Model::prepare on the host, k_gauss_increment on the device, the same Words operations in the same
   order, so both ledgers are bitwise equal; verified every step in the parity tests)
   For every accepted step with endpoint continuity defect c_n = (after + bank_conductor + bank_exterior - before)_n -
   div(plasma)_n dt (the per-node value that the endpoint continuity guard bounds and admits; exactly 0 in field-only
   steps) and Gauss increment mismatch delta_n = (Q1 - Q0)_n - (div plasma + div metal + div outer)_n dt (the per-node
   value that the Gauss increment guard bounds and admits):
       L_n(new) = L_n(old) + (delta_n - c_n)            (signed, compensated words)
   The ledger enters the candidate only after both step guards have admitted the step; a refused step leaves no trace.
 1.3 Why D changes by exactly delta - c in one step
   D(new) - D(old) = (Q1 - Q0) - (after - before) - (conductor change) - (exterior change), with the conductor change =
   metal_delta dt + bank_conductor and the exterior change = outer_delta dt + bank_exterior (the ledger updates of the
   step). Substituting Q1 - Q0 = delta + (div plasma + metal_delta + outer_delta) dt gives
       D(new) - D(old) = delta - [(after + bank_conductor + bank_exterior - before) - div(plasma) dt] = delta - c.
   In floating point the difference is the rounding of the Words additions, |rho| <= ~4 eps^2 (|Q| + |net| + |banks|) per
   step; over 2e6 steps <= ~4e-25 of the local magnitudes, far inside the 256 eps term of B_n (5.7e-14 of the same
   magnitudes). Hence D_n(k) - L_n(k) = D_n(0) + r_n(k) with |r_n| negligible: the ledger carries the whole admitted drift.
 1.4 The admission (admit_constraints on the host, k_admit_nodes on the device)
       |D_n - L_n| <= B_n        for every node n
   Signs and magnitudes: L_n is signed and is subtracted from the signed defect BEFORE the magnitude is taken. The approved
   form |D_n| <= B_n + |L_n| follows by the triangle inequality (|D_n| <= |D_n - L_n| + |L_n|) and is reported as the
   implied ratio; the implemented form is the stricter one: it admits a defect only on the side, and by the amount, that
   the admitted steps moved it. The bound expression B_n is unchanged; no other guard or tolerance changed.
 1.5 Why it cannot admit charge that did not come through an admitted step
   Charge q placed at node n outside a step (in the free charge, a bank or the field) changes D_n by -q or +q and leaves L_n
   unchanged, because L_n changes only inside prepare by the values the step's own guards admitted. The admission then
   needs |D_n(0) + r_n -/+ q| <= B_n, i.e. |q| <= B_n + |D_n(0)| (about 2 x 1.73e-23 C): exactly the window that the v3.1
   guard allowed at bootstrap, independent of |L_n|. (With the literal |D_n| <= B_n + |L_n|, a charge up to 2|L_n| of the
   opposite sign would pass; the subtracted form refuses it: test "charge 2L at the node of largest |L|".)
   L_n itself cannot be supplied from outside: a bootstrap requires a zero ledger ("bootstrap admitted Gauss ledger must
   be zero"); the ledger is part of the field/charge state digest (host streaming record and device Merkle record) and of
   the checkpoint, so a changed ledger word changes the digest (refused at preflight and at restore). A checkpoint whose
   ledger, charge and field are edited consistently and resealed is outside what any state guard detects (as for a
   consistent edit of E and the charge ledger); the checkpoint seal is the protection there.
 1.6 Identities
   Host state digest hall-v2-field-state-v1: the charge ledger record is followed by text
   "hall-v3.2-admitted-gauss-ledger-v1" and the ledger words (nn words; an empty ledger encodes as zeros).
   Device Merkle digest hall-v3.2-field-state-merkle-v1: arrays E, B, E_low, B_low, free_net_hi/lo, free_gross,
   conductor_hi/lo, exterior_hi/lo, gauss_admitted_hi/lo, [residue]; descriptor hall-v3.2-resident-field-descriptor-v1.
   Configuration identity: binds admitted_gauss_ledger_rule, motion_work_rule and the motion forms text; hence every
   configuration, state, candidate and proposal digest differs from v3.1 (physical values do not, section 5).
   Checkpoint: key "admitted_gauss_C" (Owner and ResidentOwner); the "charge" record is unchanged.
   Constraints: max_gauss_C is the checked value max |D - L| (= the v3.1 value where L = 0); max_gauss_bound_C unchanged;
   new max_gauss_defect_C (|D|), max_gauss_ledger_C (|L|), max_gauss_node_ratio (|D-L|/B), max_gauss_defect_over_bound_plus_ledger.

2. Motion work tolerance (motion-work-gross-v1)
 Specification, derivation, refusal margins and measurements: MOTION-WORK-GROSS-V1.txt (native v7 amendment A1; native
 implements the same bits). Code: source/motion_work_gross_v1.hpp (pf_v7::ExactSum only), used by the CPU motion assembly
 (motion_provider.cpp), the device motion fold (gpu/push_fold_v32.hpp, a copy of the sealed hall-gpu-resident-v1 fold with
 this tolerance), MotionAssembly::assemble, the joint particle/Maxwell work bound (|K - plasma work| <= the same tolerance)
 and the complete energy bound (adds the admitted work tolerance; triangle inequality). The receipt keeps the scales
 (MotionReceipt::work_scales(); in no record). Refusal texts unchanged.

3. Native v7 motion forms (Model(..., MotionForms{native_v7 = true, neutral_stream_dt_max_s}); default: v3 forms)
 HALL-LINK.txt 1-3, 4 and A1 (3e):
   census v2: per proposal (MotionProvider census v2 cache, filled by the first admission after the binding checks, so the
     v3.1 refusal order is kept): packet keys, state indices (StepInputs::census_states = native's state table; empty =
     the sorted states of the census), canonical check (codes mapped to the v3.1 texts), digest with the forms library;
   kinetic_change_J and field_work_J: exact forms (pf_v7::ExactSum / pf_v7::field_work_J); a nonzero status refuses
     ("canonical motion work mismatch");
   motion record v7: MotionAssembly::seal_v7 / record_v7 (header on the host, section block lists of endpoints, paths,
     sinks and births with the forms library); MotionReceipt::record_v7_sha256(); no v2 JSON seal, no v3 record;
   resident-neutral stream: pf_v7::neutral_stream with the configured bound (CPU); the GPU provider refuses a step above
     the bound with a resident neutral carrier ("neutral authoritative interval");
   preflight (Owner, ResidentOwner): census v2 and record v7 recomputed and compared;
   A1 admission inputs: JointProposal / ResidentJointProposal::trial_E(), trial_B() (the accepted iterate's midpoint,
     total_field_sha256 = trial_field_sha256) and current_A_m().
 The forms are bound into the configuration identity. With the v3 forms the motion records are bit-identical to v3.1.
 Remaining host work per carrier on the v7 path: census v2 records and digest (once per proposal), the section block lists
 (accepted iterate), the exact sums of K and the work scales (host ExactSum in the fold). The device entry points of the
 forms library (Forms::census_check/block_list/kinetic_aligned/sum_products) are not used yet (open item).

4. Resident checkpoint (ResidentOwner::save / ResidentOwner::restore)
 checkpoint.json (schema hall-compensated-field-owner-v3, checkpoint hall-v3.2-resident-checkpoint-v1): configuration_sha256,
 dt_max_s, solver name and descriptor, motion forms, motion_mesh_sha256, field words, charge record, bank_residue_C,
 admitted_gauss_C, field_sha256 (resident Merkle digest), heat_sha256 with heat/ slices, kinetic ledgers and digest.
 Restore: the recorded configuration must equal the Model's, or Model::configuration_with(recorded dt_max, solver name,
 solver descriptor, motion forms) (same mesh, drive, material, rules, limits; the field state does not depend on these
 operands); census charge = ledger free charge (same rule as Owner::restore); host and device admission; the Merkle digest
 under the recorded configuration must equal field_sha256; heat and kinetic-ledger seals must match. The restored owner
 continues bit-identically (tests: same process, fresh process, Model with another dt bound). Not carried: the increment
 history of the controller's third-difference indicator (diagnostic; it restarts after two steps).

5. Parity with v3.1 (physical state)
 Where the guard is not binding the field words, free/conductor/exterior charge and heat are those of v3.1 bit for bit
 (the ledger is an added state array; the values of every other array are unchanged): vacuum 10,000 steps final physical
 digest equal to the v3.1 value 32f31874...; Hall v3 component tests and coordinator v5 T1-T7 with the physical state
 compared to the Hall v3 run (results: outputs RESULTS.txt). Identity changes: section 1.6 and the motion forms.

6. Layout
  source/                  v3.1 sources with the changes above; motion_work_gross_v1.hpp; gpu/push_fold_v32.hpp,
                           gpu/work_scales_v32.hpp
  tests/test_v32.cpp       a09 | ledger-run | restore-continue | dense | prerun | compare
  tests/test_resident.cpp, tests/test_hall_v3.cpp  v3.1 tests (state comparison includes the ledger)
  tests/plot_series.py     PNG plots of a prerun series (standard library only)
  run_hall_v32.py          build | run | gpu (receipted, frozen launcher, leases incl. loaned cores 0/1/3, sealed-input
                           hashes incl. the native v7 forms sources)
  run_coordinator_parity.py coordinator v5 T1-T7 on v3.2 vs the Hall v3 run (run-p3-221714)
  run_acceptance.sh        acceptance batch
  MOTION-WORK-GROSS-V1.txt  specification of the work tolerance
 Results of record: /home/treefor/projects/.workspace/outputs/hall-corrective-v3.2-20261031/RESULTS.txt

7. Open
  - Device paths of the v7 forms (census check/digest, section block lists, exact sums on the device) not used yet.
  - Native v7 admission of Hall v3.2 receipts: through the coordinator v7 binding pair (Hall v3.2 + native v7).
  - Deposit-rounding term of the work tolerance is rigorous only up to ~520 deposits per edge (MOTION-WORK-GROSS-V1.txt 4c).
  - The heat/source owner (electrode-heat v3) and other owners keep their own interval bounds.
  - Cross-architecture bit identity (sm_89, sm_120) not shown; a resident checkpoint restores on another architecture only
    through the configuration transfer of section 4 (descriptor differs), continuation identity there not shown.
