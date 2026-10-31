Coordinator v8, stage 1 -- source owner v7, owner-templated Hall step, whole-run restart, run interface -- 2026-10-31

Coordinator v7.1 (pic-water-common-commit-v7.1-20261031, build c6) is the base, read in place; v8 is a versioned copy
(namespace raikou::native_water::v8, files *_v8). Every owner operation, the step order, the step controller, predicate v7
and every identity text of v7.1 are kept; with benchmark options v8 reproduces v7.1 bit for bit (section 1). Stage 1:
  1. heat/source owner v7 in place of v6.1; its trial interval bound and interval policy as run options; the source-owner
     interval bound as a controller criterion; every coupled owner's cap configurable; dt calibration;
  2. the Hall owner type (host Owner or device-resident ResidentOwner) chosen per run, through one protocol;
  3. a whole-run atomic checkpoint and restore (A05/A06);
  4. a launch interface over pf-launch, a Hall-only pre-run in the coordinator Model and genesis at a Hall checkpoint.
Hall v3.4 and native v8 (both in progress elsewhere) bind as one table entry each (section 7).

Layout
  source/bindings_v8.hpp              Hall/native owner names, field-owner cap, coupled maximum min(Hall cap, source cap)
  source/bindings_adapters_v8.hpp     one section per pair: (Hall v3.2, native v7) parity pair, (Hall v3.3, native v7.1)
                                      default pair; templates over the Hall proposal type
  source/pair_hall_native_v7_forms.hpp  motion record v7 adapter (v7.1's, templated over the proposal type)
  source/hall_owner_v8.hpp            HallTraits<Owner|ResidentOwner>: accepted host view, midpoint, make, restore;
                                      bind_source (source v7 owner on a Hall slice, either owner type)
  source/identity_v8.hpp              v7.1 identities and setup with source owner v7; PilotSettings surface preset
  source/joint_step_v8.hpp            OwnersV8<H>, TrialV8<H>, prepare_v8 / can_publish_v8 / try_publish_v8
  source/source_interval_v8.hpp       accepted rate samples, source-owner interval bound
  source/step_controller_v8.hpp       v7.1 controller + criterion source_interval
  source/genesis_v8.hpp               genesis at the first C19 knot or at a Hall checkpoint, templated over H
  source/run_checkpoint_v8.hpp        whole-run checkpoint files, manifest, commit, verification, recovery
  source/common_commit_v7.hpp, exact_weight_sum_v8.hpp (v7.1 copies; the exact sum gains its restart record)
  tests/driver_v8.cpp (+ driver_v8_run.inc, driver_v8_checks.inc, driver_v8_fixtures.inc)
                                      modes run | checks | vacuum | prerun | hallcompare
  run_v8.py                           build LABEL [--hall v3.2|v3.3] [--native v7|v7.1] | run | baseline (v7.1 c6) | syntax
  tools/pf-coupled-v8                 launch interface (section 4)
  tools/compare_v8.py                 identity / physics / checks comparisons
  tools/calibration_v8.py             dt calibration (coordinator v7's procedure, standard library)
  tools/restart_tests.py              restart, crash, tamper and rank tests
  tools/lanes.sh, calibration-jobs.txt, gpu_owner_parity.sh, gpu_stage1_jobs.sh, gpu_launch_rev0.sh   run batches
  HALL-V35-REQUEST.txt                the Hall v3.5 interface needed for genesis at a pre-run checkpoint

1. Source owner v7 (benchmark parity)
  identity_v8.hpp includes source owner v7 (pic-electrode-heat-source-v7-20261031, verified against its sealed copy
  Power-Flow/artifacts/source-owner-v7-20261031 in every build) and its binding hall_source_binding_v7; the arrival adapter
  is pf_arrival_v11::source_v7. Run options: --source-dt-max (OwnerOptions::interval.cap_s, <= 1e-12),
  --source-emission-error / --source-desorption-error / --source-safety / --source-count-floor, or
  --source-interval-preset benchmark (no targets, v6.1) | production (1e-3 / 1e-3, safety 0.5, count floor 1e-3).
  The source Field of every step carries the configured bound (s7::Field(..., cap)).
  Parity with coordinator v7.1 (build p1 = Hall v3.2 + native v7 + source v7 at benchmark options, CPU backends;
  parity/p1-*-vs-v71-c6.json; tools/compare_v8.py identity against outputs/pic-water-common-commit-v7.1-20261031/c6-*):
    pilot 2000 steps at 10 fs (frames every 1000)        bit-identical
    qual 20 steps (births, ionization, sinks; frames)     bit-identical
    dense adaptive 0.5 ps (blob 8e10, 80 steps)           bit-identical
  (genesis.json, ledger-schema.json, deterministic.bin, all 113 deterministic ledger columns, step-digests.txt,
  step-receipts.jsonl, the control ledger without its timing column, probes.bin, frames.bin.)
  checks (p1-checks-cpu): 48/48, the v7.1 cases with the same names and outcomes; evidence differs only in the refusal
  prefix "coordinator v8:".
  Default pair (build s1/s2/s3: Hall v3.3 + native v7.1): checks 48/48 (s1-checks-cpu, s3-checks-cpu). Against the parity pair the deterministic ledger
  differs at rounding level only (Hall v3.3's compensated B update and double-word seed): probe fields within 1e-11 of the
  field magnitude, W_E within 2e-14 relative, every count (births, sinks, events, carriers) equal over pilot 2000 and qual
  20 (parity/s1-*-vs-p1-physics.json); identities differ by design (Hall configuration, native source identity).

2. Surface physics options (not enabled)
  PilotSettings::surface_preset / --surface: benchmark (default: the v6.1 surface; adds nothing to any identity),
  work-function (source v7's oxidized-SS304 work function with adsorbed-water coverage, sv7::ss304_work_function; flows
  through prepare_patch; a physics change, used only after approval), work-function+enhancement and full (refused at
  launch: they need per-(patch, class) arrival streams and sample_birth births (source v7 INTEGRATION D4) and native/Hall
  wall-response births and wall heat (B1-B4, C1-C4)). Defaults equal v6.1 behaviour.

3. Caps raised, the source-owner criterion and the dt calibration
3.1 Caps. --hall-dt-max (Hall Model dt_max, <= 1 ns), --source-dt-max (<= 1e-12, the source owner's hard limit),
  --neutral-dt-max (Hall MotionForms::neutral_stream_dt_max_s and native FullTimeMotionPolicy, the same value),
  --dt-cap (controller). Coupled maximum = min(Hall cap, source cap) = 1 ps at most. --charged-only-particle-cap requires a
  neutral bound >= the controller cap (refused otherwise); with neutrals needing the particle cap, the particle cap must
  lie below the neutral bound.
3.2 Criterion source_interval (step_controller_v8.hpp, source_interval_v8.hpp; active only with error targets): after
  genesis/restore and after every publication each slice owner evaluates Owner::rate_sample on the accepted state
  (epoch length = controller cap, sets the count floor); the step bound is min over slice owners of
  Owner::interval_bound({previous, current}) -- the owner's own function, recomputed by the owner in prepare from the same
  evidence (every slice owner prepares a trial when a policy is set, so each owner checks the step itself; a refusal
  "source interval bound" is a step-size refusal and is retried). Start-up (one accepted sample: genesis, first step):
  dt <= the particle cap (10 fs, the v7.1 step at which no bound applied) and the evidence's previous sample is the
  current one re-timed by that step (no measured change); recorded in the receipt ("source_interval": startup, owner
  records). Without targets the controller text, decision record and control ledger are those of v7.1 (parity). With the
  production preset the decision record, receipts and the control ledger gain the bound (source-interval.jsonl per step).
  Cost: +8 % step time in the pilot (rate samples of 16 slice owners).
3.3 Calibration (CPU backends, build s1; calibration/pilot-20ps.json, calibration/qual-2ps.json; tools/calibration_v8.py,
  the procedure of coordinator v7 README 6.2: errors against the dt/4 run, relative to the largest reference magnitude).
  Caps for the controller runs: Hall, source and neutral bound 1 ps, --charged-only-particle-cap, --dt-cap 1e-12,
  landing every comparison time (--sync-every).
  pilot, particle-free, 20 ps (0.7498-0.7698 ns), compared every 1 ps:
                          steps   dt mean/max (fs)   probe field   W_E       circuit work  emitted e-   wall s   dt_eq(W_E)
    fixed 10 fs            2000   10 / 10            1.09e-7       4.88e-8   2.78e-8       1.82e-5      142.9
    fixed 5 fs             4000   5                  2.19e-8       9.73e-9   5.55e-9       6.06e-6      311.2
    fixed 2.5 fs (ref)     8000   2.5                -             -         -             -            597.3
    lte_tol 1e-8           884    22.6 / 40.0        8.23e-7       3.67e-7   9.06e-8       6.38e-5       66.6    26.6 fs
    lte_tol 1e-8 prod.     884    22.6 / 40.0        8.23e-7       3.67e-7   9.06e-8       6.38e-5       84.2    26.6 fs
    lte_tol 1e-7           421    47.5 / 87.2        3.72e-6       1.66e-6   4.19e-7       1.42e-4       38.3    56.2 fs
    lte_tol 1e-6           203    98.5 / 190.3       1.66e-5       7.41e-6   1.89e-6       3.05e-4       19.6   118.4 fs
    lte_tol 1e-6 prod.     203    98.5 / 190.3       1.66e-5       7.41e-6   1.89e-6       3.05e-4       22.2   118.4 fs
    observed order: field quantities 2.00-2.01, expected emission 1.00 (explicit exposure). Limiting criterion: lte on
    every step but the start-up (particle cap), growth and landing steps; the source bound stays at the cap (the pilot's
    rates do not change measurably above the count floor) and the production preset gives the same steps. The default
    lte_tol 1e-8 takes 22.6 fs on average: 2.26x fewer steps and 2.1x less wall time than 10 fs, at 7.5x the 10 fs
    field error (dt_eq 25-27 fs; second order). lte_tol 1e-6: 9.9x fewer steps, 7.3x less wall time, field error 1.7e-5.
  qual, few particles (births every step, 20 keV electrons, sinks, ionization), 2 ps, compared every 0.1 ps:
                          steps   dt (fs)   probe field   W_E       circuit work  emitted e-   wall s
    fixed 10 fs            200    10        1.24e-8       2.02e-8   8.39e-7       2.94e-6      106.8
    fixed 5 fs             400    5         2.58e-9       4.03e-9   1.68e-7       9.82e-7      198.3
    fixed 2.5 fs (ref)     800    2.5       -             -         -             -            397.0
    caps 1 ps, particle cap 10 fs    200  10    1.24e-8 (deterministic ledger identical to fixed 10 fs)   105.1
    caps 1 ps, particle cap 100 fs    20  100   1.38e-6       2.16e-6   8.94e-5       3.77e-5       17.1
    (same, production preset)         20  100   identical                                           17.8
    With charged carriers the particle cap limits every step (100 fs: omega_c dt 0.018, 0.047 cells, no refusal); counts
    at 2 ps within Poisson (electron births 9534 vs 9565, H2O births 8711 vs 8754, events 41 vs 25 +- 5-6, sinks 18 vs 19).
    The particle cap stays 10 fs by default; this is the first carrier series above 10 fs.
  Limits that remain: the source owner's 1 ps hard limit; LTE (lte_tol) and drive (c_drive) govern particle-free steps.

4. Owner-templated step (host Owner and ResidentOwner)
  --hall-owner owner|resident. OwnersV8<H>, TrialV8<H>, the genesis and the run loop are templates over the Hall owner;
  HallTraits holds the four differences (accepted host view: the ResidentOwner downloads slot 0 once per accepted state
  for the source field, controller bounds, probes, frames and charge totals; midpoint accessors; bootstrap; restore). The
  resident Model carries ResidentBackend(descriptor of this device); --field cpu runs the resident step with the sealed CPU
  PCG on the downloaded right-hand side (host_cpu_pcg, the solve of the host Owner's cpu-pcg), --field gpu the device CG.
  Parity (Hall v3.3 resident owner, CPU motion, build s3; parity/op-*-owner-vs-resident.json): qual 20 steps, dense
  adaptive 30 steps, pilot 200 steps: frames, probes, control ledger (decisions, dt bits), heat, kinetic-ledger, RNG and
  arrival digests and every state ledger column bit-identical; the field digests differ by design (block-Merkle vs streaming
  record, hence native prefix and source seal); receipt scalars computed by device reductions in another order differ by
  rounding: field/plasma/Joule/circuit energies <= 3.5e-14 relative, divB bound 1.6e-16, terminal voltage 1.3e-16 (the
  energy residual columns, themselves rounding quantities, differ by up to 7.5e-8 of their size). The production
  configuration (resident + device CG + GPU push, qual 20 steps): push parity 40/40 iterates, ledgers closing.

5. Whole-run checkpoint and restart (A05/A06)
5.1 Content (run_checkpoint_v8.hpp): native.json (PicArena::save), hall/ (Owner::save or ResidentOwner::save: field words,
  charge ledger with the admitted Gauss ledger, heat slices, kinetic ledgers, Hall seals), source-<i>.bin (16 source v7
  owners), arrival.bin (arrival v1.1; its RNG is counter based: next step/ordinal), chemistry.json (scheduler clock and
  product RNG states), coordinator.json (step, accepted source seal and owner digests, receipt chain head, end time, run
  state with exact sums as integer records, controller scalars, source rate samples, statistics, genesis record),
  controller-ring.bin (the LTE history: four accepted (t, E, B)). Owner formats are the owners' own.
  manifest.json = {"body", "body_sha256"}: schema, step, time, ranks recorded (not required at restore), bindings, run
  options and controller fields (applied at restore), every owner's version (native source identity; Hall binding, owner
  type, configuration; source v7 configurations; arrival policy; chemistry configuration; coordinator version), every
  owner's state digests after the step, every file with size and SHA-256, the genesis checkpoint reference.
5.2 Atomicity: files written into ckpt-<step>.tmp and fsync'ed; manifest to manifest.json.tmp, fsync, rename; directory
  fsync; rename ckpt-<step>.tmp -> ckpt-<step> (commit point); parent fsync. Recovery (--resume ROOT) takes the newest
  committed directory and lists skipped temporary ones; its manifest and every file must verify or the restore refuses.
5.3 Restore: run options from the manifest; version checks; direct restore of every owner -- native first through
  PicArena::restore_full_time with the recorded boundary. Native v7.1 restores only a genesis-state (step 0) checkpoint
  and refuses a nonempty one (restore_deferred, A05/A11); then the verified replay runs: the run's genesis checkpoint is
  restored directly, the recorded steps are re-run, and the replayed owners must equal the manifest's digests and receipt
  chain head and every owner's own checkpoint restored beside them (Hall restore, 16 source loads, arrival restore,
  chemistry record, native file bytes, coordinator state and LTE ring). When native restores nonempty checkpoints
  (native v8, identity schema v8) the direct path is taken automatically.
5.4 Results (tools/restart_tests.py; host owner: build s2, outputs rs2-*, rs2-restart-results.json; resident owner:
  build s3, rs3r-*, rs3r-restart-results.json): 25/25 on each owner type.
    restores into fresh processes, continued bit-identically (deterministic ledger, digests, receipts, control ledger,
    probes): qual fixed 10 fs at steps 0 (direct), 10, 20 (replay); dense plasma blob (adaptive, mid-plasma) at 10 and 20;
    pilot adaptive with caps at 1 ps and the production source preset (LTE ring, rate samples) at 20 and 40.
    crash injected (exit 86) at four stages of the step-20 checkpoint write (owner files written, before the manifest,
    manifest temporary written, before the directory rename): recovery takes checkpoint 10 and continues bit-identically.
    refused: truncated native.json, flipped byte in hall/checkpoint.json, missing source-3.bin, arrival owner removed from
    the manifest (resealed), mismatched native version (resealed), mismatched Hall binding (resealed), flipped manifest
    byte, flipped arrival.bin byte with its digest updated (resealed: the arrival owner's checksum refuses), flipped LTE
    ring byte with its digest updated (resealed: the replayed coordinator state differs), flipped byte in the genesis
    checkpoint.
    rank independence: a restore launched with 3 ranks recorded (checkpoint written with 1) continues bit-identically.

6. Run interface, pre-run, genesis at a Hall checkpoint
6.1 tools/pf-coupled-v8: pf-launch options (--site, --cores/--cpus, --gpus, --ranks-per-gpu, --cores-available, --wait,
  --dry-run, --receipt) plus --hall-checkpoint, --resume/--restore-checkpoint, --mesh, --collision-grid, --dt-cap (raises the
  Hall cap, the source bound and the neutral bound with it), --particle-cap, --preset (production default), --hall-owner
  auto (resident with a GPU), --field, --push, --fixture, --window/--steps/--end-time, --checkpoint-every. The coordinator
  is one process: pf-launch runs it with --launcher exec on the selected cores (leases held, frozen root, receipt); the
  plan's rank count is recorded in every checkpoint. Inside a SimQueue job pf-launch uses the job's allocation.
6.2 prerun mode: the Hall v3.3 pre-run (DC seed, declared charge, landing grid m x 10 ps, equal sub-steps <= 100 fs) in
  the coordinator Model. To 20 ns on the resident owner (device CG): 194,026 steps, 631 s (SimQueue job 182, exclusive
  GPU), checkpoint prerun-20ns/checkpoint-20.000000ns. hallcompare (hallcompare-20ns/hall-compare.json): restored with the
  coordinator Model and compared with the Hall v3.3 20 ns checkpoint restored with its own pre-run Model (its recorded
  configuration 289715c1... reproduced): E, B, E_low, B_low, free/conductor/exterior charge and the admitted Gauss ledger
  (14,812 nodes), the clock and the 8 slice thermal states (temperatures, 19 cumulative energies): identical word for word.
6.3 Genesis at a Hall checkpoint (--hall-checkpoint DIR; genesis_v8.hpp): the Hall owner is restored (provisional with the
  predicted census, real with the actual census; equal seals), the native genesis binds the restored field/heat digests and
  time, sources/arrival/chemistry bind to it; genesis.json records the checkpoint. From the 20 ns checkpoint the genesis is
  refused by the sealed owners: Hall and native tie the Hall revision to the native accepted step (Hall
  prepare_closed_proposal and preflight_final, native motion admission, V3 owner_values_match), and a native genesis is at
  step 0 while the restored owner is at revision 194,026 ("genesis owners disagree with the accepted native boundary").
  Decided (2026-10-31): Hall v3.5 restore_as_genesis (HALL-V35-REQUEST.txt). Stage 1 validates the path at revision-0
  Hall checkpoints (the hall/ part of a run's genesis checkpoint), launched through pf-coupled-v8 -> pf-launch ->
  run_v8.py (parity/start0-*):
    host owner, qual fixed 10 fs 30 steps (start0-qual-cpu vs rs2-qual-ref) and pilot adaptive, caps 1 ps, production
    preset, 60 steps (start0-pilot-cpu vs rs2-pilot-ref); resident owner, the same pair (start0-qual vs rs3r-qual-ref,
    start0-pilot vs rs3r-pilot-ref): deterministic ledger, control ledger, probes, field, heat, kinetic-ledger, RNG and
    arrival digests identical to the fresh-genesis runs; native prefix and source seal differ (genesis provenance:
    "coordinator-v8/genesis-at-hall-checkpoint", genesis.json "hall_checkpoint"). Ledgers close: water inventory and H/O
    nuclei exact every step, energy residual <= 5.0e-15 J, charge change <= 2.3e-41 C, Gauss <= 0.0099 and divB <= 2.5e-5
    of their bounds. Resident owner with the device CG (start0-pilot-gpu, 60 adaptive steps): ledgers close likewise.
    A whole-run checkpoint of a checkpoint-started run restores (step 20, replay from its own genesis checkpoint) and
    continues bit-identically (parity/start0-qual-cpu-restore-20.json).

7. Stage 2 (what Hall v3.4/v3.5 and native v8 must provide; the coordinator side is one table entry each)
  Hall v3.4: HALL_BINDINGS 'v3.4' + PAIRS ('v3.4','v7.1') with its adapter section (the device step behind the same
  ResidentOwner protocol; its resident checkpoint schema with the divB gross ledger is restored through HallTraits
  unchanged); Session& and Forms& of the resident owner and the device census inputs, so the coordinator attaches native
  v7.1's FullTimeDevice lane on the same stream (HALL-LINK-7.1.txt); block-form digests (field_digest_form blocks_v71).
  Hall v3.5: restore_as_genesis (HALL-V35-REQUEST.txt): seal under the recorded configuration, heat-only declaration,
  geometry and limits transfer for vacuum states, revision rebased to 0 with the clock kept, provenance bound into the
  genesis identity. Then the 20 ns start and the pilot from it.
  Native v8: NATIVE_BINDINGS 'v8' + the pair; nonempty restore_full_time (done in candidate008, identity schema v8) makes the
  direct restore path the default (the replay remains for v5-schema checkpoints); a test of the direct path at steps 10/20
  goes with the pair entry. Chemistry v3, wall-response births (source v7 INTEGRATION B1-B4), metal species and
  merge/split enter through that pair as well; the surface presets of section 2 need B1-B4 and D4.
  Source v7 INTEGRATION D2 (sink -> WallHit) and D6 ledger columns follow with the wall model.

8. Decisions taken (engineering; none changes physics at benchmark options)
  - Parity is shown on the v7.1 owner pair (Hall v3.2 + native v7) with source v7 at benchmark options; the default pair is
    Hall v3.3 + native v7.1 (identities and rounding-level values differ by design).
  - Source criterion: evidence from accepted rate samples with epoch length = controller cap; start-up rule above (dt <=
    particle cap, previous sample re-timed); every slice owner prepares a trial when a policy is set. Production preset
    1e-3/1e-3/0.5/1e-3 (approved 2026-10-31); the driver default stays benchmark (v7.1 parity), the launch tool default
    is production.
  - v8 additions to records only when active (decision record, control ledger columns, criterion code 15), so benchmark
    runs keep the v7.1 bytes.
  - Restart: owner formats unchanged; manifest sealed by SHA-256 of the body; the rename of the checkpoint directory is the
    commit point; recovery never falls back past a corrupt committed checkpoint (it refuses); verified replay from the
    genesis checkpoint for native v7.1 (approved 2026-10-31); ranks recorded, not required.
  - Resident owner: accepted host view by slot-0 download once per accepted state; ResidentOptions sample_seed fixed
    (digest leaf sampling only).
  - Launch interface: one coordinator process (--launcher exec) on the selected cores; ranks = the plan's rank count.
  - Pre-run comparison restores the reference with its own recorded configuration and the plain foil planes of the Hall
    test fixtures (its motion mesh digest).

9. Builds and results of record (outputs /home/treefor/projects/.workspace/outputs/pic-water-common-commit-v8-20261031/)
  p1  Hall v3.2 + native v7 (parity pair)      p1-pilot-fixed-2000, p1-qual-fixed-20, p1-dense-ctl, p1-checks-cpu
  s1  Hall v3.3 + native v7.1 (default pair)   s1-*, cal-* (calibration), prerun-20ns
  s2  s1 + restore of the run settings          rs2-* (restart tests, host owner)
  s3  s2 + hallcompare mode                     hallcompare-20ns, op-* (owner parity), rs3r-* (resident restart), start0-*
  s4  build of record: final sources (one comment changed after s3); binary byte-identical to s3 (build receipts);
      s4-qual-fixed-20 bit-identical to s1-qual-fixed-20. tools/pf-coupled-v8 uses s4 by default.
  parity/*.json, calibration/*.json, rs2-restart-results.json, rs3r-restart-results.json; build-*/build-receipt.json.
  Every compile and run inside run-frozen; CPU runs on the lent cores 0/1 (lease, taskset, nice 15) and cores 6/7
  (any-core.sh); GPU runs through SimQueue (job 182 exclusive pre-run; jobs 184, 186 shared).
