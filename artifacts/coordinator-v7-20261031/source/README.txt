Coordinator v7 -- adaptive step length for the coupled water-breakdown model -- 2026-10-31

Coordinator v6.1 (pic-water-common-commit-v6.1-20261031, build a2) is read in place. v7 is a versioned copy
(namespace raikou::native_water::v7) that adds a step controller, common-commit predicate v7 and the records that
go with them (plan section 8, artifacts/coupled-pilot-plan-20261030/PLAN.md). Owners, step order, every owner
operation, every identity text and every refusal point are the v6.1 ones. The Hall and native owners are selected
by one owner binding per build (section 0); the build of record links Hall corrective v3.1 (configurable step cap,
default 1e-14 s, default identities equal to Hall v3) and native candidate006 (native v6).

0. Owner bindings
  source/bindings_v7.hpp (+ bindings_adapters_v7.hpp) is the only coordinator source that names the Hall and the
  native owner; run_v7.py build LABEL --hall NAME --native NAME chooses the owner pair (HALL_BINDINGS: v3 =
  hall-corrective-v3-20261030, v3.1 = hall-corrective-v3.1-20261031; NATIVE_BINDINGS: v6 = candidate006; PAIRS: the
  admitted pairs with their adapter section in bindings_adapters_v7.hpp), sets the include paths, owner objects,
  libraries and the pair's adapters, verifies the Hall sources against their sealed artifact copy, refuses a
  link that mixes headers of two Hall bindings, and records the binding in the build receipt, the binary
  (PF_HALL_BINDING, PF_NATIVE_BINDING), every run receipt and run-summary.json. The coordinator reaches the field-owner
  step cap only through bind::field_owner_dt_max(model) and bind::make_hall_model(dt_max, ...), which detect the
  configurable-cap interface (Model::dt_max(), Model(..., ThermalDeclaration, double dt_max)); without it the cap is
  the v3 literal and a requested dt_max other than 1e-14 refuses. The coordinator's owner maximum is
  bind::coupled_dt_max = min(Hall cap, heat/source owner bound 1e-14). The planned Hall v3.2 + native v7 pair drops in
  as one HALL_BINDINGS, one NATIVE_BINDINGS and one PAIRS entry plus one adapter section (section 7).
  --hall-dt-max sets the Hall model's configured cap (v3.1: enters the Hall configuration identity when != 1e-14).

Layout
  source/bindings_v7.hpp, bindings_adapters_v7.hpp   owner bindings (Hall and native owner selection, cap interface)
  source/step_controller_v7.hpp   controller: inputs from the accepted state, criteria, dt ladder, decision record
  source/common_commit_v7.hpp     predicate v7 (V3/native-v6 conditions with dt <= configured cap)
  source/joint_step_v7.hpp        prepare_v7 (dt <= configured cap, decision bound into the trial), can_publish_v7
                                  (predicate v7), try_publish_v7; otherwise the v6.1 step
  source/identity_v7.hpp, genesis_v7.hpp   v6.1 copies (namespace and refusal prefix only)
  tests/driver_v7.cpp             modes run | checks | vacuum; fixtures pilot, qual, dense
  run_v7.py                       build | run | baseline (the v6.1 a2 binary, sha256-verified)
  tools/compare_v61_v7.py         identity (files, deterministic ledger columns by name, digests, frames) and checks
  tools/calibration.py            dt series comparison (errors, observed order, equivalent fixed step)
  tools/verify_receipts.py        sealed step-receipt chain and its cross-checks against the ledgers
  tools/refusals.py               owner step-size refusals per attempted omega_p dt
  tools/ensemble.py               seed-ensemble tally of the stochastic counts
  tools/write_summary_v7.py       qualification/QUALIFICATION-SUMMARY.json
  tools/ledger_csv.py, ledger_view.py   v6 tools, unchanged

1. Step controller (step_controller_v7.hpp)
  Inputs, computed before every step from the accepted owners only (start-of-step state) and from the last accepted
  step's records:
    omega_p,max  max over Hall nodes of sqrt(sum w q^2/(eps0 m) / V_node): the operand of the Hall motion plasma
                 guard (same CIC point deposit, same factor order) over every charged carrier of the census; the
                 terms of a node are summed in ascending order of value, so the result does not depend on how the
                 census is partitioned (ranks) or ordered.
    omega_c,max  max over charged carriers of |q/m| |B|, |B| bounded per cell by the face maxima of the accepted B
                 (per component) over the 3x3x3 cell neighbourhood of the carrier's cell.
    v_max, a_max max |u|/gamma over all carriers; max |q/m| |E| over charged carriers (edge maxima, same
                 neighbourhood).
    chemistry    largest unresolved upper-bound parent-event rate of the last step's exposures (events_upper /
                 exposure duration, coverage ledger) and the coverage policy's per-exposure limit (0.01).
    drive        C19 current at t (Hall model drive).
    lte          max over E and B of max_i |x'''_i| / max_i |x_i|, x''' = 6 x the third divided difference of the
                 last four accepted (t, E, B) states.
    births       expected electron births per second of the last step (expected count / w0 / dt).
  Criteria (predicates on a candidate dt; + - * / sqrt only, so the decision bits are the same on any IEEE-754
  machine and independent of rank count and field/push backend):
    plasma        omega_p,max dt <= c_plasma                         (Hall motion guard: 0.1)
    gyro          omega_c,max dt <= c_gyro                           (Hall push guard: 0.05)
    displacement  dt min(c, v_max + a_max dt) <= c_disp h_min        (Hall push guard: 0.2 cell)
    chemistry     chem_rate dt/2 <= c_chem x 0.01                    (coverage exposure refusal: 0.01 events)
    drive         |I(t+dt) - I(t)| <= c_drive max(|I(t)|, 1 kA)
    lte           particle-free steps only: lte dt^3/12 <= lte_tol (before four states exist: dt <= particle cap)
    births        particle-free steps only: births rate x dt <= birth_threshold, else dt <= particle cap
  Bounds: dt <= cap (configured, <= the coupled-owner maximum), dt <= particle_cap while carriers exist (charged, and
  neutral while neutrals_need_particle_cap),
  t + dt <= end of the C19 drive.
  Decision: dt is the largest rung of the ladder particle_cap x 2^(-j/8) (fixed binary64 mantissas; rung 0 is the
  particle cap bit for bit, rung 8 half of it, negative rungs above it up to the largest rung <= cap) admitted by
  every criterion. Growth: at most 8 rungs (x2) per step, and a rung above the
  previous one must pass every physical criterion at 1.05 dt (hysteresis). Reduction is immediate. Optional landing
  times (--sync-every): a step that would pass a landing time ends on it (two equal steps when one full step would
  leave less than half a step). The limiting criterion is the one that refused the next larger rung (or cap,
  growth, hold, sync); its guard value at the chosen dt and its limit are recorded.
  Refusal: a preparation refused by an owner step-size guard (plasma guard, gyro/displacement guard, neutral
  displacement guard, device motion budget, work-path convergence, excessive path, coverage exposure limit,
  reaction resource limit, thermal capacity, the Hall midpoint fixed point within 16 iterations) is retried 8 rungs
  (one octave) lower inside the same step transaction (a refused preparation leaves every owner unchanged; T3 and
  the V7 checks); up to 6 retries; every attempt (dt, decision digest, refusal text) is recorded; no growth for the
  next 4 accepted steps. Any other refusal stops the run as in v6.1.
  Fixed mode (default; --dt): dt is the given value on every step, no retries (the v6.1 behaviour). Calibration
  only: --fixed-grid-retry retries a step-size refusal at half the step and completes the grid interval with the
  next step.

2. Common-commit predicate v7 (common_commit_v7.hpp)
  accepts_closed_step(e, current, candidate, clock) is the V3 contract with "0 < dt <= cap" in place of
  "0 < dt <= 1e-14", where StepClock = {cap, coupled-owner maximum, decided dt, decision digest} requires
  0 < cap <= coupled-owner maximum, the decided dt bit-identical to the expected step's dt and a digest for the
  decision record. accepts_closed_or_birth_step(..., clock) is the native v6 birth-step admission (candidate006
  adapters/common_commit_v6.hpp) with that closed-step predicate. Every other condition is unchanged and in the
  same order. With cap = 1e-14 and a step clock that carries the decision of the step, the decision equals the V3/v6
  predicate's for every dt <= 1e-14 (V7 checks evaluate both on the same observations). prepare_v7 refuses dt above the configured cap and a configured cap above the
  coupled-owner maximum (bind::coupled_dt_max: min(Hall Model::dt_max(), heat/source owner bound 1e-14)).

3. Records
  ledger.bin, ledger-schema.json, deterministic.bin, step-digests.txt, genesis.json, frames.*: the v6 layouts and
    schemas, unchanged (dt_s is the step's dt).
  control-ledger.bin, control-schema.json: one row per step: dt, rung, landing, criterion code (control-schema
    criterion_codes), criterion value and limit, the first decision before retries, attempts, every controller input,
    predicted guard values at dt (omega_p dt, omega_c dt, displacement, chemistry events, LTE, births), the
    continuous bound of every criterion, the guard values the Hall step measured (plasma_dt, gyro, displacement)
    and the controller time (t_control_ms, not deterministic).
  step-receipts.jsonl: the sealed step receipts. Receipt k = canonical JSON of {step, previous receipt digest,
    controller configuration digest, decision record (inputs, dt, rung, criterion, value, limit, bounds) and its
    digest, attempts, common-clock digest of the step (which binds the dt bits and the accepted before-state),
    dt and its bits, criterion, value, limit, after-state digests (native prefix, field, heat, kinetic ledgers,
    source seal, RNG, arrival)}; its digest = sha256(previous digest + receipt text). The genesis link binds the
    controller digest and the genesis native/field/source digests; head and genesis link are in run-summary.json.
    The decision digest is also bound into the trial (StepClock) and checked by predicate v7 before publication.
    tools/verify_receipts.py recomputes the chain and cross-checks dt bits and digests against the ledgers.
  probes.bin, probes-schema.json: accepted end-of-step E and B at ten fixed mesh entries (gap, foils, feed, blob)
    and the field energies W_E, W_B; row 0 is genesis.
  The controller configuration is not part of the owner configuration identity (so that a 10 fs run keeps the v6.1
  identities); it is bound into every receipt and the genesis link.

4. Fixtures
  pilot   empty plasma, oxide-water sources, w0 = 1e6 (pf-pilot-1 setup).
  qual    qualification fixture (births every step, ionization events, ten sinks in step 1).
  dense   qual plus a seeded plasma blob: --blob-pairs co-located electron/H2O+ pairs (default 400) of weight
          --blob-weight (default 1.6e11 per carrier), uniform in native x 2.6-3.4 mm, y 5.1-5.9 mm, z -3..3 mm
          (collision cells without authored water), electrons Maxwellian at --blob-ev (default 3 eV), ions at
          0.025 eV (fixed generator seed). Peak omega_p dt at 10 fs: 0.0646 at 1.6e11, 0.0457 at 8e10, 0.030 at 3.45e10. The blob's
          Debye length (about 0.1 um) is far below the 200 um cell; it is a calibration load for the plasma
          criterion, not a resolved plasma.
  --seed-offset N adds N to the arrival, collision-clock and product seeds.

5. Build and run (cores: lent 0, 1, 3 or Power-Flow 6-7; one single-threaded job per core)
  flock .../pic-profile/testing-core<K>.lock taskset -c <K> nice -n 15 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 runtimes/arch-base-20261028/run-frozen /usr/bin/python -I run_v7.py build LABEL
  ... run_v7.py run LABEL run NAME --fixture pilot --steps 2000                       (fixed 10 fs, v6.1 behaviour)
  ... run_v7.py run LABEL run NAME --fixture dense --dt-mode adaptive --window 2e-12 --sync-every 1e-13
  ... run_v7.py run LABEL checks NAME
  Controller options: --dt-mode fixed|adaptive, --dt, --dt-cap, --hall-dt-max, --particle-cap, --charged-only-particle-cap, --c-plasma, --c-gyro, --c-disp,
  --c-chem, --c-drive, --lte-tol, --birth-threshold, --growth-rungs, --max-retries, --sync-every, --window,
  --fixed-grid-retry. run_v7.py refuses new work on cores 0-3 once RETURN-CORES-0-3 exists.

6. Results of record (/home/treefor/projects/.workspace/outputs/pic-water-common-commit-v7-20261031/)
  Build of record b5 (build-b5/build-receipt.json; Hall v3.1 + native v6). b1-b3 link Hall v3, b4-b5 Hall v3.1 (its
  CPU path and default identities equal Hall v3). Builds differ only in: b2 --fixed-grid-retry; b3 the
  neutral-carrier particle cap, core 2 admission, the end-time tolerance, the calibrated c_plasma default (0.03);
  b4 the owner bindings (Hall v3.1, --hall-dt-max); b5 binding pairs, the coupled-owner maximum and the ladder anchored
  at the particle cap. The calibration runs cal-* used b2 (cal-denseM-* b3), the statistics runs st-* and ens-* b3;
  none of their decisions depends on a later change (cap = particle cap = 1e-14; c_plasma explicit or irrelevant;
  the b2 fixed runs only took one extra step past the window, not used), and b5-denseM-ctl reproduces cal-denseM-ctl
  (6.1). qualification/QUALIFICATION-SUMMARY.json
  collects the records below; calibration/*.json are the full comparison tables (tools/calibration.py).

6.1 Parity with coordinator v6.1 (CPU backends; compare_v61_v7.py; build b5 = Hall v3.1 + native v6)
  pilot, 2000 steps, dt fixed 10 fs (b5-pilot-fixed-2000) vs v6.1 a2-pilot-cpu-2000         bit-identical
  pilot, 2000 steps, adaptive controller (b5-pilot-adaptive-2000; cap chosen on 2000/2000)   bit-identical
  qual, 20 steps, fixed (b5-qual-fixed-20) vs a2-qual-cpu-a                                 bit-identical
    (genesis.json, ledger-schema.json, deterministic.bin, step-digests.txt, all 113 deterministic ledger columns by
    name, frames 3/3 for the pilot). The same holds for b3 (Hall v3) and b4 (Hall v3.1): b3-/b4-pilot-fixed-2000,
    -adaptive-2000, -qual-fixed-20 (qualification/earlier-builds), b1-pilot-cpu-2000-pinned, b1-qual-cpu-a-pinned.
  Across Hall bindings: the default dense-medium controller run on Hall v3.1 (b5-denseM-ctl) is bit-identical to the
    same run on Hall v3 (cal-denseM-ctl, b3), step digests and control ledger included. With the Hall cap configured to
    100 fs and the coordinator cap at 10 fs (b5-qual-hallcap-adaptive-20) the deterministic ledger equals the fixed
    10 fs run (only identities differ). A controller cap above the coupled-owner maximum is refused at configuration
    (b5-cap-refusal: "controller cap must lie in (0, largest step of the coupled owners]").
  checks (b5-checks-cpu; also b3, b4): 45/45. The 31 v6.1 cases come first with the same names, outcomes and evidence; only three
    refusal texts differ (prefix "coordinator v7:"; the dt guard text "dt within the configured step cap" in place of
    "the unchanged native full-time dt guard"). 14 new V7 cases: predicate v7 against the V3/native-v6 predicate on
    the same observations (they agree on every evaluation; they differ only where v7 adds a condition); refusal of a
    dt differing from the decision by one ulp, of a cap above the field-owner maximum, of a missing decision digest,
    of dt above a configured 5 fs cap and of a cap above 1e-14 (owners unchanged); 9.17 fs and 5 fs steps publish; the
    ladder; a step-size refusal after the Hall proposal leaves every owner unchanged and the retried 5 fs step equals
    a direct 5 fs step bit for bit; non-step-size faults are not retried; controller inputs and decision bit-identical
    for the census in 1, 2, 3 and 8 parts in reversed order; the controller's omega_p dt equals the Hall motion
    plasma guard value of the step bit for bit (0.0296236 at 4.585 fs, dense fixture); a real Hall plasma-guard
    refusal at 10 fs on an 8x heavier blob (omega_p dt 0.183) leaves owners unchanged and the controller's dt
    (1.62 fs, omega_p dt 0.0296) publishes.
  Step receipts: tools/verify_receipts.py verified every chain checked (qualification/receipts-*.json); a run-loop
    retry is in smoke-qual-retry (b1: attempt 1 refused at 10 fs, retried at 5 fs, then 4 steps without growth).
  Adaptive runs that chose the cap on every step (cal-pilot-ctl 2000 steps, cal-qual-ctl 200 steps) are bit-identical
  to the fixed 10 fs runs.
  Controller cost: 0.75 ms per step particle-free, 2.2 ms with 5.6e3 carriers, 4.4 ms with 1.9e4 carriers (about 1 %
  of the step); memory: four accepted (E, B) states for the LTE estimate (2.8 MB).

6.2 Calibration: same physical window at fixed 10, 5, 2.5 fs and with the controller (CPU backends)
  Relative error = largest deviation from the 2.5 fs run over the comparison times / largest magnitude of the 2.5 fs
  value (a probe field component: / largest magnitude of its field vector at that probe). p = observed order from
  e(10)/e(5) against the 2.5 fs reference. dt_eq = fixed step with the controller's error under the fitted law.
  pilot (particle-free, 20 ps, comparison every 1 ps; steps 2000 / 4000 / 8000; controller 2000, cap every step)
                                    10 fs      5 fs       p
    probe field (all entries)       1.10e-7    2.19e-8    2.0
    field energy W_E                4.88e-8    9.73e-9    2.01
    cumulative field change         4.18e-8    8.35e-9    2.00
    cumulative circuit work         2.78e-8    5.55e-9    2.00
    cumulative Joule heat           8.7e-12    2.6e-12    1.23
    cumulative thermal change       5.1e-7     1.2e-7     1.73
    expected emitted electrons      1.82e-5    6.06e-6    1.00   (explicit start-of-step source exposure: first order)
    LTE estimate: max|x'''|/max|x| = 7.18e33 s^-3 (median) in all three runs (dt-independent: it measures the
    solution); relative LTE per step at 10 fs median 6.0e-10, max 3.4e-9; with lte_tol 1e-8 the particle-free bound
    is 14.3-25.6 fs in this window (not reachable while the cap is 10 fs).
  qual (births every step, ionization, sinks; 2 ps, every 0.1 ps; steps 200 / 400 / 800; controller 200, cap)
    probe field (all entries)       1.24e-8    2.58e-9    -
    field energy W_E                2.02e-8    4.03e-9    2.00
    cumulative circuit work         8.39e-7    1.68e-7    2.00
    expected emitted electrons      2.94e-6    9.82e-7    1.00
    expected ionization events      3.1e-2     1.3e-2     0.54   (includes stochastic secondaries)
    census kinetic energy           4.6e-4     2.1e-4     0.18   (dominated by the ionization realization)
    counts at 2 ps (10 / 5 / 2.5 fs): electron births 9513 / 9540 / 9565, H2O births 8750 / 8691 / 8754, ionization
    events 35 / 25 / 25 (expected 36.2 / 34.8 / 35.3), sinks 18 / 26 / 19; every difference to the 2.5 fs run within
    1.7 sigma (Poisson). Charge collected -1.602e-18 C in every run (the ten authored sink electrons of step 1; later
    sinks are emission electrons of weight 2e-23).
    Seed ensemble (calibration/ensemble-qual-1ps.json; 8 seeds, 1 ps): ionization events / expected events
    140 / 137.2 = 1.02 +- 0.09 at 10 fs and 130 / 136.4 = 0.95 +- 0.08 at 2.5 fs (difference 0.56 sigma); electron
    births 38144 / 38144, H2O births 34842 / 34824, sinks 86 / 85: no dt dependence of the stochastic counts.
  dense-light (blob 3.45e10 per carrier, omega_p dt = 0.030 at 10 fs; 0.815 ps, every 0.1 ps; landing every 0.1 ps)
                                    10 fs      5 fs       c=0.02     c=0.01     p
    steps                           88         169        137        255
    blob probe field                1.98e-4    3.97e-5    7.17e-5    8.84e-6    2.00
    field energy W_E                5.11e-5    1.02e-5    1.85e-5    2.26e-6    2.00
    cumulative plasma work          4.87e-5    9.72e-6    1.75e-5    2.17e-6    2.00
    census kinetic energy           2.37e-6    4.74e-7    8.54e-7    1.06e-7    2.00
    dt_eq                                                 6.33 fs    3.23 fs
    (c=0.05 chose the cap: identical to 10 fs.) Error law: the extrapolated error at 10 fs, 2.1e-4 over omega_p t =
    2.45 rad, is 0.096 (omega_p dt)^2 per radian, close to the implicit-midpoint frequency error (omega_p dt)^2/12.
  dense-medium (blob 8e10, omega_p dt = 0.0457 at 10 fs; 0.5 ps, every 0.05 ps; default controller c_plasma = 0.03)
                                    10 fs      5 fs       controller p
    steps                           50         100        80 (70 at rung 5 = 6.48 fs, 10 landing steps)
    blob probe field                1.37e-4    2.75e-5    4.97e-5    2.00
    field energy W_E                9.07e-5    1.82e-5    3.28e-5    2.00
    cumulative plasma work          8.09e-5    1.62e-5    2.92e-5    2.00
    census kinetic energy           4.20e-6    8.41e-7    1.52e-6    2.00
    dt_eq of the controller 6.33 fs; counts within 1.6 sigma of the 2.5 fs run.

6.3 Owner step-size refusals against the attempted omega_p dt (calibration/refusals.json; 27 dense runs, 2667
  attempts, 14 refused): the Hall closed proposal "did not reach the identical authoritative midpoint within 16 iterations"
  (an ulp-level cycle of the bitwise fixed point) at
    omega_p dt < 0.02: 0 / 1114;  0.02-0.03: 2 / 435;  0.03-0.04: 3 / 624;  0.04-0.05: 7 / 304;  0.05-0.06: 2 / 190;
    0.065 (fixed 10 fs, heavy blob): 1 in 5 steps.
  The Hall plasma guard (0.1) was never reached by a controller decision. In adaptive runs every such refusal was
  retried at half the step inside the transaction and the retry published.

6.4 Safety factors (defaults of ControllerConfig)
  c_plasma = 0.03   frequency error (omega_p dt)^2/12 <= 7.5e-5 per radian (6.2: error law confirmed), fixed-point
                    refusals about 0.5 % of attempts (each costs one extra preparation and 4 steps without growth),
                    a factor 3.3 below the Hall plasma guard. c_plasma = 0.05 would allow 1.67x longer steps at 2.8x
                    the error and about 2 % refusals.
  c_gyro = 0.025    not exercisable here (|B| about 1 T: omega_c dt = 0.0018 at 10 fs; the neighbourhood bound
                    over-predicts the measured Hall gyro value by 10 %); half the Hall guard 0.05, with the same
                    second-order phase error law as the plasma criterion ((omega_c dt)^2/12 <= 5.2e-5 per radian).
  c_disp = 0.1      half the 0.2-cell guard; never binding below about 67 fs (c dt / h_min).
  c_chem = 0.5      half of the per-exposure coverage refusal limit (0.01 unresolved events); no unresolved events
                    occurred in the calibration windows.
  c_drive = 1e-3    relative drive change per step; bound 3.8 ps in the pilot window.
  lte_tol = 1e-8    particle-free only; 14-26 fs in the pilot window (6.2); acts only once the cap exceeds 10 fs.

6.5 Step counts
  With the Hall v3 cap of 10 fs the controller can only shorten steps: in the pilot and qual regimes it chose 10 fs
  on every step (no extra steps, bit-identical to the fixed runs); in the dense-medium regime it took 80 steps
  (about 78 without the calibration landing times) where a fixed 5 fs step needs 100 and a fixed 2.5 fs step 200, and
  avoided the 10 fs step (50 steps, 2.8x the error, about 2 % midpoint refusals at omega_p dt 0.046). Savings against
  10 fs need the larger particle-free cap (section 7).

6.6 Limits found in sealed owners (not step-size limits; the controller cannot remove them)
  - Hall motion work guard: |kinetic change - field work| <= 1e-20 J + 1e-10 (|kinetic change| + |field work|) is
    relative to the net exchange. Every dense run stopped with "canonical motion work mismatch" at the same phase of
    the blob's energy exchange (omega_p t about 2.6-2.8 rad: 0.82-0.88 ps light, 0.40-0.46 ps heavy, 0.23-0.26 ps 4x
    heavy blob), where the net plasma work per step passes through zero; a smaller dt fails slightly earlier (10 fs
    0.88 ps, 5 fs 0.85 ps, 2.5 fs 0.82 ps). The defect is consistent with the push path fixed-point tolerance
    (1e-18 m) times q w |E| per heavy carrier, which does not scale with dt. A tolerance relative to the gross exchange
    (sum over carriers of |w q E_path . D|) is a Hall-owner change needed before any dense-plasma window.
  - Hall absolute Gauss bound: a 10 eV blob (1.6e11) was refused at 0.38 ps ("absolute Gauss constraint disagrees
    with actual declared charge") at both c_plasma 0.035 and 0.045 (the fixed allowance of pf-pilot-1 NOTES).
  - Hall midpoint fixed point: section 6.3 (handled by the retry).
  - Plan section 8 quotes omega_p dt about 0.56 at 1e12 cm^-3 for 10 fs; the Hall operand gives 5.6e-4 at 1e12 cm^-3
    (omega_p = 5.64e10 s^-1). c_plasma = 0.03 at 10 fs corresponds to a node density of 2.8e15 cm^-3 and the 0.1
    guard to 3.1e16 cm^-3.

7. Hall v3.1 linked (builds b4, b5) and what a longer step still needs
  - The build of record links Hall v3.1 through the binding pair (Hall v3.1, native v6). With its default cap the
    coordinator is bit-identical to v6.1 (6.1). --hall-dt-max sets Model(..., dt_max); the coordinator takes the
    Hall cap from Model::dt_max() (bindings_v7.hpp) and StepClock carries the coupled-owner maximum; predicate v7 is
    unchanged. With dt_max = 1e-13 and dt = 10 fs (v31-pilot-fixed10, b4, 2000 steps) the deterministic ledger is
    identical to the Hall v3 run (cal-pilot-fixed10); genesis and step digests differ because v3.1 binds dt_max into
    its configuration identity.
  - A coupled step above 10 fs is still refused by two sealed owners that keep the 1e-14 literal:
      the heat/source owner: electrode-heat v3 kernel TrialFrame::validate (pic-electrode-heat-feedback-20261029,
        mfem_electrode_heat_source_v3/trial_source_binding.hpp:18, "invalid immutable midpoint/common interval"),
        used by source v6.1 on every step: with Hall dt_max = 1e-13 a fixed 40 fs step (v31-pilot-fixed40) was refused
        at step 1, and the adaptive controller (v31-pilot-ctl-lte1e-8, -1e-7, -1e-6) took three 10 fs steps (the LTE
        estimate needs four states) and was refused at step 4 when it grew to 19.3 fs;
      the neutral stream of resident neutral carriers: Hall 8dddec0 kinetic_boundary.cpp:57 ("neutral authoritative
        interval", also in v3.1).
    b5 therefore takes the coordinator's owner maximum as min(Hall Model::dt_max(), source-owner bound 1e-14)
    (bind::coupled_dt_max): a cap above it is refused at configuration (b5-cap-refusal) instead of mid-run, and the
    controller's particle cap stays applied to neutral carriers (neutrals_need_particle_cap).
  - b5 anchors the dt ladder at the particle cap (rung 0 = 1e-14 exactly, negative rungs above it up to the cap), so
    a larger cap keeps 10 fs as a rung (b4, ladder anchored at the cap: 9.6 fs steps with charged carriers under a
    100 fs cap); with cap = particle cap the ladder and every decision above are unchanged.
  - For longer particle-free steps the remaining work is in those owners: a source owner whose trial interval bound
    is configured (as Hall v3.1 did for the field), and a neutral stream with a configured bound. Then cap_s can be
    raised (lte_tol and c_drive govern: 14-26 fs at lte_tol 1e-8 in the pilot window, the bound scales as
    lte_tol^(1/3); the v3.1 dt study gives the field error to 1 ps) while particle_cap_s stays 1e-14 until a carrier
    dt series above 10 fs is calibrated.
  - Hall v3.2 (Gauss mismatch ledger and the native v7 producer forms) with native v7 is planned as one binding pair:
    one HALL_BINDINGS, one NATIVE_BINDINGS and one PAIRS entry in run_v7.py and one adapter section in
    bindings_adapters_v7.hpp (native v7 admission predicate, motion-record adapter); the coordinator step code does not
    name either owner. v3.1's per-step controller receipt can be recorded next to v7's predictions.
  - The absolute Gauss allowance (v3.1 RESULTS section 6: vacuum refusal at 1.236 ns independent of dt; fixed in
    v3.2) and the work guard of 6.6 are not step-size limits; they bound long and dense windows regardless of dt.
