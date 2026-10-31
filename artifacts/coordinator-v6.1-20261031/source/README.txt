Coordinator v6.1 -- coordinator v6 with less host work per step, same accepted results -- 2026-10-31

Coordinator v6 (pic-water-common-commit-v6-20261031, build q3, binary eb1936f1...) is sealed and read in place.
v6.1 is a versioned copy (namespace raikou::native_water::v6_1) that uses heat/source owner v6.1
(pic-electrode-heat-source-v6.1-20261031, API-v6.1.txt) in place of source v6. Every other owner, the step
order, every guard, every refusal point and every identity text are the v6 ones. Nothing physical changed.

Layout
  source/identity_v6_1.hpp     v6 identities and setup; source owner v6.1, binding hall_source_binding_v6_1
  source/joint_step_v6_1.hpp   observe, prepare_v6, can_publish_v6, try_publish_v6 (changes below)
  source/genesis_v6_1.hpp      v6 genesis (unchanged apart from the namespace)
  tests/driver_v6_1.cpp        modes run | checks | vacuum (v6) and profile (sealed v6 vs v6.1 on a live state)
  run_v6_1.py                  build | run | baseline (the sealed v6 q3 binary, sha256-verified, same admission)
  tools/compare_v6_v61.py      identity (file bytes + every deterministic ledger column by name + frames),
                               timings, checks; tools/write_summary_v61.py; v6 tools copied unchanged

Changes relative to v6 (all inside prepare_v6 / can_publish_v6 / observe; publication unchanged)
  1. Source field. v6 copied the Hall mesh and re-hashed it (about 9 MB of encoding, 10 ms) in every source Field
     and in every Field::matches call (each slice trial prepare and preflight). v6.1 builds the Field from the
     immutable mesh identity of the live Hall mesh (s6::MeshIdentity::of: one snapshot hashed once, reused while
     the live mesh is exactly equal to it); matches() compares frame, mesh (exactly) and E/B (bit for bit). The
     one-sided extraction field at each patch is sampled from the source Field's own trial field (v6 built a
     second AxialTrialField from a mesh copy with the same values).
  2. Emission exposures. Source v6.1 prepare_patch with the exact fast spectrum (bit-identical to the sealed
     surface_emission::evaluate; API-v6.1.txt 2.4). The prepared flux of each exposure is kept and reused to
     sample the electron births of that stream (v6 prepared the patch again on a third trial field).
  3. observe(). The Hall model's motion mesh digest and the frame digest are recomputed only for another model
     object or a changed native options byte image (the model is an immutable shared object); every comparison
     of observe() is still made. Source owner state digests are cached per owner generation and thermal clock
     (source v6.1).
  4. Sub-phase timings: subphase-timings.bin / subphase-schema.json (source field digest and build, extracting
     and retarding exposures, can_publish parts: native, chemistry, observe, predicate, sources, arrival, Hall).

Emission bound (task item 2) -- decision
  A bounded "negligible exposure" entry in place of a full evaluation was not adopted: arrival v1.1 hashes the
  exact expected count N, its log and N/w0 into the step plan and accumulates N exactly (expected_sum) and
  f(1-f) (variance_sum) in its state digest, so the arrival state, the source transition seal and every later
  step digest (and the ledger columns arr_expected_electrons_step, arr_cum_expected_electrons,
  arr_max_log10_expected_electrons_patch) would change although the decision (no births) would not. The
  bound is used inside the kernel quadrature instead, where it is exact: an energy node is skipped only when a
  rigorous upper bound of its log weight (monotone correctly rounded operations of the kernel's own action)
  proves that its weight underflows to +0. 53 % of the extracting-field nodes are skipped on the pilot.

Outputs that are bit-identical to v6 (same build inputs, same options)
  genesis.json, ledger-schema.json, deterministic.bin (all deterministic ledger columns), every deterministic
  column of ledger.bin, step-digests.txt (native state prefix, Hall field/heat/kinetic-ledger digests, accepted
  source seal, chemistry RNG digest, arrival state), frames-layout.json, frames.bin, frames.jsonl; every owner
  state (the checks' owner snapshots equal the v6 T1 reference); source owner configuration, state, trial,
  candidate, prepared-patch, field and mesh digests and checkpoint bytes (source v6.1 tests).
  Different by design: timing columns (t_*, rss_kib) of ledger.bin, run-summary.json (wall times; added fields
  coordinator, source_owner, emission_statistics, subphase_median), receipts, refusal message prefix
  "coordinator v6.1:"; new files subphase-timings.bin, subphase-schema.json.

Build and run
  core-lease-20261030/any-core.sh --wait /usr/bin/python -I run_v6_1.py build LABEL
  core-lease-20261030/any-core.sh --wait /usr/bin/python -I run_v6_1.py run LABEL run NAME --fixture qual --steps 20
  flock GPU locks + any-core.sh ... run_v6_1.py run LABEL run NAME --fixture pilot --field gpu --push gpu --parity-every K
  run_v6_1.py run LABEL checks NAME | profile NAME --fixture pilot --steps N
  run_v6_1.py baseline run NAME [options]       (sealed v6 q3 binary, for before/after)
  Pilot (SimQueue): as v6 with run_v6_1.py run a2 run <out>/run ... (same options).
  Build of record: a2 (binary sha256 in outputs/pic-water-common-commit-v6.1-20261031/build-a2/build-receipt.json).
  Build a1 has the same coordinator, source owner and kernel sources; a2 only adds sealed-owner timings to the
  profile mode of the driver. The a1 runs (a1-*) are kept as a first, identical evidence set.

Qualification of record (outputs/pic-water-common-commit-v6.1-20261031/, qualification/QUALIFICATION-SUMMARY.json)
  Bit identity against the v6 records (compare_v6_v61.py identity: files, 113 deterministic ledger columns by
  name, frames):
  (builds a1 and a2, runs a1-* and a2-*):
    qual CPU 20 steps, two fresh processes (vs q3-qual-cpu-a/b)                           identical
    qual GPU field + GPU push 20 steps (vs q3-qual-gpu-a)                                 identical
    pilot CPU 20 steps, pilot GPU 20 steps (vs q3-pilot-cpu-20, q3-pilot-gpu-20)          identical
    pilot GPU 2000 steps (vs the pf-pilot-1 record, rows 1-2000; frames at 0/1000/2000)  identical
    pilot CPU 2000 steps (vs base-v6-pilot-cpu-2000, the sealed v6 binary run alongside)  identical
  Checks (T1-T4, T6, spot): 31/31, the same cases, outcomes and refusal points as q3-checks-cpu.
  Source v6.1: tests 11/11, the sealed v6 suite 9/9 on v6.1 (source-v6.1/s4).

Per-phase timings, median ms per step (20-step windows, build a2; v6 = sealed q3 binary run on the same cores)
                              qual CPU        qual GPU        pilot CPU       pilot GPU
                              v6     v6.1     v6     v6.1     v6     v6.1     v6     v6.1
  source field                11.2   1.6      11.1   1.7      11.0   1.7      11.0   1.7
  emission exposures          39.2   10.1     39.4   10.2     39.3   10.0     39.5   10.1
  emission births             22.1   0.03     22.2   0.03     0.004  0.002    0.004  0.003
  source trials               77.3   2.1      78.2   2.2      0.002  0.002    0.002  0.002
  Hall proposal               49.7   50.1     55.0   53.9     41.3   41.0     42.0   41.9
  can_publish (preflight)     91.1   15.8     91.7   16.2     10.9   10.9     11.0   11.0
  prepare total               223.9  89.3     229.9  95.2     101.0  61.9     102.4  63.5
  step total                  318.1  105.3    321.8  111.8    112.1  72.8     113.5  74.6
  2000-step pilot windows: CPU backends 112.3 (v6, run alongside) -> 73.4 ms/step; GPU backends 110.6 (pf-pilot-1
  record, rows 1-2000) -> 72.8 ms/step. Emission exposures 39.6 -> 10.3 ms, source field 11.0 -> 1.7 ms; 53 % of
  the extracting-field quadrature nodes (2.6e8 of 4.9e8) were skipped as exact zeros; no sealed fallback.

Remaining host costs (pilot, GPU backends, v6.1, about 73 ms/step; a2-prof-pilot/profile.json for the pieces)
  Hall proposal 41.9 ms: field prepare 28.7 (moving to the GPU in Hall v3.1), motion 3.2, heat 2.3, census
    admission 1.7, seals 1.5, solve 1.2.
  Preflight 11.0 ms: Hall preflight_final 8.8 (sealed Hall v3: its independent recomputation of the accepted
    field-state digest, 0.95 ms although equal to its cached value, of the candidate digests, the midpoint arrays
    and their total-field digest, 0.58 ms, and the motion record), native v6 predicate 1.1, observe 0.8, chemistry
    can_publish 0.3. Caching Hall's accepted-state digests per generation is a Hall-owner change.
  Emission exposures 10.1 ms: 16 extracting patches x 0.6 ms (about 7000 contributing quadrature nodes each; the
    per-node libm calls log/exp/log1p dominate), retarding patches 0.4 ms in total.
  observe 1.5 ms in prepare and 0.8 ms in preflight: the sealed chemistry owner recomputes its accepted RNG digest
    on every call (0.71 ms), and prepare builds the chemistry checkpoint JSON to read next_interval (0.62 ms).
  The total-field digest of the midpoint arrays is computed three times per step (coordinator binding check,
    can_publish, Hall preflight; 0.58 ms each) and that of the accepted arrays once (source frame name).
  Chemistry exposures 2.5 + 2.7 ms, seal 1.2, native motion admission 0.6, source field build 1.1 (sealed v3
    trial-field copy and fingerprint of E/B, content digest, exact live-mesh comparison 0.2).
