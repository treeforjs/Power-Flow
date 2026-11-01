# Coupled PIC water-breakdown: handoff

Last updated 2026-10-31 23:30 MDT. Keep this file current at every milestone so
the work can be picked up from it alone. Snapshots of the in-progress workers'
own notes are in [`docs/coupled_pic_handoff/`](coupled_pic_handoff/) (copied at
each update); the live versions are the `HANDOFF.txt` files named below.

## Takeover and performance qualification (2026-10-31 23:47 MDT)

The original audit corrections are committed and pushed through `e168eb8`.
Continuation now starts from those corrected components, with separately
versioned successors. Existing external workers remain active; do not edit
or stop their builds or numerical jobs. New numerical work goes through
SimQueue with explicit CPU, memory and GPU claims. Functional checks can
borrow a scheduler-lent idle core and freeze when its owner needs it.
Performance claims require an uncontended physical core, including its SMT
sibling; a core lease alone does not exclude the current external SMT builds.

Working successors, all unqualified until their receipts pass:

- `pic-water-common-commit-v8.0.2-cpu-20261031`: exact product-error fast path
  with environment/range fallback, and immutable cross-section query caching.
  Compare against the actual v8.0.1 c2 binary, including deterministic physical
  histories and continuation. Study timings are not release timings.
- `native-full-time-v8.2-birth-memo-20261101` under the native runtime:
  reuse a validated host birth-scope proof. Preserve the audited identity and
  chemistry guards. Compare exact physical records and ledgers from both real
  archives; source-identity hashes necessarily differ.
- `swarm-validation-audit-v1-20261101`: reject incomplete/duplicate/nonfinite
  swarm output, account for the writer's exact serialized field labels,
  bind outputs to their tested source and prevent incomplete seals.
- `gas-target-owner-v1.1-guards-20261101`: refuse nonrealizable feedback instead
  of clamping temperature, make caller updates atomic, and require two ordered
  chemistry stages. This is a guarded successor of unfinished gas-owner work,
  not a qualified coupled energy/momentum owner.

Current live chemistry v3.2 has overlapping output writers and lacks immutable
completion provenance. Its v3.2.1 adoption has not yet regenerated and qualified
all tables. Existing live validation/seal scripts can pass incomplete evidence;
no production promotion follows from those files. Changed swarm comparison
criteria and model-derived cross-section extensions need separate scientific
assessment. Technical completeness does not establish model validation.

The isolated swarm integrity correction is now qualified for its technical
scope: frozen SimQueue job 312 passed 31 synthetic regression tests. Publication
is `artifacts/swarm-validation-audit-v1-20261101`; it includes the exact-base
adoption patch and bound evidence. It has not been applied to active writers,
and explicitly does not qualify any physics model or production run.

The next integration milestone is a conservative all-owner pilot after gas
recoil/target feedback, current chemistry binding, thermal/field/source coupling
and restart have passed their component and coupled conservation checks.
A 20 ns field-only checkpoint cannot supply omitted early plasma history.
The old electron-only production run remains paused. Full breakdown, a full
continuous 60 ns history and physical convergence are not yet qualified.

## What this is

A coupled particle-in-cell simulation of water-vapour breakdown at the foil
electrodes, driven by the C19 (Mykonos LTD, shot 17133) current. Independent
owners exchange sealed, digest-bound state through a coordinator:

| Owner | Role |
|---|---|
| Hall corrective | electromagnetic field (E, B, J), Gauss ledger, field-only pre-run |
| Native reaction transactions | carriers, collisions (MCC), births/sinks, motion record, restore |
| Water chemistry | e/ion + H2O cross-section tables, coverage policy (refuse vs bounded) |
| Source owner | electrode surface physics: emission, desorption, wall response |
| SS304 phase owner | electrode heating, melt, vapour |
| C19 circuit | drive current / port coupling |
| Coordinator | step control (adaptive dt), binding, restart, launch (`pf-coupled-v8`) |

Every component is built and qualified in its own versioned directory, sealed,
and then copied (filtered) into `artifacts/<component>-<version>-<date>/` here.

## Where things live

- Builds: `/home/treefor/projects/runtimes/<component>-<version>-<date>/` (not in git;
  each has README.txt and usually SEAL.json or SOURCE-SHA256.json, RESULTS.txt,
  INTEGRATION.txt and a HANDOFF.txt while in progress).
  Native lives under `runtimes/pic-reaction-transactions-20261029/native-full-time-candidateNNN/`.
- Run outputs and receipts: `/home/treefor/sim-workspace/powerflow/outputs/<component>-<version>-<date>/`
  (same directory as `/home/treefor/projects/.workspace/outputs/`).
- Integration backlog: `/home/treefor/sim-workspace/powerflow/outputs/integration-backlog.txt`
  (copy: [`coupled_pic_handoff/integration-backlog.txt`](coupled_pic_handoff/integration-backlog.txt)).
- Repo copies: `artifacts/` here. Water-chemistry tables derived from LXCat stay
  local (LXCat terms); only the specs, digests and generators are committed.
- Tools: `tools/pf-launch` (per-run core/GPU selection with site files),
  `tools/simqueue` (node scheduler source), `tools/pilot-render`.

## Status (2026-10-31 23:30)

| Component | Sealed / committed | In progress | Next |
|---|---|---|---|
| Hall | v3.4 (61b04a1): host work removed (device deposit support, preflight, census; exact-sum residual; block-form digests). Hall-only step 7-19x faster than v3.3 at 1e5-1e6 carriers; per-kernel FP64/bandwidth/latency split in RESULTS.txt sections 12-13 | v3.5 (`runtimes/hall-corrective-v3.5-20261101`, HANDOFF.txt): `restore_as_genesis` first, then GPU node-walk parallelization (up to ~70 % of the 200 um device step) and radix-sort passes | graded 35-45 mm transverse mesh, circuit port-current mode, SS304 v1.3 thermal owner, source v7 wall heat, merge/split relocation deposit, GMG option, open issues 10.1-10.6 (v3.5 or v3.6) |
| Native | v8 audit correction (f09bda2), source 5f7d171c: conservation, checkpoint identity, restore and device lifetime corrections qualified | v8 (`candidate008`, HANDOFF.txt): independent device-birth/sink and Stage B development; corrective adoption patch supplied | adopt audit corrections before promotion; swap chemistry binding to v3.2; Raikou kinetic-core adoption |
| Water chemistry | v3.1.1 (ea2b9ee): thermal coverage, positive parcels, strict unknown-channel policy and debit corrections | v3.2 (`runtimes/pic-water-chemistry-v3.2-20261031`, HANDOFF.txt): model-derived extensions and swarm qualification | adopt the checked v3.2 correction patch, then qualify and bind into native v8 |
| Coordinator | v8.0.1 audit release: 53 checks, 12 continuation cases and 5 integrity faults passed; corrected native dependency pinned | - | stage 2 after Hall v3.5 + native v8: bind all owners, direct native restore, qualify any omitted early plasma history |
| Source owner | v7 (faa7274) component | - | integrate at coordinator v8 stage 2 |
| SS304 phase owner | v1.3 (d39e33d): transactional trial, pressure and accumulator corrections | - | thermal owner in Hall v3.5/v3.6 |
| C19 circuit | v1 (2feff69); Xyce 7.10 source verified, not built | - | port-current mode in Hall v3.5/v3.6 |
| SimQueue | 0.6.2 live (78767d7): fair GPU hand-off, reservation line, published state `pic-profile/simq-gpu-line.json`, `gpu-lease` helper; rollback command in `tools/simqueue/release.json` | - | 0.6.3 at the CPU swap: zero-core borrow-only `mfem-pic` partition with lender list, 16-core node layout |
| CPU performance | - | study (`runtimes/cpu-isa-study-20261031`, HANDOFF.txt): x86-64-v3/v4/znver5/-O3 builds must give bit-identical digests; host profiling; proven patches collected in PATCHES.md for the next versions; per-site ISA flavor via pf-launch site files | apply patches in Hall v3.6 / native v8.x / coordinator stage 2 |

The 20 ns field-only pre-run inside the coordinator model matched the Hall v3.3
20 ns checkpoint word for word (194,026 steps, 631 s). Starting the coupled run
from it also omits the earlier particle and source history. Hall v3.5
`restore_as_genesis` addresses a revision/accepted-step binding problem; that
software fix alone does not qualify a cold plasma start at 20 ns.

## Audit corrections (2026-10-31)

Corrective releases below supersede the audited baselines for their stated
component cases. They do not qualify a full plasma-breakdown run.
SS304 phase v1.3 is qualified and pushed (d39e33d): 200204 CPU checks and
73 GPU checks passed, with three historical defects separately reproduced.
Water chemistry v3.1.1 is qualified and pushed (ea2b9ee): 11 regressions,
CPU process/event controls, ASan/UBSan and two-process CPU/GPU bit identity
passed; all 623 frozen data bindings are unchanged. Its exact-base v3.2 adoption
patch is under `artifacts/water-chemistry-v3.1.1-20261031/adoption/v3.2/`.
Native candidate008-auditfix is qualified and pushed (f09bda2), source
`5f7d171c62708cf3e2e03f0fd843e7c906ed65a5db7d5955624773dee3516f2b`, archive
`098c1aaa14a533877ba4474f1ffae432994329db809b12b2f1391c6e2a65f9d7`.
Its full CPU suite and GPU lifetime/flow/chemistry controls passed, including
long coverage labels and v8.1 journal restore. Coordinator v8.0.1 pins that
archive and passed frozen job 293: 53 checks, 12 continuation cases and five
integrity fault cases. Its default binding remains Hall v3.3, parent ionization
and v5 native configuration; nonempty coordinator restore still uses replay.

Use `artifacts/coordinator-v8.0.1-20261031/RESULTS.txt` and
`docs/audit-fixes-20261031.txt` for the exact evidence and limitations. The
historical sealed native v7.1 device-verification path remains affected by the
audited lifetime defect. Independently edited active runtimes are unchanged:
fold the checked adoption patches into native v8/chemistry v3.2 and requalify
before promotion. Old chemistry checkpoints cannot silently acquire missing
parcel, sink or reachable-support identity.

## Open model choices

- Coordinator v8 stage 2 production presets: enabling source v7 surface physics
  in production, the melt inertia model, sub-cycling K.
- Mixed precision for the field solve (FP32 factor + FP64 refinement): faster on
  the RTX 3080 Ti but changes results within tolerance; bring numbers first.
- GPU purchase deferred (prices); Hall's per-kernel split says node walks are
  latency-bound, so software comes first.

## Decisions in force

- Pilot path first: a 1-2 ns coupled electron + ion + H2O-ionization pilot with
  conservation ledgers once the all-owner step works; long-history restart work after.
- Any rank count N, configurable at launch, results independent of N (keyed RNG,
  ordered reductions, rank-independent checkpoints). Portable site files.
- Bit-identical results are non-negotiable for every optimization (CPU ISA
  flavors, GPU kernel changes): an optimization is kept only if digests match.
- Adaptive time step from the physics guards replaces the fixed 10 fs; the first
  step after a start stays at 10 fs (source v7 needs two rate samples). Particle cap 10 fs by default.
- Neutral births use the quantized electron route.
- Parent ionization: BEB to 10 keV, provisional RBEB above until measured data.
- Hall v3.2 Gauss-ledger guard; native v7 exact superaccumulator sums and
  block-Merkle motion digest; Hall v3.4 D1-D3 (exact-sum residual, block-form
  digests, coordinator takes Hall's digest in the v7.1 lane).
- Native v8: chemistry integration choices 1-7; coverage policy refuses
  unbounded/bounded-only exposure above max(1 event, 1 % of applied events);
  C3 max_events 64, max_candidates 4096; nonempty restore included.
- Chemistry v3.1: recombination 2-30 eV from digitized Rosen 2000 / Jensen 2000
  figures, bounded above; slow-ion join rules; H- detachment zero below 0.75 eV.
- Chemistry v3.2 (data gaps never stop runs; estimates are sampled, tagged
  model-derived with uncertainty, and the geometric/total bounds stay documented
  with their ratios):
  - H+ + H2O: Wedlund 2019 capture/ionization fits, Gobet 2004 product channels,
    Rudd 1992 secondary-electron spectra; momentum transfer ZBL (PSTAR for H+).
  - Singly charged dressed ions (O+, OH+, H2O+, H3O+, H2+): ionization = He+
    sigma11 at equal speed (validated against the Juhasz 2019 O+ 6.6 keV point);
    capture = proton sigma10 at equal speed; O++ from He2+ data; anions
    q^2 x proton; projectile breakup/detachment from the free-collision estimate
    (1.5-5x high against He+ and H- checks).
  - H3O+ capture, O- rearrangement, H2O+ and H- proton transfer at keV: sampled
    estimates (energy-defect-scaled / equal-speed forms).
  - Electron rotational J=0->2, 0->3 above 5 eV: J=0->1 energy dependence scaled
    at 5 eV (factor 3); ion-pair formation above its windows: Munoz total-inelastic
    shape (factor 10); dissociative attachment zero above its resonance tables.
    Electron swarm re-validated against measured drift/diffusion at v3 tolerances.
  - Ions 3 eV-1 keV: log-log bridges between low-energy data and the 1 keV values
    (factor 3-10); H- detachment rises from zero at threshold.
- Coordinator v8: restart by verified replay now, real native restore in v8;
  emission/desorption error targets 1e-3; Hall v3.5 `restore_as_genesis` for the 20 ns start.
- Shared with Raikou: nuclei ABI v4 (16 slots: 0 H, 1 O, 2 Fe, 3 Cr, 4 Ni, 5 Mn,
  6 Si, 7 C, 8 N, 9-15 spare), 134-digit superaccumulator as the standard exact
  sum, Philox purpose registry (arrival 0x01-0x03, chemistry 0x30-0x3F, test
  0xF1-0xF2, Raikou P8 0x40-0x4F; salted source v7 wall 0x21/0x22, merge/split).
- C19: 5-stage LTD, 36 bricks/stage, 40 nF per capacitor, +/-80 kV charge,
  36 nH internal; feed ~40 nH, target ~25 nH. Coupled runs need a >= 35-45 mm
  graded transverse box (the 9 mm box underestimates port inductance ~32 %).

## Exact commands

All builds and runs use the frozen runtime and the core leases; outputs go to
`/home/treefor/sim-workspace/powerflow/outputs/<component>/`.

- Scheduler: `simq status`, `simq partitions`, `simq show ID`;
  `simq submit -p powerflow -c N [--gpus 1 --gpu-mem 2G] --mem X -t T -n NAME -C DIR --shell "CMD"`.
- Direct GPU use: `G=/home/treefor/projects/runtimes/simqueue-0.6.2-20261031/bin/gpu-lease;
  $G --core any -- CMD` (or `--core both`, `--core K`, or GPU only `$G -- CMD`); see `USAGE-GPU.md` there.
- Core lease for CPU work: `flock <pic-profile>/testing-core7.lock taskset -c 7 nice -n 15 CMD`
  (core 6 uses `testing-core.lock`), or `runtimes/core-lease-20261030/any-core.sh --wait CMD`.
- Hall (in a version directory): `python3 run_hall_v34.py build|run|gpu ...` (v3.5: `run_hall_v35.py`);
  acceptance `run_acceptance.sh`; coordinator parity `run_coordinator_parity.py`; 20 ns pre-run
  `run_prerun.sh` (SimQueue, exclusive GPU ~11 min, notify the other sessions first).
- Native v8 (in `candidate008`): `scripts/lease.sh env CUDA_VISIBLE_DEVICES= python3 build-components.py normal`;
  CPU suite `scripts/cpu_suite.sh OUTDIR`; timing table `python3 scripts/timing_table.py HOSTDIR DEVDIR`;
  compatibility `scripts/compat.py`.
- Water chemistry (in its version directory): `taskset -c 6,7 nice -n 15 python3 -I run_checks.py LABEL
  [inputs tables build core sanitized gpu swarm]`; table rebuild: `tools/make_production_spec.py`,
  `tools/build_tables.py`, `tools/report_processes.py`, `tools/make_coverage_policy.py` (see its README).
- Corrected coordinator (`pic-water-common-commit-v8.0.1-20261031`): `python3 run_v8.py build LABEL [--hall v3.3] [--native v8-auditfix] | run | baseline`;
  launch `tools/pf-coupled-v8` with pf-launch options (`--site --cores/--cpus --gpus --dry-run --receipt`)
  plus `--hall-checkpoint DIR`, `--resume`, `--mesh`, `--dt-cap`, `--preset production`, `--hall-owner auto`,
  `--window/--steps/--end-time`, `--checkpoint-every` (README section 6).
- Per-site ISA flags (all components): `PF_BUILD_ISA_FLAGS="-march=..."` (empty by default; a build with
  flags set needs its own label and a bit-identity check against the generic build).
- Repository tests: `/home/treefor/projects/runtimes/arch-base-20261028/run-frozen env PYTHONPATH=$PWD
  /home/treefor/projects/Power-Flow/.venv/bin/python -m pytest -q tests` from the repo root, under a core lease.

## Machine and sharing rules (zenobia)

- Ryzen 7 9850X3D (8 cores; logical 0-7 cores, 8-15 SMT), RTX 3080 Ti 12 GB, 30 GB RAM
  (~12 GB free under load). A 16-core Ryzen 9 9950X3D2 is due around 2026-11-05;
  afterwards repartition and requalify by fingerprint comparison.
- Cores (development split): 0-3 Spheral, 4-5 raikou-mfem, 6-7 Power-Flow; after the 16-core swap Spheral gets one
  chiplet (8 cores), raikou-mfem and Power-Flow 4 each on the other, shared mfem-pic work borrows idle cores there.
  A properly optimized run that needs all ranks may use all cores; during a Spheral full-run window Power-Flow keeps
  only logical CPU 15 (nice 19) for one GPU host thread. Lending Power-Flow cores needs your direct OK.
  The GPU is shared with raikou-mfem and the Raikou kinetic-core work.
- Always pass `--mem`; never take the GPU lock files with raw flock (use the scheduler or `gpu-lease`).
- Lock order: GPU locks first (`testing-gpu.lock`, then `testing-gpu-<UUID>.lock`),
  then the core lease. Lock files: `/home/treefor/sim-workspace/powerflow/outputs/pic-profile/`.
- Notify the other project sessions before any exclusive or >30 min GPU run, and when it ends.

## Resuming

1. Read this file, the snapshots in `docs/coupled_pic_handoff/`, and the integration backlog.
2. `simq status` and `simq partitions` for running work.
3. For an in-progress component, read its runtime `HANDOFF.txt` and the newest
   receipts under its outputs directory; rebuild and rerun tests before trusting
   anything that was mid-edit.
4. Check `git status` and `git log` before committing: other work (audits) may
   have changed files; do not revert changes you did not make.
5. Commit sealed components here (filtered copies, no large checkpoints) and push
   to `master`. Work is committed directly to `master`; no long-lived side branches.

## Repository history

On 2026-10-31 every branch was merged into `master` and then deleted. The
superseded `research/flash-publication-foundation` docs were recorded as merged
with the current text kept; merged files received the workspace-path
normalization (`artifacts/path-normalization-20261031`). Last commit of each
deleted branch, for restoring one (`git branch <name> <commit>`):

| Branch | Last commit |
|---|---|
| `electrode-heat-feedback-20261029` | `373c91a83c9c` |
| `electrode-heat-source-v3-20261029` | `51a867d96b78` |
| `electrode-heat-source-v4-20261030` | `c65ac51874f0` |
| `neutral-spectral-observation-wip-20261028` | `1da3b3a6653d` |
| `research/flash-publication-foundation` | `18d155b919e7` |
| `solid-material-feedback-pilot-20261029` | `19eee7e4d15e` |
| `ss304-melt-ablation-proposal-v1-20261030` | `1c116564991e` |
| `ss304-native-metal-transfer-20261029` | `ef5e28f47f67` |
| `ss304-spatial-metal-transfer-20261029` | `f2d3f93e494d` |
| `ss304-thermal-expansion-20261029` | `7ff9faad2ff3` |
| `ss304-thermal-properties-20261029` | `1445d9cda4fd` |
| `ss304-thermal-properties-v2-20261029` | `ca74b667f19a` |
| `surface-emission-20261029` | `52be614fa3dc` |

## Known issues

- Dissertation evidence serialization: retained the manifest-bound CRLF
  records and added `.gitattributes` preservation rules. The complete suite
  then exposed three CSVs whose recovery generators and own metadata require
  LF, while the aggregate pilot manifest recorded CRLF copies. Those three
  now use the generators' exact LF bytes and original per-artifact hashes;
  the aggregate manifest references were corrected. No scientific values,
  per-artifact expected hashes or tests changed. The byte-level change record
  is `artifacts/source-audit-20261031/eol-repair.json`. Full-suite SimQueue job
  255: 3650 passed, 39 skipped; the existing tests and expected hashes are intact.
