# Coupled PIC water-breakdown: handoff

Last updated 2026-10-31 21:50 MDT. Keep this file current at every milestone so
the work can be picked up from it alone.

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
  each has README.txt and usually SEAL.json, RESULTS.txt, INTEGRATION.txt or a HANDOVER file).
  Native lives under `runtimes/pic-reaction-transactions-20261029/native-full-time-candidateNNN/`.
- Run outputs and receipts: `/home/treefor/sim-workspace/powerflow/outputs/<component>-<version>-<date>/`.
- Integration backlog (next versions, kept current):
  `/home/treefor/sim-workspace/powerflow/outputs/integration-backlog.txt`.
- Repo copies: `artifacts/` here. Water-chemistry tables derived from LXCat stay
  local (LXCat terms); only the specs, digests and generators are committed.
- Tools: `tools/pf-launch` (per-run core/GPU selection with site files),
  `tools/simqueue` (node scheduler source), `tools/pilot-render`.

## Status (2026-10-31 evening)

| Component | Sealed / committed | In progress | Next |
|---|---|---|---|
| Hall | v3.3 (111953c): exact-form magnetic update, 20 ns pre-run | v3.4: remove remaining host work (device deposit support, device preflight/census, exact-sum residual, block-form digests for native v7.1) | v3.5: `restore_as_genesis` (start the coupled run from the 20 ns vacuum checkpoint), graded 35-45 mm transverse mesh, circuit port-current mode, SS304 v1.2 thermal owner, source v7 wall heat, merge/split relocation deposit, GMG option |
| Native | v7.1 (2c37382): device admission, copy-on-write state | v8 (candidate008): chemistry v3.1 binding, device chemistry, nonempty restore, wall births, metal species, merge/split; nuclei ABI v4 | Raikou kinetic-core adoption |
| Water chemistry | v3.1 (848e3e9): recombination and slow-ion coverage | v3.2: keV ion/negative-ion + H2O coverage (data where available, physics bounds elsewhere) | swap into native v8 |
| Coordinator | v8 stage 1 (93e3936): source v7 wired, raised dt caps, restart, launch from checkpoint | - | stage 2 after Hall v3.5 + native v8: bind all owners, production presets, direct native restore, start from 20 ns |
| Source owner | v7 (faa7274) component | - | integrate at coordinator v8 stage 2 |
| SS304 phase owner | v1.2 (66fa599) | - | thermal owner in Hall v3.5 |
| C19 circuit | v1 (2feff69); Xyce 7.10 source verified, not built | - | port-current mode in Hall v3.5 |
| SimQueue | 0.6.1 live (d819b50) | 0.6.2: fair GPU handoff (reservation line, published state, direct-use helper), live upgrade with backup/rollback | partitions for the 16-core CPU |

The 20 ns field-only pre-run inside the coordinator model matched the Hall v3.3
20 ns checkpoint word for word (194,026 steps). Starting the coupled run from it
is blocked only by the revision/accepted-step binding between Hall and native;
Hall v3.5 `restore_as_genesis` resolves it without changing any check.

## Decisions in force

- Pilot path first: a 1-2 ns coupled electron + ion + H2O-ionization pilot with
  conservation ledgers once the all-owner step works; long-history restart work after.
- Any rank count N, configurable at launch, results independent of N (keyed RNG,
  ordered reductions, rank-independent checkpoints). Portable site files.
- Adaptive time step from the physics guards replaces the fixed 10 fs; the first
  step after a start stays at 10 fs (source v7 needs two rate samples). Particle cap 10 fs by default.
- Neutral births use the quantized electron route.
- Parent ionization: BEB to 10 keV, provisional RBEB above until measured data.
- Hall v3.2 Gauss-ledger guard; native v7 exact superaccumulator sums and
  block-Merkle motion digest; Hall v3.4 D1-D3 (exact-sum residual, block-form
  digests, coordinator takes Hall's digest in the v7.1 lane).
- Native v8: chemistry integration choices 1-7; coverage policy refuses
  unbounded regions when exposure exceeds max(1 event, 1 % of applied events);
  C3 max_events 64, max_candidates 4096; nonempty restore included.
- Chemistry v3.1: recombination 2-30 eV from digitized Rosen 2000 / Jensen 2000
  figures, bounded above; slow-ion join rules; H- detachment zero below 0.75 eV.
- Chemistry v3.2: close the keV ion + H2O gap with literature data or physics bounds.
- Coordinator v8: restart by verified replay now, real native restore in v8;
  emission/desorption error targets 1e-3; Hall v3.5 `restore_as_genesis` for the 20 ns start.
- Shared with Raikou: nuclei ABI v4 (16 slots: 0 H, 1 O, 2 Fe, 3 Cr, 4 Ni, 5 Mn,
  6 Si, 7 C, 8 N, 9-15 spare), 134-digit superaccumulator as the standard exact
  sum, Philox purpose registry (arrival 0x01-0x03, chemistry 0x30-0x3F, test
  0xF1-0xF2, Raikou P8 0x40-0x4F; salted source v7 wall 0x21/0x22, merge/split).
- C19: 5-stage LTD, 36 bricks/stage, 40 nF per capacitor, +/-80 kV charge,
  36 nH internal; feed ~40 nH, target ~25 nH. Coupled runs need a >= 35-45 mm
  graded transverse box (the 9 mm box underestimates port inductance ~32 %).

Open decisions for coordinator v8 stage 2: enabling source v7 surface physics in
production, the melt inertia model, sub-cycling K.

## Machine and sharing rules (zenobia)

- Ryzen 7 9850X3D (8 cores; logical 0-7 cores, 8-15 SMT), RTX 3080 Ti 12 GB, 30 GB RAM
  (~12 GB free under load). A 16-core Ryzen 9 9950X3D2 is due around 2026-11-05;
  afterwards repartition and requalify by fingerprint comparison.
- Cores: 0-3 Spheral, 4-5 raikou-mfem, 6-7 Power-Flow. The GPU is shared with
  raikou-mfem and the Raikou kinetic-core work.
- Use SimQueue: `simq submit -p powerflow -c N [--gpus 1 --gpu-mem 2G] --mem X -t T -- cmd`;
  `--borrow` for idle cores of other partitions. Always pass `--mem`.
- Lock order: GPU locks first (`testing-gpu.lock`, then `testing-gpu-<UUID>.lock`),
  then the core lease. Lock files: `/home/treefor/sim-workspace/powerflow/outputs/pic-profile/`.
- Notify the other project sessions before any exclusive or >30 min GPU run, and when it ends.

## Resuming

1. Read this file and the integration backlog.
2. `simq status` and `simq partitions` for running work.
3. For an in-progress component, read its runtime README/HANDOVER and the newest
   receipts under its outputs directory; rebuild and rerun tests before trusting
   anything that was mid-edit.
4. Commit sealed components here (filtered copies, no large checkpoints) and push
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

- 14 tests in `tests/test_pilot_comparison_bridge.py` and
  `tests/test_pilot_comparison_evidence.py` fail on Linux checkouts:
  `data/dissertation/pilot_comparison_manifest.json` records the sha256 of
  `data/dissertation/figure_4_20_profiles.csv` with CRLF line endings
  (53a0ac56...), while the repository stores LF bytes (551a6cbb...). Fails
  identically before and after the 2026-10-31 merge; the rest of the suite
  passes (3636 passed, 39 skipped).
