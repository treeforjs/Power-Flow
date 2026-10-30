# V4 stable handoff — sealed heat/source v2 and signed mechanics

Implemented, qualified, preserved and pushed on `ss304-thermal-expansion-20261029`.
Source: `85ba2a1640a56dfd11d9534ec3e5de881af3b0ec`. Compact results: `6fe6a0679d6f086533f9e149a0798460432a87c0`.
Final results and verified remote: `7ff9faad2ff3cd46d02b7d67b415505c5c44a775`. Working tree is clean.

All 62 source changes are new files. V3, the baseline SS304 record and the heat
owner's live source are unchanged. The new package preserves all 25 sealed v2
snapshot files and all 12 working kernel/material/binding inputs. The exact
source kernel is compiled once against the frozen MFEM archive.

## Qualified behavior

Signed mechanics, native conduction, actual electron birth/continuation and
actual adsorbed-water transfers close local and global discrete energy. The
nineteen thermal cumulative ledgers preserve electron kinetic/reference/image
and field-work ownership. Only electron wall heat enters the thermal equation;
gas H2O reference remains zero and desorption raises the negative adsorbed
reservoir. Joule and authored heat are zero.

The joint balance is
`ΔH + ΔUelastic + ΔKsolid + ΔKe + ΔRef_e + ΔImage + Kwater + Uint_water + ΔUads − Wfield = 0`.
Native nodal heat and mechanical cell identities are separately checked; the
nodal-to-cell homogenization is not claimed as a common cellwise energy PDE.
All heat/layer/gas/path/inventory state commits and restarts with mechanical u/v
and source-stage/draw provenance. Lost acknowledgements require a matching
whole-candidate digest or checkpoint recovery.

Passed: **437 distinct Python tests**, **166 native v4 checks**, and **8,731
unchanged sealed-v2 regression checks against the working derivative**.
The legacy test invocation skips the 11 new tests, which separately pass.
The geometry test receipt transcribes its observed successful stdout explicitly.
Both signed four-step histories pass complete rollback and bitwise restart
continuation with an active source layer. Maximum joint energy residual:
**4.4072e-21 J**; maximum native nodal residual: **1.3538e-23 J**. Maximum
component proof ratio is below0.057. Full receipts are in `qualified-final`.

## Actual material-feedback pilot

The audit preserves 244 accepted states with independent upper/lower native
checkpoints and verifies canonical node ordering using MFEM DOF certificates.
The pilot's accepted thermal run eventually reaches2975.287K; that is outside
this mechanics qualification. The conservative first smooth-coefficient prefix
contains60 states, spanning0.749788055–12.428157577ns. All59 intervals pass the
32x4x2-per-foil conditional mechanics calculation. The final native range is
300.569674–370.569113K; the next actual stage crosses373.15K and is excluded.

| Support assumption | Minimum midplane gap change | Maximum sampled stress | Maximum proof ratio |
|---|---:|---:|---:|
| Free | −22.541192nm | 32.811453MPa | 0.210156 |
| Both end planes clamped | −23.284622nm | 27.107708MPa | 0.173624 |

Each foil retains nominal reference mass0.0002kg. Sampled Jacobian, volume/area
closure and mass guards pass. The comparison uses a separately hashed8000kg/m³
at300K density scenario, explicit Poisson ratio0.29, and prescribed Q1 corner
average temperatures. It does not infer missing supports, reconstruct missing
prehistory, project native caloric energy, debit heat, apply magnetic force or
feed geometry into EM. It is not spatially converged. The owner's completion
binary/material hashes match the observed snapshot; a formal source-to-build/run
seal is not established. See `actual-pilot-geometry-summary.json` and
`actual-pilot-geometry-mesh32-z2/conditional-geometry-history.png`.

## Limits and preservation

Mechanics remains298.2–773.15K, within one open coefficient segment per step,
and subject to nominal proof, strain, positive capacity and Jacobian guards.
The source surface is a fixed stationary reference approximation: emitted
support impulse is diagnostic, with no moving-wall recoil/work model.
Prescribed fields and caller event weights are component inputs. No physical
compressible EOS, plasticity, full PIC/gas transport, deformed field solve,
measured foil fixture or production readiness is claimed.

All numerical work used frozen CPU7/nice15, one thread, the shared lock and short
leases; geometry chunks contained at most8stages and13s inside a15s timeout.
Dependencies were not rebuilt; no GPU or production work was performed.

`source-manifest.json` seals106 inputs and four binaries. `integration.patch`
applies to the immutable v3 base. The compact Git evidence is
`artifacts/ss304-native-thermoelastic-v4-20261029`. Backup:
`/home/treefor/sim-workspace/source-backups/reconstructed-source-20261029/ss304-thermal-expansion-20261029/foil-support-v4`. It contains a verified selected-source/evidence archive, all four
binaries,1,195 verified raw evidence files and an incremental Git bundle.
The bundle requires v3 base`d491d934cc6d163a3a2536c2e0b45d9e3556314c`, available
in the unchanged sibling `foil-support-v3/source-and-results.bundle`.

Read `preservation.json`, `push-receipt.json` and `handoff.json` for machine
receipts. The initial preservation attempt found three owner logs excluded by
the repository's global `*.log` rule. They were explicitly committed, and the
final archive verification includes them. No numerical source changed.
