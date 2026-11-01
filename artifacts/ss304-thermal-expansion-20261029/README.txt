SS304 thermal expansion and geometry benchmarks — 2026-10-29

Implemented integrated mean-alpha eigenstrain with explicit stress-free T,
bounded E(T), static/free/restrained and transient rod mechanics, affine foil
areas/normals/gap/rho(J), conservative remapping and local geometry transactions.
Source and implementation: docs/ss304-thermal-expansion.md; material record:
data/materials/ss304-v1.json. Source commit and file hashes: source-manifest.json.

Validation:229 passed (101 new mechanics/source/mapping/port checks plus128
existing material-transport checks), frozen CPU7/nice15/one thread/shared lock.
White figures: uniform-displacement-gap.png, transient-displacement-gap.png,
source-admissibility.png. Full input/output hashes and refinement values saved.

Exact304/S30400/1.4301; selected mechanics-data range293.15..773.15K, not melting.
New matched256x32 T/Az/rate run used exclusively. Last common temperature-safe
frame17.413910165ns; first invalid surface17.622896529ns. Full history max
surface13300.599K. Lower bulk is explicitly mirrored, surface-checked, unverified.

Mechanical supports are unspecified. Source-driven curves are independent
free thickness-bar ELASTIC TRIALS, not connected-foil predictions. Edge trial
crosses annealed-bar proof proxy about10.08ns; plasticity is absent. All physical
moving-field handoff remains blocked. No production launch or paused-run signal.

The source package enforces interval-end heat stages, revision/hash/mass identity
and rollback/replay protection in-process. No checkpoint persistence or parent
global atomic coordinator is claimed. Heat debit and mechanical work are null
until a consistent caloric/work law is bound; energy diagnostics are not closure.
Separate high-T thermal-properties work owns enthalpy/transport/phase data.
No PDFs are included here; producer source bytes remain in the outside-project
report/material-sources directory, with frozen hashes/provenance in the record.
