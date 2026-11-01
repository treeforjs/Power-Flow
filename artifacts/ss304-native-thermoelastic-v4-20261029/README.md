# Sealed heat/source v2 with signed mechanics

Source commit: `85ba2a1640a56dfd11d9534ec3e5de881af3b0ec`.
Base: `d491d934cc6d163a3a2536c2e0b45d9e3556314c`. All source changes are
new v4 files. V3, the baseline material and the heat owner's live source are
unchanged. `source-manifest.json` seals 106 source/dependency files and four
native executables; `files-sha256.json` seals this compact evidence package.

The v4 component combines signed mechanical exchange, native conduction,
actual sealed electron birth/continuation and actual water desorption. The
nineteen thermal ledgers retain electron reference/image/field-work ownership,
zero gas reference and the negative adsorbed reservoir. Joule/authored heat is
zero. Complete heat/layer/gas/path/mechanical rollback and restart are tested.
The two signed four-step histories have a maximum joint residual of
4.4072e-21 J and a maximum native nodal residual of 1.3538e-23 J. Full receipts
and active-layer common checkpoints are in `qualified-final`.

Qualification passed 166 native v4 checks, 8,731 sealed-v2 checks compiled
against the working derivative, 11 new joint Python tests, 418 existing Python
tests and 8 geometry projection/admission tests. The 11 skips in the existing
suite are the separately executed v4 tests. No tests are silently omitted from
the 437 distinct Python passes. All numerical work used locked CPU7, nice15,
one thread and the frozen runtime; no GPU or production run was used.

The actual material-feedback pilot comparison uses independently preserved
upper/lower temperature checkpoints and native coordinate certificates. All
59 intervals in the first smooth-coefficient prefix passed through 12.428 ns,
ending at 300.570–370.569 K. The next stage crosses 373.15 K and is excluded.
Minimum gap reductions are 22.541 nm for free foils and 23.285 nm for clamped
end planes. Maximum sampled stresses are 32.811 and 27.108 MPa. See
`actual-pilot-geometry-summary.json`, its text companion, and the figure.

The pilot comparison prescribes reference-cell average temperature and assumes
nominal density 8000 kg/m³ at 300 K, Poisson ratio 0.29 and explicit support
brackets. It does not infer the fixture, project a native caloric ledger,
apply magnetic forces or feed geometry back into EM. It is not spatially
converged. The owner's completion binary/material hashes match the observed
snapshot; a formal source-to-build/run seal is not established.

The source surface remains a fixed stationary reference approximation; emitted
support impulse is diagnostic. Moving-wall momentum/work coupling, transport
after PIC handoff, physical compressible EOS, plasticity and production field
handoff remain outside this qualification. Mechanics is restricted to
298.2–773.15 K, a single open coefficient segment per step, positive sampled
Jacobians, strain and nominal proof-screen guards.

Full raw native checkpoints, all 59 geometry states and runtime binaries are
preserved outside Git under
`/home/treefor/Documents/Codex/2026-10-28/we/outputs/ss304-thermal-expansion-20261029/foil-support-v4`.
The final `HANDOFF.md`, `handoff.json`, preservation receipt and push receipt
there bind this evidence to the results commit and backup. The backup contains
a selected-source archive and an incremental Git bundle based on the already
preserved v3 commit, plus complete raw evidence and all four binaries.
