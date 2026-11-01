# Native thermoelastic v3 validation

Source commit: `b85475df8594cb224d1347f951558b783e9847f1`.

418 Python tests and 76 native MFEM assertions passed. Both manufactured
four-step conductive cases conserve internal plus kinetic energy within their
recorded tolerances and reproduce continuation bit for bit after restart.
Joule, authored, wall and species heat are exactly zero. Temperature spans
298.9–300.1 K initially; mechanics remains capped at773.15 K.

The native executable SHA is in source-manifest.json. It links the unchanged
frozen MFEM archive. Native binary backups are outside the repository.

The independent review recommendations were resolved: source/density/stage
bindings, local representability bounds, proof-screen and initial-mechanics
identity, staged restart admissibility and atomic rejection. The review is
retained as evidence of the review, not an unresolved finding list.

This is a fixed-reference-grid microscale numerical binding benchmark. It
is not a production foil, moving mesh, deformed EM, physical EOS, plasticity,
or fixture qualification. It is based on sealed heat API v1; live pending
v2 is not imported. See compatibility.json for owner handoff.
