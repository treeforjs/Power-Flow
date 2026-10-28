# MFEM Joule backend switch

MFEM `joule.cpp` replaces FLASH as the selected electrode-heating backend.
Use the [current research plan](../../docs/mfem_h2o_publication_plan.md),
[backend requirements](../../docs/mfem_joule_heating.md) and
[surface contract](../../docs/mfem_surface_contract.md).
The earlier FLASH and neutral-comparison records retain their original scope;
they are not rewritten to claim evidence for MFEM.

The local MPI MFEM 4.10 installation was reused. An isolated consumer compiled
the unchanged upstream Joule driver, solver and helpers, then ran the stock 3D
rod on two MPI ranks. No installed MFEM library or shared source was changed.
The saved-field reader loaded initial/final VisIt fields `T,w,E,B,Phi,F`.
They were finite, with positive final integrated temperature and Joule heating.
Native final time was 0.2; final temperature minimum was -0.00128638 and maximum
0.178637. The negative undershoot is retained. No Kelvin/SI interpretation,
positivity, energy convergence or physical foil qualification follows from this
demonstration.

The SI table adapter preserves absolute-versus-increment Kelvin semantics,
one-sided material-trace/sample-depth association, x/z patches, true face
geometry and MFEM provenance. It feeds existing planar and Temkin conservative
release without another depth/current-length factor. A separate manufactured
table exercises this handoff; it is not the native rod relabeled as a foil.
The analytic heating comparator admits explicit MFEM provenance separately from
earlier FLASH profiles.

[completion.json](completion.json) records exact commands, results and checksums.
Compact accepted reports are copied here; full native fields, logs and HDF5
outputs remain under ignored `runs/`. Current native evidence originated in
`runs/mfem_joule_native_demo/`, including `installed_identity.json`, which records
the independent library/configuration identity. A matching source/build
attestation is not inferred from file hashes alone.

Remaining physical work is the experiment-specific finite-width SS304 driver,
shot drive/current closure, material and thermal-boundary qualification,
temperature feedback where needed, analytic/conservation/convergence checks,
and qualified one-sided planar face evaluation. Stock Joule has no fluid motion
or Hall induction. Raikou state/rate/source acceptance and experimental reduction/
collection/clock requirements remain unchanged. FLASH acquisition is no longer
a blocking dependency.
