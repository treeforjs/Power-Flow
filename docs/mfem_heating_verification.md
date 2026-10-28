# MFEM analytic heating-profile comparison

The existing `HeatingProfile` comparator now accepts MFEM depth profiles without
labeling them as FLASH. Numerical reference formulas, cell-average semantics and
error norms are unchanged. See the independent
[magnetic diffusion reference](flash_heating_verification.md).

An MFEM profile declares `source_kind=mfem`, `mfem_version` and a full
`source_revision` commit, with empty `flash_version`. Its CSV sidecar uses
`mhdlab-heating-profile-v2` and `heating_backend=mfem_joule`. An explicitly
`synthetic` profile may carry the same backend metadata for an authored
interchange fixture; that label remains synthetic. Neither case authenticates
an external execution or material model.

Earlier FLASH and generic synthetic profiles retain their exact v1 metadata
format. Both formats use the same fixed-SI CSV columns and SHA256 linkage.
Missing MFEM identity, mixed FLASH/MFEM provenance, an inconsistent backend or a
v1 downgrade of an MFEM profile reject. Native dimensionless Joule fields are
not converted to SI by this comparator.

```python
from mhdlab.heating_comparison import load_heating_profile_csv, compare_heating_profile

profile = load_heating_profile_csv("/path/to/qualified_mfem_depth_profile.csv")
report = compare_heating_profile(profile, matching_constant_property_reference)
```

The reference must match the numerical case's coefficients, magnetic drive,
initial state, computational depth and sampling. Integrate B, J and Joule power
separately over each cell; squaring a cell-averaged current is not the cell
average of J squared. A one-sided DG surface trace is a different quantity.

Tests exercise synthetic MFEM roundtrip, missing/mixed identity rejection,
downgrade rejection and exact comparison with analytic cell averages. They do
not constitute an MFEM calculation of the analytic diffusion problem. A qualified
physical solver benchmark and its independent current/energy/convergence evidence
remain required before the finite-width foil case.
