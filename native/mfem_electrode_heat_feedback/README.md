# Evolving electrode heat and adsorbed inventory

Native MFEM geometry/load quadrature with bounded nonlinear enthalpy, monotone
orthogonal dual-cell conduction and mass-lumped storage. Actual weighted surface
collection/solid-volume transfers and actual emitted packet energies change T.
Returned accepted boundary T, normal, area and sigma(T) feed the next source/field
stage. This is a component, not an independent PIC or full breakdown driver.

See [API-v1.txt](API-v1.txt) for transaction, lifetime, two-limb inventory, event
identity, source and energy contracts. Parent preflights all participants, then
commits accepted states together. No callback/RNG/production ownership is taken.

## Conditional material

`ss304_heat_material.hpp` consumes the property owner's frozen native callbacks.
The tested selection is explicit `nominal_synthesis_1800`, fixed-reference-volume
Cp-as-Cv and caller reference density8000kg/m3 (source-model approximation, not
measured foil density). Transport pools304/304L and is source-synthesized above
1500K; mushy/phase assumptions are retained. No vapor/deformation claim. The
separate mechanical773.15K limit is not lifted. All range failures abort trials.
Source property/scenario data and original seals are retained under `material/`.

## Reproduce

Use the frozen launcher and MFEM4.10 installation in the receipt. Configure this
directory with `-DMFEM_DIR=.../install-mpi/lib/cmake/mfem`, build with one job, then
run `electrode-heat-feedback-tests OUTPUT` under CPU7/nice15, one-thread env and
`outputs/pic-profile/testing-core7.lock`. `--matched` runs one bounded native
256x32 EM/thermal reference interval. No GPU or long heating rerun is required.

Evidence includes1740 component assertions plus4 native matched-reference checks,
exact restart/rollback, ABA placement reconstruction, canonical energy-owner
deduplication, compensated inventory, proper-velocity/relativistic packet K,
thermal/inventory coupling-interval and mesh/time refinement, and direct material
phase/transport callbacks. Local temperature differences from the older consistent
mass formulation are reported; integrated native Joule input is identical.

Reference files are unchanged copies of the current native edge source used only
by the tests; their exact hashes are in the report input manifest. The controlled
reference current and current-controlled heatplate are component fixtures, not
the measuredC19 application history. The new application drive remains the
measuredC19 total-current trace supplied by the field owner; this adapter imposes
no voltage waveform. No physical H2O law or four-reservoir defaults are installed
by this generic API. Source-law and external-work binding are a separate version.

Stable read-only source snapshot/API and full logs are in the report prefix
`outputs/pic-electrode-heat-feedback-20261029`. No production launch occurred.
