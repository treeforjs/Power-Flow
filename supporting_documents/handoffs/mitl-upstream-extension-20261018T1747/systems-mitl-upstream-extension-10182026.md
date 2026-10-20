# MITL upstream extension handoff

Date: 18 October 2026

This addendum supplies coverage-dependent release, a Hall-induction verification
reference and bounded stainless-property evidence. It is pinned to Power-Flow
commit `dc33c35a944457eaf685b28268a3fe3192cca3e9` on
`research/flash-publication-foundation`. The three local commits are separated
by feature and are included in `upstream-extension.patch`. No remote push,
merge or publication has occurred.

The preceding baseline/native handoff remains at
`/workspace/handoffs/mitl-heating-source-handoff.tar.gz`. This addendum preserves
that frozen package. The `source/` subset here is independently runnable and
contains 29 byte-identical public-project files plus named extraction wrappers.
Its initializer imports no legacy production runner.

## Delivered code and evidence

1. **Separate coverage-dependent first-order release.** Each pool uses the
   remaining population in the prefactor and a supplied kinetic barrier
   `Ed*(1-alpha*N/Nref)`. The reference density and maximum supported inventory
   are explicit and independent of initial coverage. All kinetic and energy
   parameters and their distinct provenance are required. There are no water
   defaults or SS304 calibration claims.
2. **Monotone frozen-temperature integration.** The kernel inverts a separated
   log-depletion integral and retains tiny-depletion, survivor and inventory
   diagnostics. Nonlinear coupling `Ed/(kBT)*alpha*N/Nref` is qualified through
   512; larger coupling raises a named unsupported-range error. The declared
   coverage domain must retain a nonnegative barrier throughout. No clipping or
   silent initial-rate substitution occurs.
3. **Distinct source history and interchange.** Midpoint temperature histories
   retain number, separate binding/internal/translational energy, numerical
   inversion diagnostics and full FLASH association metadata. The separate
   HDF5 schema declares law, time model and numerical qualifications, rejects
   external data dependencies, and checks local kinetics and energy ledgers.
   Equivalent qualified numerical roots are portable across library iterations.
4. **Two different controls.** Matched uniform-width moments redistribute the
   resolved source after release while preserving interval/state number and
   energy. The optional constant-barrier control separately sets alpha to zero
   while retaining the initial inventory and wall-energy costs. It is a kinetic
   sensitivity, not the spatial null hypothesis.
5. **Periodic Hall-induction reference.** An exact helical mode checks Faraday
   and Ampère signs, workless Hall rotation, scalar resistive damping, independent
   cell-mean Joule power and perturbation magnetic/Joule energy exchange. Point
   and finite-volume fields, manufactured export configuration, CSV sidecars and
   hashes are included. This reference prescribes stationary ions and qualifies
   the induction subsystem only; it is not a full Hall-MHD or gap solution.
6. **Material-source evidence.** Eight source records distinguish measured 304,
   separately identified 304L estimates/abstracts, extrapolated compilations,
   unavailable full text and a rejected wrong-alloy SRM substitution. Ten cached
   artifact hashes and PDF identities were checked. A verified 1700 K liquid
   enthalpy conflict in ANL-75-55 is preserved. No constitutive table is adopted.

The run guides are `source/docs/coverage_dependent_release.md`,
`source/docs/hall_induction_verification.md` and
`source/docs/ss304_property_sources.md`. The acquired PDFs remain outside this
packet; their official or explicitly identified mirror URLs and hashes are in
the source manifest. No complete source texts or figures are redistributed.

## Interfaces and physical limits

`TemkinReservoir`, `TemkinReleaseTrial` and `TemkinReleaseHistory` are distinct
types. The ordinary first-order source loader and packet sampler refuse them.
A truncated exponential emission-time CDF is generally wrong for this new law.
**No law-specific packet-time adapter is supplied yet.** Consume the integrated
moments only under the identified law; do not fabricate a baseline history to
bypass the gate. This temporal adapter and a current-resolved planar extension
are the next upstream integration tasks.

Each pool's coverage controls its own barrier. Splitting one physical adsorbate
layer among gas-state records does not automatically preserve that law. Shared
coverage/site constraints and molecular-state allocation need a documented
physical model from their owner, without duplicating inventory. Boltzmann,
collisional-radiative and kinetic-gap implementations are outside this work.

Kinetic barriers do not define surface internal energy. This extension retains
a separately supplied constant adsorbate-to-gas-reference cost, assigned state
energy and mean outgoing translation `2*kB*T`, counted once. The reported wall
debit is diagnostic under prescribed temperatures; FLASH cooling and coupled
thermal feedback are not applied. A future finite-packet host must balance its
actual weighted launch energy and retain sampling defects separately.

Positive-depth FLASH samples still require an explicit boundary-proxy reference.
Whole-current averages require explicit uniformity evidence before nonlinear
release. This new history is width-resolved only. It does not replace the
baseline's already delivered planar source or make averaged temperatures local.
Moving boundaries, readsorption, surface reactions and bulk diffusion remain
outside these kernels.

The Hall reference has constant scalar resistivity/permeability/electron density,
no electron pressure or inertia and prescribed stationary ions. Free ions would
respond to its Lorentz force; its mode cannot qualify their coupled dynamics.
Background magnetic energy stays separate from perturbation energy. Compare
independently averaged energy and Joule power, not products of mean fields.

## Verification and host handoff

The public-project targeted runs passed 141 Temkin tests and 94 Hall tests with
warnings treated as errors. They include independent ODE/integral, field-law,
cell-quadrature and energy checks, strict domain/schema tests and complete export
workflows. Exact commands, exits, logs, hashes and the isolated extraction result
are recorded in `verification.json` and `evidence/`. Existing passing whole-project
suites were not rerun for this independent addition.

The isolated source extraction also passed all **235 tests** with warnings
treated as errors. Reference-output hashes match both the project and delivered
modules. Delivery manifests identify all frozen files and extraction wrappers.

The Temkin and Hall example outputs are explicitly manufactured numerical
fixtures. They establish implementation behavior, not measured water release,
foil properties or spectra. The previously recovered thermal profiles used shot
17133; the electrical trace is shot 15210. Their separate evidence limits remain
unchanged, and neither this addendum nor the previous packet restores raw data.

No installed licensed FLASH source or configured target workstation was available
for execution. No actual FLASH run is claimed. No Quark source, build, runner or
production job was changed or used. This packet is an owner handoff for review;
native registration and host verification require a separate integration pass.
