# MITL thermal boundary and material-input handoff

Date: 18 October 2026

This packet extends the upstream source interface through law-specific emission
times and paired three-dimensional launches. It also recovers a bounded primary
304 transport table and supplies its explicit input adapter. Public source is
pinned to Power-Flow commit `ca986cf4114ff06adf8415feb2fbc0d3f8185daf` on
`research/flash-publication-foundation`. Five local commits separate input
validation, emission timing, planar sources, table recovery and material inputs.
They are included in `boundary-extension.patch`.

The previous frozen source/Hall handoff is
`/workspace/handoffs/mitl-upstream-extension-20261018T1747.tar.gz`; the original
baseline/native packet remains at
`/workspace/handoffs/mitl-heating-source-handoff.tar.gz`. This packet preserves
both. Its source extraction contains 55 byte-identical public files plus named
packaging wrappers; it does not import a legacy runner.

## Delivered interfaces

- `FrozenTemkinEmission` conditions on the actual coverage-dependent emitted
  population. CDFs and quantiles use the separated nonlinear law, retaining
  tiny fractions, early events in stiff intervals, and near-unity probabilities.
  Positive times with unsupported float64 precision fail explicitly;
  `log_quantile` preserves their logarithms. This completes the temporal adapter
  that was absent from the preceding packet.
- `sample_temkin_emission_packets` conserves each positive reservoir/width
  component with unequal weights and at least one packet. It uses the qualified
  conditional time distribution and outward flux-Maxwellian thermal velocities
  at each component's temperature and supplied exact-state mass.
- `TemkinPlanarReleaseHistory` applies local nonlinear release on current/width
  patches before area integration. Its schema retains geometry, proxy/depth
  evidence, nested raw extraction records, independent-pool assumptions and
  numerical qualifications. HDF storage and shapes are checked before reads.
- `sample_temkin_planar_emission_packets` allocates across those actual patches,
  uses each row's length once, and produces paired 3D coordinates. Resolved and
  uniform-width hypotheses retain identical y, z, state, mass, time, momentum,
  velocity and weight; only width positions differ. Parent metadata is retained
  once, with qualified row references rather than duplicated raw manifests.
- The surface input helper now rejects masked arrays. History validation checks
  the compensated conservation residual and replayed diagnostic regime, and
  requires a zero survivor for the declared underflow regime. Equivalent root
  iteration counts and qualified residual values remain portable.
- A guarded recovery supplies all ten published smoothed Flynn 304 transport
  rows, preserving printed precision and exact Celsius/SI conversions. PDF,
  text, rasters and reviewed renderer crops are hash-identified. Documented
  OCR/unit repairs and their evidence are explicit; no missing points are added.
- `BoundedTransportTable` and the named Flynn adapter preserve grade, sample,
  solid phase, status and source identity. Resistivity and thermal conductivity
  interpolate within the declared interval; electrical conductivity is derived
  afterward as the reciprocal of resistivity. Exact-ratio arithmetic avoids
  losing finite contributions through underflowed interpolation weights.

The guides are `source/docs/temkin_boundary_interface.md` and
`source/docs/material_transport_inputs.md`. The source subset includes runnable
width and planar workflows, packet exporters and the recovery tool. Its README
describes independent installation and testing.

## Physics and integration responsibilities

Every Temkin pool uses an explicit fixed reference coverage and its own local
remaining population. No shared-site competition or state-population allocation
is inferred. Splitting a physical adsorbate layer among molecular states changes
this law unless the physical occupancy model preserves the intended coverage.
Exact species/state identities, masses and internal energies remain inputs from
the molecular owner. No collision, radiative or kinetic-gap code is changed.

The kinetic barrier sets the rate; it is not an additional heat cost. The
separately supplied constant desorption cost, assigned gas-state internal energy
and expected translational energy `2*kB*T` are charged once. Sampled velocities
retain their actual weighted kinetic energy, the expected mean and their
sampling difference. A coupled host must balance its applied launch energy
against the wall or use a verified moment-conserving representation. It must
also define its matching mass/momentum coupling. These kernels do not apply wall
cooling, momentum feedback or motion.

Both spatial hypotheses retain the full source-temperature and emission-time
mixtures. Do not rethermalize the uniform control at a mean temperature or
resample its states/times independently. Unequal weights need native support or
a separately verified adapter. The old first-order source loader/sampler still
rejects Temkin records; the new adapter uses actual Temkin history types.

Source histories freeze linearly interpolated midpoint temperatures. Refining
this integration alone cannot recover the missing FLASH snapshot cadence or
spatial structure. The core's dimensionless coverage coupling remains qualified
through 512. Planar controls reject positive SI moment densities below float64's
normal range, where subnormal rounding can change conservation ratios; exact
zero remains supported. These are declared numerical domains, not physical rate
floors or calibrated surface parameters.

The optional source-run summary must hash-match its exact HDF5 source before
the manufactured/sourced/unvalidated parameter label enters packet metadata.
Without that associated record the label remains unreported. Heating-kind and
mass labels are retained independently. NPZ files contain numeric/Unicode arrays
and can be read with `allow_pickle=False`.

The Flynn table applies to one supplied solid plate sample at **253.15–423.15 K**,
with unreported chemistry and thermal history. Its values are published smoothed
measurement results, not the dissertation foil's raw data or universal 304
defaults. Source error estimates retain their original interpretation; no
independent/Gaussian uncertainty model is created. The adapter supplies no
density, heat capacity, enthalpy, phase transition or high-temperature extension.
An illustrative 700 K query is rejected. Source PDFs are not redistributed.

## Verification and delivered evidence

The new boundary/recovery scopes passed **235 tests** under warnings-as-errors;
the bounded material contract separately passed **128**. Affected baseline and
history compatibility checks passed 106 tests. Exact commands, exits, logs,
artifact metrics and extraction results are in `verification.json` and
`evidence/`. Unchanged whole-project suites were not repeated.

The packet includes manufactured width and planar source outputs, safe launch
archives, a conditional-time illustration and bounded material-input output.
The 4096-packet planar example retains 144 source components. Its expected and
sampled thermal energies remain separate; their relative sampling difference
is approximately −0.54%. This is numerical sampling behavior, not a measured
energy discrepancy or water fit. The timing illustration likewise compares
manufactured laws; it is not an experimental desorption result.

Source, implementation, output and delivery hashes identify the frozen bytes.
The isolated subset is tested with the external source PDF deliberately absent,
so the optional recovery integration can skip without preventing data use.
That isolated check passed 362 tests, with only that one optional test skipped.

No installed licensed FLASH source/configured workstation was available, so no
actual FLASH run is claimed. No Quark source, build, runner or production job
was changed or used. No remote push, merge, publication or external message
occurred. Native port/registration remains an owner handoff.

Next upstream work is a resolved x/z/time convergence study of the complete
Temkin source/launch chain, followed by additional grade-specific high-temperature
property evidence and the remaining heating verification ladder.
