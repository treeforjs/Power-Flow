# MITL heating and boundary-source handoff

Date: 18 October 2026

Status: independent upstream reference implementation. No Quark source,
integration branch, build, runner or production job is changed or used.
Boltzmann transport and collisional-radiative development have separate
ownership. This packet supplies heating inputs, recovered evidence and
conservative boundary-source interfaces to those workstreams.

## Delivered contracts

- Explicitly configured FLASH leaf-output extraction, with physical axis/unit
  mapping, material purity, strict pointwise coverage and source-file hashes.
- Persistent distinction between boundary temperatures, positive-depth material
  samples and explicitly justified boundary proxies. Current-direction averaging
  cannot silently become a resolved nonlinear desorption source.
- Finite-inventory first-order binding reservoirs, exact frozen-temperature
  depletion, compensated inventory residuals and separate rate/energy parameters.
- Width- and current-resolved planar sources. Release precedes area averaging;
  physical patch areas determine extensive number and energy totals.
- Uniform-width controls preserving each time/reservoir/current-patch number
  and energy budget. Full-face controls are separate explicit hypotheses.
- Paired thermal launch packets retaining source states, times, velocities and
  component inventories while randomizing only width position.
- Recovered published temperature, field/current-density and electrical drive
  curves with hashes, calibrated coordinates and reconstruction limits.
- Analytic constant-property magnetic diffusion, finite-volume field/current/
  Joule comparisons, and independent magnetic/Joule/Poynting energy accounting.
- Standalone C++20 surface-release kernels and physical verification groups,
  delivered as native candidates without host registration.

The mathematical model is generic. Manufactured examples do not supply physical
water desorption parameters, a complete molecular state registry or spectra.
No real FLASH calculation is claimed by the packet.

The implementation is pinned to Power-Flow commit
`3b5e7d12a5d993122197330a33053136dbd390d0`. The isolated Python package passed
560 tests. The exact copied modules, tests and data are identified by
`source-files.sha256`; extraction-specific wrappers are listed in the
verification record. Native candidate bytes are pinned separately.

## Integration responsibilities

1. Treat configured FLASH extraction as a supported-layout reference, not a
   generic moving-material surface finder. Verify actual output fields, lattice
   temperature, materials, axes, units and sample-plane geometry locally.
2. Keep physical species/state identities and their gas-reference energies
   consistent with the molecular owner. A surface pool split among levels must
   conserve its original inventory.
3. Debit binding, assigned internal energy and launch translation exactly once.
   Activation barriers determine rates and are not additional heat costs.
   Prescribed-temperature diagnostics do not apply cooling to FLASH automatically.
4. Preserve local surface temperature/inventory correlations. Average released
   extensive quantities after local nonlinear events, rather than averaging
   temperatures before them. Retain patch areas and sampling-depth provenance.
5. Boundary packets use unequal component-specific weights to conserve each
   reservoir/patch source. A downstream collision solver must support them or
   supply a verified weight adapter. No adapter to the old equal-weight gap
   reference is included.
6. Thermal packet energy fluctuates about its expected flux-Maxwellian moment.
   Keep sampled-versus-expected diagnostics visible. The paired source controls
   share identical samples; never independently rethermalize their velocities.
   A finite-packet host must balance the actual weighted launch energy with
   its applied wall counterpart, or use a verified moment-conserving source.
   Keep the analytic mean and sampling difference as separate diagnostics;
   the difference is not an additional physical heat source.
7. Accept source inventory, applied wall-energy exchange, packet creation and
   cumulative transfer diagnostics transactionally in the native host.
8. Native registration and runtime verification remain with the integration
   owner. No production runner work is requested by this handoff.

## Numerical and physical limits

Each time interval freezes an interpolated midpoint temperature. Convergence
requires refining this interval, the underlying FLASH snapshot cadence and the
surface partition. Float64 inventories retain compensated residuals when tiny
decrements lie below their stored spacing. The current release law does not
include coverage-dependent barriers, bulk diffusion, readsorption or evolving
surface chemistry.

The published thermal profiles used shot 17133; the recovered electrical trace
is shot 15210. They are separate benchmarks, not a same-shot reconstruction.
Field/current-density literal-legend conversions remain unverified and cannot
establish a quantitative Joule-heating benchmark. Recovered raster/vector
coordinates do not restore original detector or oscilloscope data.

The accompanying source manifests and verification record pin the delivered
implementation and evidence. Local commands and file layouts are described in
the source run guides. They invoke independent reference programs only.

## Package layout

- `source/`: isolated upstream modules, tests, configurations, run guides and
  recovered figure/drive data. Its initializer imports no legacy runner.
- `source/examples/paired_planar_packets.py`: executable three-dimensional paired
  boundary example using an existing planar release output.
- `native/`: independent C++20 surface-release candidate and standalone harness.
- `evidence/`: verification logs and selected manufactured/reconstruction outputs.
- `verification.json`: pinned source commit, versions, commands and observed results.
- Source and delivery hash manifests identify the exact packaged bytes.

The full project suite passed 721 tests with two optional tests skipped. The
native surface candidate passed all 11 runtime verification groups and its
standalone CTest target. The refined manufactured width study reports no
non-asymptotic warnings; its temporal orders approach two. Numerical convergence
of that fixture does not validate a physical desorption parameter set.
