# Width-dependent desorption and VUV emission with MFEM Joule

MFEM `miniapps/electromagnetics/joule.cpp` is the selected electrode-heating
backend, replacing FLASH. The scientific target is unchanged: determine whether
finite-width heating and water release lead to observable spectra different
from uniform desorption, then estimate neutral inventory and temperature where
identifiable. Use the [backend requirements](mfem_joule_heating.md) for inspected
capabilities, limitations and the qualification ladder.

The chain is MFEM solid-electrode heating → conservative parent-water release →
owner gap transport and state-resolved Raikou → optical collection/gate/response →
recovered dissertation observations. Boltzmann, CR and kinetic-gap development
remain with their owners. The retired reduced MHD prototype is not restored.

## Evidence and experimental comparison

Shot 15750 is the pilot, with recovered H/H2 gap profiles and partial wavelength
lineouts. Shot 15765 provides the first additional wavelength comparison. Shot
17133 is a separate published heating benchmark; its drive cannot silently
replace a spectroscopy shot's drive. Original detector files are unavailable.
Use the [pilot manifest](../data/dissertation/pilot_comparison_manifest.json)
and [recovery/scenario contract](dissertation_pilot_comparison.md).

Derived transition-density displays retain literal units. Unresolved dimensional,
reduction and wavelength-response choices are separate named reconstruction
scenarios. The 41-point `10**linspace(-2,2,41)` multiplier applies only to
admissible nominal radiance and is an exploratory sensitivity range. It is not
noise, a confidence interval or a cure for dimensional ambiguity. Use maximum
asymmetric readout-envelope departure without artificial weight from dense
digitization. Preserve missing segments, clocks, spatial support and normalization.

Compare these three questions separately:

1. **Source localization:** compute the finite-width joint source and redistribute
   the same complete time/temperature/state/number/energy mixture uniformly over
   each emitting face. Keep excitation, transport treatment, observation and
   calibration common. Do not average temperature before nonlinear release.
2. **Heating model:** compare finite-width and independently specified uniform
   heating under matched drive. This may change total release and its timing.
3. **Excitation feedback:** when qualified owner models permit it, compare
   self-consistent excitation feedback with common prescribed-electron cases.

The measurement integrates along foil width. Under fixed linear full-width
collection, equal total upper-state histories give identical spectra. A hot-edge
image alone does not establish spectral discrimination. Report competing and
non-identifiable solutions, and use held-out shots or regions when possible.

## Electrode calculation and source handoff

Start from the unmodified MFEM 4.10 stock rod as a native execution/readback
demonstration. Its zero-temperature initialization, fixed coefficients,
sinusoidal voltage and artificial exterior are not a foil experiment. A physical
driver must provide the finite-width geometry, both conductors, current path,
edge radius, contacts, sourced SI SS304 properties and shot-specific registered
drive. Check direct Joule power, total conductor current and boundary/magnetic/
thermal energy accounting. Resolve skin depth, mesh, timestep and element order.

The initial physical scope is a stationary solid electrode. Stock Joule has
one-way electromagnetic-to-thermal coupling, no material-property feedback,
fluid motion or Hall induction. Qualify the valid temperature/time window before
using it for desorption. Extensions for property feedback, motion, phase response
or Hall physics need their own evidence. A numerical exterior must not become
the electron population of the chemistry model.

The [MFEM surface contract](mfem_surface_contract.md) feeds existing
`PlanarSurfaceHistory` and conservative release code. Record SI units, absolute
versus temperature-rise interpretation, clock, exact face geometry, true area,
one-sided metal DG boundary trace or positive-depth sample, and field/source
identities. Do not silently treat a bulk temperature as a boundary value or
average current-direction variations before nonlinear release.

Use [finite-inventory release](surface_release_pipeline.md),
[coverage-dependent release](coverage_dependent_release.md),
[planar patch handoff](mfem_surface_contract.md) and
[joint boundary quadrature](boundary_quadrature.md). Preserve each parent pool,
source-time window, molecular state, translational/internal/removal energy and
joint covariance. Wall debits must be applied once if two-way thermal coupling
is later implemented; the present handoff does not evolve MFEM temperature.

## Gap kinetics, radiation and inference

Retain three velocity components, source angular distributions, fragmentation
recoil and physical reaction energy partition. Check ballistic flight times,
wall events, nuclei/charge/energy budgets and magnetized charged transport.
Resolve rovibrational states for every retained molecular neutral, ion and anion.
Atomic H 2s and 2p remain separate. Missing rates are explicit blockers, not zeros
or invented closures. Use the [Raikou contract](raikou_receiver_contract.md) and
[external-model comparison plan](external_plasma_yacora_comparison.md).

Obtain an actual physical receiver receipt and complete time/cell/state histories.
Gas geometry determines ray-cell lengths. Apply normalized collection weights,
gate integration, line profiles and detector bins once. The original 5-mm foil
width convention is not a newly fitted slab depth. Justify optically thin
admission or use resolved transfer. Shared observation treatment is required for
Raikou and Yacora comparisons.

The [finite ensemble engine](neutral_ensemble_comparison.md) varies one parent-water
inventory and one explicitly defined neutral thermal coordinate. Complete forward
cases regenerate populations, transport and energy; no solved-population scaling
substitutes for new density/temperature physics. Translational, rotational and
vibrational interpretations stay distinct. Derive coverage from parent inventory
and a sourced monolayer convention, without fitting another density amplitude.

## Completion evidence

Advance through native smoke/readback, analytic electrical/thermal benchmarks,
converged physical foil/source cases, accepted physical Raikou ensembles and then
dissertation comparison. Preserve passing evidence and rerun affected checks.
Require independent mesh/time/spectral/state refinement below 1% in relevant
projected features and numerical uncertainty below 10% of any claimed contrast.
Successive differences alone do not establish that uncertainty bound.

The earlier [FLASH specification](flash_h2o_publication_plan.md) and setup records
are retained as design history. FLASH download is not required for this plan.
Primary implementation sources are the official
[MFEM description](https://mfem.org/electromagnetics/),
[v4.10 driver](https://github.com/mfem/mfem/blob/v4.10/miniapps/electromagnetics/joule.cpp)
and [solver](https://github.com/mfem/mfem/blob/v4.10/miniapps/electromagnetics/joule_solver.cpp).
