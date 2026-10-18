# Width-dependent water desorption and VUV emission in planar MITLs

Working research specification, 18 October 2026. The agreed goal is to test whether
finite-width electrode heating changes water release and fragment transport enough
to produce distinguishable recorded spectra relative to a uniform-desorption model.
The new spectroscopy model must resolve the H2 rovibrational structure that the
previous CR analysis did not resolve.

## What exists and what remains

This branch provides a dissertation recovery inventory, an initial recovery of
published profiles, evaluated water collision data where available, and a validated
SI surface-history interface for external FLASH calculations. These are preparation
and data products. There has been no FLASH run, no full water-chemistry solve, no
rovibrational population solve, and no demonstration that the spectra distinguish
the hypotheses. The old `MHD_heating.py` path still executes the legacy prototype;
the new surface interface is not connected to that runner.

The user's FLASH source and workstation are available; workstation setup is planned
for Wednesday, 21 October. Hardware sizing and a version-specific FLASH setup remain
to be determined from a small benchmark. No licensed FLASH implementation is included
in this public-code repository.

## Publication question and comparison

The dissertation already reports preferential edge heating with ALEGRA in section
4.5, Figs. 4.27–4.30, printed pp. 93–98. A change of MHD code or another hot-edge
image alone does not establish the proposed advance. The useful advance is the
quantitative chain from electrode temperature to water flux, fragment distributions,
resolved populations, and the measured diagnostic signal.

Use three distinct comparisons:

1. **Spatial-source test:** compute the finite-width source Gamma1(x,t) from FLASH
   surface temperature and a common surface model. Construct the uniform null
   Gamma0(t) = integral[Gamma1 dA] / A, separately for each emitting face. This
   preserves the release rate at every time while changing its spatial distribution.
   Keep electron closure, chemistry, instrument model, and nuisance-parameter priors
   identical. This is the primary test of spatial localization.
2. **Heating-model test:** compare finite-width heating with an independently defined
   planar/uniform-heating calculation under the same drive. This can change total
   release and timing. It answers a different question; averaging temperature before
   applying an exponential desorption law does not preserve release.
3. **Excitation feedback test:** if the gap model is self-consistent, allow the changed
   neutral source to alter electron kinetics and fields in both hypotheses. Compare
   this result with a common prescribed-electron calculation to separate geometry
   and plasma feedback.

Report edge and center release per unit area, onset times, released inventory, and
the fraction of emission contributed by each width interval. Define edge intervals
before inspecting the result; check sensitivity to their width. A hotter edge does
not guarantee greater *integrated* release after a finite inventory is exhausted.
Do not insert an edge temperature or coverage enhancement to force the result.

The VUV images resolve distance across the A–K gap and integrate emission along the
foil width. They do not directly resolve a width-dependent source. Fit the projected
signal, test recoverability with synthetic observations at the actual resolution,
and report non-identifiability if a uniform source can reproduce the observations.
Use held-out shots or spectral/spatial regions and joint parameter constraints;
avoid separately adjusting arbitrary scale factors for every lineout.

## Electrode calculation in FLASH

Use a Cartesian cross-section with x across the foil width, y normal to the surface
and through the gap, and z along the current. Include both conductors and their
finite-width magnetic geometry. Resolve foil thickness, magnetic diffusion length,
and physical edge radius. Contacts, dielectric support, finite length, and current
feed effects need either representation or an explicit approximation with sensitivity
tests. The thesis ALEGRA example used shot 17133 and omitted the dielectric support;
it is not the measured drive for every VUV shot.

Required physics includes resistive magnetic diffusion and Joule heating, appropriate
cold-metal electrical/thermal properties, heat capacity, and material response as
melting becomes relevant. For a solid surface, desorption responds to the surface
lattice temperature, not automatically the electron temperature in a three-temperature
model. Validate room-temperature initialization before using FLASH's tabular EOS or
plasma transport closures for SS304. Do not use Spitzer resistivity for a cold foil.

At each time verify integral[Jz dA] = I(t), with Jz = (dBy/dx - dBx/dy)/mu0 in SI.
Distinguish electrical resistivity eta [ohm m] from magnetic diffusivity eta/mu0
[m2/s]. FLASH field normalization and CGS-to-SI conversion must be checked explicitly
against the installed version. A uniform magnetic boundary that removes finite-width
fringing is not a suitable boundary for the edge-heating test.

Start with magnetic-diffusion and Joule-heating benchmarks. The documented implicit
resistivity solver generally uses a uniform grid; verify version-specific AMR support
before committing to it. Test vacuum-floor density/resistivity, interface treatment,
domain extent, edge rounding, timestep, and mesh refinement. Inspect deposited energy
and temperature timing as well as images. An artificial vacuum must not supply the
electron population for the chemistry.

The [initial coupling contract](flash_surface_contract.md) is `mhdlab/flash_surface.py`: a stationary planar face,
one face per file, SI temperature histories with provenance and explicit width-bin
edges. It is an interchange format, not a raw FLASH AMR reader. A local exporter must
locate the physical material face, convert units, and project onto these bins without
mixing in vacuum or bulk-averaged temperatures. Moving surfaces, normals, electromagnetic
fields, and deposited particle energy require additional explicit contracts; this
initial temperature format is not sufficient for fully coupled discharge modeling.

## Surface release and kinetic transport

Begin with a finite water inventory N(x,t) [molecules/m2] on each face. An illustrative
first-order thermal law is Gamma = nu N exp[-Ea(theta)/(kB Ts)], where theta is
dimensionless coverage. Coverage-dependent binding, multiple site energies, and
water/hydroxyl surface chemistry should be included only with sourced data or
identified uncertainty. Distinguish intact water release from surface dissociation
and bulk-hydrogen outgassing; the latter are alternative or additional sources.

Conserve the inventory: the released amount during a step cannot exceed remaining
coverage. For a constant first-order rate, the exact depletion update uses
DeltaN = N * (1 - exp(-k dt)); do not clip coverage after emitting too many particles.
Use actual face area to convert surface flux to particles or volume sources.
Adsorption, reflection, and recycling must be boundary events with probabilities
per impact, rather than a sticking loss applied once to every occupied cell per step.

For a thermally accommodating surface, inject an outward, flux-weighted Maxwellian
distribution, proportional to vn exp[-m v^2/(2 kB Ts)] for vn > 0, with the correct
normal and flux normalization. Retain three velocity components even if spatial
transport is two-dimensional. Test ballistic transport first, then include collisions
when mean-free-path and collision-time estimates require them. Fragment recoil and
reaction energy partition need their own distributions; do not reset every product
to the surface temperature.

Check the time-of-flight budget at the measured distances and gates. If thermally
released intact water cannot reach a measured region in time, identify whether earlier
release, light fragments, nonthermal desorption, or a different excitation source is
needed. Do not turn a discrepancy into an unreported neutral drift velocity.

## Water breakup and electron kinetics

Here, 'full H2O breakdown' means a documented, state- and process-resolved model with
converged coverage of the relevant conditions. Listing many species is not sufficient.
The core neutral set is H2O, OH, H, O, H2, with O2 and other products retained if
production and sensitivity justify them. Include parent/fragment ions and negative
ions where supported by the channels and conditions. Candidate ions include H2O+,
OH+, O+, H+, H2+, H3O+, H3+, H-, O-, and OH-. Distinguish molecular O2+ from atomic
O2+ written as O^(2+); use unambiguous identifiers in machine-readable data.

The reaction inventory must cover, as relevant:

- elastic momentum transfer and rotational/vibrational energy exchange;
- neutral dissociation into specified products, including H2O -> OH + H and
  subsequent OH breakup, with other branches only where established;
- dissociative excitation into emitting fragment states, including direct H(2p)
  production from water;
- parent and fragment ionization, with secondary-electron energy distributions;
- dissociative attachment, detachment, recombination, charge transfer, and
  ion–molecule chemistry;
- neutral reactions, quenching, radiative cascades, and surface losses or recycling.

Every channel needs provenance, units, threshold, energy/temperature validity range,
uncertainty where supplied, product states/branching, and energy partition. Enforce
H/O-nucleus and charge conservation, with explicit energy and photon budgets. Treat
missing rates as missing; a zero rate is a physical assertion. Stiff time integration
must preserve nonnegative populations without disguising conservation errors by clipping.

Use evaluated electron–H2O data such as Song et al. (2021) as a starting point.
Partial ion yields alone do not specify exclusive complete reaction branches.
NIST BEB total ionization is not a fragmentation or state-population model. A
swarm-fitted effective energy-loss channel must not be assigned literal fragment
or photon yields. Maintain distinct, consistent transport and state-resolved reaction
datasets and prevent duplicate counting of the same collision in both layers.

Reaction coefficients depend on the electron distribution:
k_r = integral[sigma_r(E) v(E) f_E(E) dE], with f_E normalized per unit energy.
A mean electron energy does not determine all these integrals. A local unmagnetized
BOLSIG E/N table is a useful benchmark or conditional closure only when its assumptions
are adequate for the magnetized, transient, nonlocal MITL gap. Surface temperature,
electron temperature, translational temperature, and rotational/vibrational population
parameters are distinct quantities.

For chemistry conditional on measured or prescribed electrons/fields, label that
assumption and do not claim prediction of electrical breakdown. To predict discharge
onset, current loss, or voltage collapse, add a validated kinetic electron/ion model,
self-consistent fields, surface emission, and circuit coupling. FLASH MHD alone does
not supply the magnetically insulated electron sheath. WarpX is a possible PIC/MCC
framework to evaluate; its documented collision models are not a complete water
fragmentation/CR package. If used, rotate coordinates and vectors explicitly: its
2D Cartesian plane is XZ, whereas the convention here is XY.

## Rovibrationally resolved collisional-radiative model

Resolve the H2 X-state populations N(v,J), and the B and C excited electronic-state
rovibrational populations/transitions required for the measured Lyman and Werner
bands. Respect rotational branches, degeneracies, ortho/para nuclear-spin weights,
predissociation where relevant, and radiative cascades. Include electron excitation,
de-excitation, dissociation, and relevant collisional redistribution/quenching. Couple
H2 formation and destruction to the water-fragment network with sourced nascent-state
distributions or explicit uncertainty.

Use a sourced transition list with upper/lower electronic, vibrational, and rotational
labels, wavelengths, and Einstein A values. A Boltzmann rotational distribution can
be a comparison closure but does not constitute a population-resolved CR solution.
Increase state/level cutoffs until the predicted recorded bands converge. Use a
reduced or quasi-steady excited-state treatment only after comparing its timescales
with the evolving chemistry and the camera gate.

See [the state-resolved data requirements](rovibrational_data_requirements.md) for
the MOLAT Lyman/Werner radiative archive and the distinction between a transition
list, a rovibrational corona reduction, and a full time-dependent population solve.

For H, distinguish 2s from 2p. The thesis's unresolved n=2 treatment cannot assign
all n=2 atoms to Lyman-alpha: 2s has two-photon decay and collisional/field mixing,
whereas Lyman-alpha is 2p -> 1s. Track excitation of existing H, direct dissociative
excitation from H2O/OH/H2 where supported, recombination, and cascades separately.
An atomic line intensity alone cannot uniquely determine the parent water inventory.

## Predict the recorded diagnostic

Use emissivity j_lambda = (hc/lambda) A_ul n_u phi_lambda / (4 pi), with a
unit-area spectral profile. Use photon emissivity A_ul n_u/(4 pi) instead when
the calibration is in photons. Integrate over the true viewing volume and camera
gate, then apply wavelength-dependent throughput, gain calibration, spatial response,
line-spread function, and pixel-bin integration. The thesis VUV setup used a 3 ns
gate, covered approximately 115–165 nm, and assumed 5 mm emitting depth in the
original uniform-width inversion. Recover shot-specific settings rather than using
the legacy prototype's 5 ns gate and generic wavelength axis.

Assess optical depth, especially Lyman-alpha. If negligible, document that test;
otherwise include absorption/trapping with dI/dl = j - kappa I and a justified
radiative treatment. Do not equate a plotted transition density with an independently
measured total ground-state neutral density, or adopt the original uniform-depth
density correction without checking the forward model.

The thesis contains calibrated spectral images and lineouts, not the original camera
counts, dark frames, complete per-shot response files, or noise covariance. Recover
curves with source pixel coordinates, axis calibrations and digitization uncertainty.
Retain reported calibration uncertainty separately. Nearby samples along a digitized
curve are correlated and are not independent experimental shots. A fit to these
products supports a reconstruction-based comparison with explicitly limited precision.

## Validation and work order

Before workstation setup:

1. Inventory thesis figures, shot settings, calibration information and contradictory
   labels. Recover selected plotted spectra/profiles and drive traces with provenance.
2. Audit collision and transition datasets by required physical process and state;
   import evaluated tables without inventing missing branches.
3. Specify the matched null, instrument projection, and parameters shared across shots.
4. Validate the SI coupling contract and define analytic/conservation benchmarks.

On Wednesday, record the local FLASH version, OS, CPU/RAM and toolchain. Build and
run a small official verification case, then a planar resistive-diffusion/Joule-heating
benchmark. Measure cost and resolve cold-material/vacuum behavior before a large foil
run. No hardware performance prediction is defensible yet.

Then advance through the following evidence gates:

| Component | Required evidence before publication runs |
| --- | --- |
| Electrode | current/energy conservation; analytic diffusion comparison; material and mesh convergence |
| Surface | finite inventory, face-area flux balance, timestep convergence |
| Transport | ballistic flight-time/angular tests, wall loss and particle-number accounting |
| Chemistry | H/O nuclei, charge and energy budgets; benchmark rates; channel/state completeness |
| CR | known line/band benchmarks; state-cutoff and time-dependent convergence |
| Instrument | photon/energy normalization, gate and LOS integration, calibration uncertainty |
| Inference | null recoverability, sensitivity to priors, held-out prediction, honest uncertainty |

A first manuscript figure set can show: geometry/diagnostic projection; recovered
observations; verified FLASH surface temperature; water/fragment phase-space evolution;
resolved band synthesis and pathway contributions; and matched-null residuals with
uncertainty. Whether the data favor a nonuniform source remains an outcome to compute.

## Primary sources

- [Smith dissertation and source record](https://deepblue.lib.umich.edu/items/9ccfa886-0b19-4bd8-b3b2-bb685b9ddefc),
  especially sections 4.3–4.5 and Appendices A–C; see `dissertation_recovery.md`.
- [FLASH MHD and resistivity](https://flash.rochester.edu/site/flashcode/user_support/flash_ug_devel/node107.html),
  [EOS caveats](https://flash.rochester.edu/site/flashcode/user_support/flash_ug_devel/node114.html),
  [Ohmic heating](https://flash.rochester.edu/site/flashcode/user_support/flash_ug_devel/node128.html),
  [vacuum treatment](https://flash.rochester.edu/site/flashcode/user_support/flash_ug_devel/node146.html),
  [circuit drive](https://flash.rochester.edu/site/flashcode/user_support/flash_ug_devel/node126.html).
- [Bennett et al. 2023, electrode formation and melt](https://doi.org/10.1103/PhysRevAccelBeams.26.040401):
  combined kinetic and MHD description; not a complete water/rovibrational dataset.
- [Song et al. 2021, evaluated electron–water cross sections](https://doi.org/10.1063/5.0035315),
  [author-hosted full text](https://kar.kent.ac.uk/88068/1/5.0035315.pdf).
- [Budde et al. 2022, water electron transport set](https://pure.tue.nl/ws/portalfiles/portal/216541528/Budde_2022_J._Phys._D_Appl._Phys._55_445205.pdf):
  effective transport channels need care when used for chemistry or radiation.
- [Electron-impact water EUV spectrum](https://ntrs.nasa.gov/citations/19850035943),
  [hydrogen 2s/2p experimental distinction](https://pmc.ncbi.nlm.nih.gov/articles/PMC6696584/).
- [WarpX collisions](https://warpx.readthedocs.io/en/latest/theory/multiphysics/collisions.html),
  [collision inputs and coordinates](https://warpx.readthedocs.io/en/latest/usage/parameters.html).

Sources describe capabilities or supply data; the proposed coupled model still
requires implementation and validation in the geometry of these experiments.
