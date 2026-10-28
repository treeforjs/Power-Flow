# External plasma model and Yacora comparison

The plasma collisional-radiative model is maintained in
[Raikou](https://github.com/treeforjs/raikou).
Power-Flow supplies heating inputs, conservative surface-release histories,
boundary-source records and experimental evidence. The
[Raikou observation handoff](raikou_observation_handoff.md) records the inspected
source reference, population-output contract and manufactured native readback.
Physical source-receiver acceptance and state/data mapping remain pending.
No physical external-plasma or Yacora comparison is established by that control.

The publication test remains whether width-dependent heating changes water
release and fragment transport enough for the resulting spectra to distinguish
it from uniform desorption. A plasma-model comparison is a separate validation
step on that path. Agreement between two models alone does not establish the
experimental calibration or a physical surface coverage.

## Interface to agree with the plasma owner

Each exchange should identify the source export hash, model/data revisions,
species and state mapping, units and time origin. Preserve full joint emission
time, temperature, parent state, number weight and patch geometry; do not replace
the boundary mixture with one temperature. The source's supplied removal,
internal and expected translational energy channels accompany the records.
Mass, velocity sampling, fragment recoil, charged-particle transport and
state-changing energy exchange are defined by their owners.

The receiver should return a versioned readback of admitted counts and moments,
including time–temperature covariance, area applied once and energy-field
rounding differences. The [boundary quadrature contract](boundary_quadrature.md)
retains exact polynomial, finite-rule and exported-field discrepancies separately.
Source-emission clocks must be registered independently to experimental photon
gates before an observation comparison.

## Staged comparison

1. **Establish common capabilities.** Identify the exact Yacora implementation,
   version, available species, electronic and rovibrational states, reaction
   channels and optical assumptions. Build an explicit overlap table with the
   external model. Do not assume Yacora covers the full water-fragment network.
   Keep unsupported species and states visible.
2. **Check common-data numerical cases.** Where both implementations permit it,
   use an owner-supplied identical state/reaction/radiative subset with the same
   rate coefficients, boundaries and initial populations. Compare transient and
   steady populations, particle balance, energy accounting and line emissivities.
   This isolates numerical/interface differences. It requires no new rate data
   or CR implementation in Power-Flow.
3. **Compare native model predictions.** Freeze each model's own data and
   approximations, then compare the same declared physical cases on their shared
   output basis. Attribute remaining differences to identified data, channels or
   approximations; do not tune a common amplitude to hide them.
4. **Apply one observation model.** Pass both sets of emissivities through the
   same line integration, gate integration, spectral-bin integration, line-spread
   function and instrument response before comparing recorded quantities.

For every case, match electron density and EEDF, background species/state
densities, gas temperature, initial conditions, time-dependent driving,
quenching, radiation field, escape/trapping assumptions and geometry. If Yacora
accepts only a Maxwellian electron temperature for a case, use a matched
Maxwellian subset. A general EEDF with the same mean energy is a different input.
Declare eV/K conversions and electron versus heavy-particle temperatures.

Map states by quantum numbers, energy references, degeneracies and nuclear-spin
conventions. If one model lumps states, compare an explicitly defined sum from
the resolved model and retain the lost resolution. Agreement of a summed band
does not validate all of its rovibrational populations.

## Rotational-band extensions and published benchmarks

Include the published rotational-band extensions as separately identified
comparison models. The following references identify an initial benchmark
family; executable variants, data revisions and benchmark inputs still need to
be supplied. An extension's capabilities must not be attributed automatically
to the standard Yacora configuration.

| Reference | Role in the comparison |
| --- | --- |
| R. C. Bergmayr, D. Wünderlich and U. Fantz (2023), *Molecular data needs for advanced collisional-radiative modelling for hydrogen plasmas*, European Physical Journal D **77**, 136. [DOI 10.1140/epjd/s10053-023-00706-5](https://doi.org/10.1140/epjd/s10053-023-00706-5) | Sections 2.2 and 4.1 distinguish electronic CR, rovibrational Fulcher corona and vibrational ground-state models and discuss their limits. |
| R. C. Bergmayr et al. (2025), *A fully ro-vibrationally resolved corona model for the molecular hydrogen Fulcher-α system*, JQSRT **338**, 109414. [DOI 10.1016/j.jqsrt.2025.109414](https://doi.org/10.1016/j.jqsrt.2025.109414) | Direct rotational-band benchmark: a Yacora-based H₂ corona model for Fulcher-α spectra. Preserve its corona assumptions and state/data version. |
| R. C. Bergmayr et al. (2026), *The plasma physicist’s guide to CR modelling based on the Yacora solver for atomic and molecular hydrogen*, Plasma Sources Science and Technology **35**, 065001. [DOI 10.1088/1361-6595/ae6df8](https://doi.org/10.1088/1361-6595/ae6df8) | Names the distinct hydrogen model variants, including the rovibrational Fulcher model, for a capability crosswalk. |

Bibliographic identities were checked against exact DOI records and the
[official IPP publication list](https://www.ipp.mpg.de/publication-search/4283191?person=%2Fpersons%2Fresource%2Fpersons288648).
The [source summary](../supporting_documents/references/yacora_rotational/summary.json)
records access and scope limits. The official
[H₂ web-model help](https://www.yacora.de/help/help-page-for-h2-model) distinguishes
its electronic-state web configuration from more detailed Yacora-based models;
its indexed public content was available while direct requests returned HTTP 503.
No model execution or full-paper validation was performed for this lookup.

For each paper, record the molecular species, electronic/vibrational band,
retained rotational quantum numbers, spin/parity conventions, branching data
and population assumptions. Compare shared `n(v,J)` populations and resolved
branch/line intensities before comparing convolved band shapes. Keep any
Boltzmann rotational-temperature fit separate from a test of non-equilibrium
rotational populations. An extension for one molecule does not establish
rotational coverage for the full water-fragment network.

Use the paper's declared plasma conditions and instrument resolution. Separate
absolute photon predictions from explicitly normalized shape comparisons.
Digitized published plots remain reconstructed observations with readout and
binning limits. Rotational/internal energy accounting also remains separate
from translational launch temperature; agreement of a band shape alone does
not validate fragment velocities or energy exchange.

Fulcher-α is a visible band. Its rotational benchmark does not establish a
resolved 115–165 nm Lyman/Werner observation model. Treat any UV corona model
as a separate named state/data configuration, with its own overlap and
validation cases. These hydrogen models also do not establish full H₂O,
fragment, molecular-ion or anion coverage.

## Outputs and acceptance

Retain populations in m⁻³ and distinguish photon emissivity in photons/(m³·s)
from energy emissivity in W/m³. State whether emission is integrated over solid
angle or given per steradian. Define photon-emissivity-coefficient normalization
and the densities divided out. Identify transitions, vacuum wavelengths,
branching, spectral density versus bin-integrated quantities and unavailable
outputs. Keep line/band sums reproducible from the contributing transitions.

Before running, freeze case inputs, crosswalks, output definitions and numerical
allowances, including scales for weak lines and differences between models.
Report convergence and balance residuals separately from differences caused by
atomic/molecular data. Experimental uncertainties, calibration scale and
covariance require evidence of their own; numerical tolerances do not supply them.

The first deliverable after the external repository is available is its concrete
receiver contract and one hash-bound readback. Yacora comparison cases follow
when both models' versions, input conventions and overlapping capabilities are
known. Physical coverage regression remains gated by qualified surface inputs
and the [absolute calibration evidence](absolute_spectral_calibration.md).
