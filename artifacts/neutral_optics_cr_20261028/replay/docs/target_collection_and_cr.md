# Collection optics and the water-coverage spectral fit

The primary physical source parameter is now **initial_h2o_monolayers**, a
free nonnegative parameter bound at each spectral forward-model evaluation.
It is not a preset experimental coverage or a multiplier on the historical
four-pool GIF. The new collection operator maps supplied CR upper-state
populations to photons at a declared detector reference plane.

## Confirmed geometry and remaining optical prescription

The user supplied an approximately 19-inch (482.6 mm) target-to-first-mirror
distance, approximately 8-inch (203.2 mm) mirror separation, and a 50-micrometre
by 1-inch-height slit. Slit height is parallel to foil thickness/gap normal y.
The first mirror is a few inches in **straight-line distance** from the slit;
this is not the second-mirror-to-slit optical path.

The [evidence record](target_optics_evidence.md) preserves these user statements
separately from dissertation and manufacturer facts. Exact mirror curvature,
tilts, projected clear pupils, final optical leg, camera magnification and
pixel registration remain unresolved. No earlier unit-magnification geometry
is retained as the actual prescription.

[target_collection_scenarios.json](../configs/target_collection_scenarios.json)
fixes the supplied distances and uses three explicitly assumed final legs of
5, 8 and 11 inches. For an unfolded, equal-effective-focal-length relay focused
at the foil-width center, solving the ABCD focus condition gives:

| Assumed final optical leg | Derived effective focal length per mirror | Source-to-slit magnification |
| --- | --- | --- |
| 127.0 mm | 64.074 mm | 0.15035 |
| 203.2 mm | 73.848 mm | 0.31646 |
| 279.4 mm | 78.583 mm | 0.49705 |

These are focus-compatible parameter scenarios, not recovered mirror part
specifications or a fit of the experimental image. The assumed +x line of
sight puts foil-width edges at different optical depths and projects them
into common y rows. Mirror astigmatism, diffraction, measured aberrations and
the final iCMOS pixel map are not supplied by this paraxial model.

## Photon amount and edge overlap

A full nominal 25.4-mm circular pupil at 482.6 mm subtends **0.00217449 sr**
on axis, collecting **1.7304e-4**, or **0.017304%**, of isotropic point-source
photons before reflection, window losses or downstream clipping. For a unit
1-W monochromatic isotropic source at 121.6 nm over 3 ns, this corresponds to
317779 photons incident on that nominal first pupil. It is a normalization,
not a measured target power or an iCMOS count prediction. Tilted projected
apertures and the actual optical train can reduce it.

`mhdlab.target_optics.build_collection_kernel` traces a finite pupil through
both mirrors, slit and spectrometer f/4.5 acceptance. It includes an opaque
foil shadow, optical-depth defocus and analytic integration along the declared
10-mm invariant current extent. That extent comes from the saved source model;
it is not automatically replaced by the 25-mm physical target height.

The kernel K(row,cell) has units m3 and already contains accepted volume and
the isotropic 1/(4 pi) factor. For resolved transition ul,

```math
N_{r,b}=\sum_{ul} A_{ul}\,T(\lambda_{ul})\,\Phi_{b,ul}
\int_{t_a}^{t_b}\sum_c n_{u,c}(t) K_{r,c}\,dt.
```

T is the separately owned wavelength-dependent photon throughput into the
specified detector plane. Phi is a supplied normalized photon-profile bin
fraction. This separable profile expression needs justification; a measured
spectral/spatial coupling requires a joint response. There is no second
volume, path length, etendue, exposure or 4 pi multiplier.

The [detector photon projector](../mhdlab/detector_emission.py) implements this
expression for supplied state populations and throughput, with exact gate
integration of their declared piecewise-linear interpolant. It requires each
line's actual upper state. Finite spectral windows are not renormalized. It
does not supply an excitation model, line width, detector response or ADU gain.

The present spatial plots contract K with the earlier **illustrative neutral
densities at unit photon rate per neutral**. They test collection/overlap;
they are not emissivities from the new common-water source. For the reference
relay at 57-60 ns:

| Neutral proxy | Near/far geometric-edge profile overlap | Collected fraction originating on side faces |
| --- | --- | --- |
| H | 0.714 | 7.97% |
| H2 | 0.645 | 7.65% |
| H2O | 0.459 | 7.35% |
| OH | 0.462 | 7.34% |

Profile overlap is the sum of minima of two individually normalized row
profiles, between zero and one. The geometric edge region includes the outer
0.5 mm of each width side and its exterior. A separate tagged-packet
calculation attributes photons to the original top/bottom/left/right emitting
face. These definitions are different: material released from a broad face
near an edge is in the edge region without being side-face-origin material.
Full CR excitation can alter every quoted fraction.

The foil preferentially hides the far side. In the equal-rate reference case,
far geometric-edge transmission relative to the unoccluded control is about
50-70%, versus 94-98% on the near edge. Emissivity, optical throughput and
spatial calibration are still required for a physical detector prediction.

Results are in `runs/target_collection_19in_final` and
`runs/target_collection_face_origin`. Images were visually inspected. The
face-tagged densities sum to the original saved densities within 4.9e-17
relative L1 error. Increasing pupil quadrature from 12x48 to 16x64 changes
total coefficients by 1.40-1.45% and row profiles by 1.42-1.57% L1. This is a
sensitivity check, not proof of convergence or optical validation.

```powershell
python tools/compute_target_collection.py --output-dir runs/new_collection
python tools/compute_collection_origin.py --transport-dir runs/neutral_desorption_current_17133_perimeter_v2 --collection-dir runs/new_collection --output-dir runs/new_collection_origin
```

## Common-water source and collisional-radiative model

The [water source binder](../mhdlab/neutral_cr_water_source.py) reuses the
existing parent-inventory and coverage-regression infrastructure:

```math
\Sigma_{\mathrm{H2O},0}=N_{\mathrm{ML}}\,\sigma_{\mathrm{ML}}\,f_{\mathrm{patch}}.
```

N_ML is supplied by the optimizer or a candidate scan. Sigma_ML is an explicit
H2O molecules/m2 definition with an area convention; it is not the Temkin
kinetic reference coverage. Desorption debits one common H2O inventory over
actual face areas. H, H2, OH, O and ions arise from supported conserved
chemistry. The new source has no independent water-derived fragment pools.
Every candidate recomputes depletion; a post-hoc spectrum multiplier is valid
only after an explicit linear-regime qualification.

The [CR model record](neutral_full_cr_model.md) and
[requirements](../configs/neutral_cr_requirements.json) identify the next
state/rate inputs. The existing Raikou native chemistry kernel has validated a
1673-state H2 radiative graph with 28466 sourced branches. A 15-branch subset
was evolved and checked against independent exponential branching, molecule
conservation and photon energy. Its authored initial population is a numerical
control, not water chemistry or a shot simulation.

Full physical CR still requires state-resolved electron excitation and
de-excitation, dissociation/formation pumping, ionization/recombination,
quenching/transfer, electron density/EEDF histories, neutral source internal
states and opacity selection. Ground-state neutral density or foil temperature
cannot substitute for those inputs. The old five-level-H/H2-band proxy is not
used as the full model.

Unresolved throughput can mimic coverage in a linear optically thin regime.
The regression must constrain that response independently or report its
degeneracy. Without a sigma_ML definition, an inventory fit constrains
N_ML*sigma_ML. No unique monolayer result is claimed from these controls.

## SS304 thermodynamics

The [ETHOS/EOS assessment](ss304_ethos_eos.md) identifies SESAME 4273 as a
relevant 304/304L candidate and preserves its published model parameters.
ETHOS supplies a useful multi-fidelity transport-table method; it did not
provide a downloaded SS304 EOS. The new bounded SI p/u table interface and
T(rho,u) inverse are tested, but the actual numeric 4273/electronic free-energy
table is still needed. No unsupported multiphase extrapolation or replacement
of the current heating properties has been installed.
