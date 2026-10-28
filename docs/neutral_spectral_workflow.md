# Neutral spectroscopy implementation

The implementation prepares the agreed MFEM Joule → water release → Raikou gap/state
history → observed spectrum → inventory/thermal comparison. It continues
session `01a14f33-a55a-713f-b2f3-1dc6d4cc10c6` from published research revision
`18d155b919e70b0b20165e22ce1eb4b7f48ac5a9`. The retired MHD prototype stays retired.

## Evidence and reconstruction

The [pilot manifest](../data/dissertation/pilot_comparison_manifest.json) preserves
shot 15750's 106 gap-profile samples and 157 wavelength samples, plus 183
wavelength samples for shot 15765. Figure 4.15's partial visible lineout remains
a separate arbitrary-intensity constraint. Original detector files are
unavailable. All VUV readouts retain their literal derived transition-density
ordinate, correlated readout envelopes and omitted support.

The [recovery and scenario contract](dissertation_pilot_comparison.md) separates
dimensional interpretation, reduction mapping and wavelength response. The
41-point 0.01–100 global scale multiplies admissible nominal radiance only.
Current derived-density displays remain at `C_abs=1`; a scale factor cannot
establish their unknown radiometric reduction. No numerical absolute correction
has been inferred from missing calibration records.

```powershell
python tools/validate_pilot_comparison.py --report runs/pilot_readiness.json
```

Evidence validation passes independently of physical readiness. Add
`--require-physical-ready` for a nonzero exit when mapping, clock/collection or
forward-model evidence is still missing. It currently returns exit code 2.

Export the literal retained arrays as an owner prediction request:

```powershell
python tools/prepare_pilot_comparison.py --observation-id 15750_4_20_wavelength --output runs/pilot_15750_request.json
python tools/prepare_pilot_comparison.py --observation-id 15765_4_21_vuv --output runs/pilot_15765_request.json
```

These packets bind array values and source rows to the manifest without inventing
clock endpoints. The bridge's optional `--ensemble-config` rejects substituted
readouts, units, support or unresolved physical admission. The generic comparison
engine's evidence labels alone do not authenticate a dissertation acquisition.
A monolayer convention remains a separate coverage requirement; it does not block
a qualified density/temperature comparison.

## Heating, source and external populations

[MFEM Joule](mfem_joule_heating.md) now supplies the selected heating backend.
The installed MPI MFEM library supports an isolated native miniapp consumer;
FLASH acquisition is no longer required. The stock demonstration is separate
from physical foil qualification. Use the [MFEM surface adapter](mfem_surface_contract.md),
existing [diffusion benchmark](flash_heating_verification.md),
[surface-release](surface_release_pipeline.md) and
[boundary quadrature](boundary_quadrature.md) interfaces after physical setup.

The [Raikou receiver contract](raikou_receiver_contract.md) checks hash-bound
source admission, inventory/energy/momentum accounting, joint moments and
covariance, canonical molecular rovibrational states, atomic H 2s/2p, rate
coverage, electron excitation, velocity closures, clocks and optical admission.
Its manufactured bundle exercises schema and arithmetic. It is not a receipt
from a physical Raikou source-admission execution.

```powershell
python tools/validate_raikou_handoff.py configs/raikou_handoff_manufactured/manifest.json
```

The bundle's simulation cells and collection rays determine optical path weights
once. Finite ray intersections include disjoint support, explicit cell-face
ownership and normalized radiance collection weights. No heating depth, foil
thickness, fitted slab length or detector etendue is inserted.

The [resolved emission projector](../mhdlab/resolved_emission_observation.py)
integrates registered upper-state histories over a contained signed photon gate,
branches each explicit state through its own transitions and integrates
unit-area Gaussian profiles over wavelength bins. It returns line-integrated
and per-nm energy radiance separately. Missing wings stay outside the window.
Only justified optically thin cases are admitted. The current fixed-Gaussian
projector does not itself derive velocity/Doppler or natural profiles from the
checked velocity closures; effective widths and their physical justification
must come from the forward owner. Thermometry requires qualified sensitivity
in that complete forward model.

## Finite ensemble comparison

The [comparison engine](neutral_ensemble_comparison.md) searches complete owner
cases indexed by one parent-water inventory and one declared neutral thermal
coordinate. It never rescales solved populations or interpolates new physical
cases. Each prediction must match the observed quantity, retained coordinates,
support, gate, collection and excitation contract. Mapping can be a qualified
reduction or an explicit conditional reconstruction; unavailable mapping rejects.
The current pilot lacks a usable mapping and registered physical forward cases.

```powershell
python tools/compare_neutral_ensemble.py --config configs/neutral_ensemble_manufactured.json --output-dir runs/new_neutral_comparison --plot
```

The example saves all case/scale alternatives, supplied column/density
diagnostics, discrete compatible parameter values and degeneracy flags. Its
figure is a synthetic recoverability demonstration. Readout compatibility is
the largest asymmetric normalized envelope departure, not a likelihood or
confidence interval. Translational, rotational and vibrational interpretations
remain separate. Derive coverage later using the existing sourced
`MonolayerDefinition`; do not fit coverage and density as separate amplitudes.

[Projected-feature qualification](projected_convergence.md) requires independent
mesh/time/spectral/state refinements below 1% and separate numerical uncertainty
bounds below 10% of any claimed hypothesis contrast. Successive differences alone
are not those bounds. The matched-source null has a synthetic equality check:
fixed linear full-width collection cannot distinguish equal total upper-state
histories simply by redistributing them across width.

## Remaining physical inputs

The local code and manufactured interfaces are reviewable. A physical density,
neutral temperature, monolayer estimate or localization result still requires:

1. A qualified MFEM foil/material/drive setup and
   converged physical surface histories for the selected shots.
2. Owner-produced physical source acceptance and complete resolved transport/CR
   ensembles, with state/rate coverage, excitation and transfer evidence.
3. A defensible scenario-specific physical-to-published reduction, clock
   registration, collection response and instrumental line shape.
4. A sourced areal-density convention before reporting monolayers.

The [implementation record](../supporting_documents/neutral_spectral_comparison/README.md)
keeps commands, checksums, verification and outstanding evidence together.
