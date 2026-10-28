# Raikou neutral-spectrum observation handoff

This continues **Model MHD Foil Desorption Spectra**, session
`01a14f33-a55a-713f-b2f3-1dc6d4cc10c6`. The target is neutral H/H2 and water-fragment
densities and internal-state distributions, followed by a surface-coverage
comparison between width-dependent MFEM heating and a matched uniform source.
It is not the legacy electron-temperature/four-hydrogen-line inversion.

The original detector data are unavailable. Use the existing dissertation
recoveries with their original quantities and support limitations. The recovered
transition-density displays are not corrected absolute-radiance arrays.

The [implemented comparison workflow](neutral_spectral_workflow.md) now includes
the versioned pilot manifest, additional wavelength readouts, a guarded
[receiver bundle](raikou_receiver_contract.md), ray-cell geometry and finite
[neutral ensemble comparison](neutral_ensemble_comparison.md). These extend the
population-output control below; physical source acceptance remains pending.

## Starting point and ownership

The local starting revision is `18d155b919e70b0b20165e22ce1eb4b7f48ac5a9` on the
published `research/flash-publication-foundation` branch. The supplied session
also records later cloud-local commits through `e356eb2`; those are not claimed
to be installed here. Historical supporting packages and their accepted checks
remain unchanged. The retired reduced MHD runner is not restored.

MFEM Joule supplies heating; desorption and the external gap model supply the evolving
emitting material. Optical path lengths come from that spatial simulation and
the declared collection geometry. Electrode temperature, heating penetration,
foil width, gas extent and line-of-sight length are distinct quantities.

[Raikou](https://github.com/treeforjs/raikou) is now the identified external
runtime. Its locally inspected source reference is
`1997f65374ee9dbdaf2f2e6688af7e28cce008fe`. A source reference and an executable
hash are separate records; a version string does not prove the build matches
that source. No Raikou, Boltzmann, CR, kinetic-gap, Quark or state-resolved data
implementation is changed by this handoff.

## Population-output contract

Raikou's JSON `chemistry.advance` operation returns `success`, `advanced_time_s`
and `state.densities_m3`, keyed by explicit mechanism species/state IDs. Those
densities are particle number per cubic metre, not molar concentrations or
photon emissivities. Native completion and exact state-ID association must be
checked before projection. Per-cell/per-time model output must retain the
mechanism, rates, time convention and physical state mapping supplied by its
owner.

The new [observation projector](../mhdlab/resolved_emission_observation.py)
consumes:

- State-resolved upper-level densities for every saved time and simulation cell.
- Explicit transition IDs, upper-state IDs, vacuum wavelengths and Einstein A.
- Simulation-derived ray-cell lengths with normalized collection weights applied
  once, represented as `path_weights_m`.
- A declared registration of saved-state time to the photon-gate clock, a
  contained rectangular gate, wavelength-bin edges and a supplied Gaussian FWHM.
- Population and geometry provenance, input status, and explicit optically thin
  transfer selection.

It does not infer an emitting length or temperature from the electrode mesh.
The supplied path weights already contain geometry; another fitted slab length
must not multiply the result. They need not have full-width support: disjoint
emitting regions retain their own cells and populations.
Collection weights are dimensionless radiance-averaging weights, not detector
etendue. The caller establishes their normalization, and geometry is fixed over
the selected gate.

For each resolved transition,

```math
I_{ul}(t)=\frac{A_{ul}hc}{4\pi\lambda_{ul}}\sum_c n_u(t,c)w_c,
\qquad
\bar I_{ul}=\frac{1}{t_b-t_a}\int_{t_a}^{t_b}I_{ul}(t)\,dt.
```

Here `w_c` is a collected path weight in metres. The projector exactly
integrates the declared piecewise-linear saved-population interpolant over the
gate; it does not establish accuracy of the underlying CR evolution or clock
registration. No extrapolation or held tail is supplied.

The Gaussian profile has unit area in wavelength, and its CDF difference gives
each bin's captured line fraction. Both bin-integrated W/(m2 sr) and bin-mean
W/(m2 sr nm) results are retained. Finite windows are not renormalized. Geometry,
gate normalization, photon energy, solid angle and bin width each enter once.

This is an explicitly optically thin observation interface. Lyman trapping,
self-reversal or absorbing cells require an owner-supplied ordered transfer
model; independently escaped slab radiances cannot simply be summed. A supplied
Gaussian FWHM does not establish Doppler/Stark thermometry or instrumental
calibration.

## Corrections established by the dissertation audit

Printed p.92 / PDF p.117 gives `A*h*nu*n_upper*L/(4*pi)` in Eq.4.8. Its dimensional
output is line-integrated radiance; per-nm radiance requires a normalized profile.
The printed Lorentz/Gaussian forms in Eqs.4.9/4.11 do not provide that
normalization. The projector uses area-normalized bin integration. The printed
equations do not establish which expressions the original reduction executed.

Printed p.91 / PDF p.116 describes reusing a total molecular band population for
each rotational transition. The new interface instead requires each line's
specific upper-state population. Existing sourced H2 radiative data provide
28,466 transitions, but radiative wavelengths/A values do not supply the missing
excited-state populations. Atomic 2s and 2p also need their own declared mapping;
an aggregate n=2 population is not automatically Ly-alpha's emitting population.

The 5 mm length in Eq.4.8 is foil width. Replacing it with a simulated emitting
distribution changes the forward projection; it is not a detector calibration
correction. Measured radiance must not be rescaled by an MHD thickness.

The existing [absolute-calibration audit](absolute_spectral_calibration.md) and
[manufacturer follow-up](calibration_primary_source_followup.md) remain in force.
The literal single-mm lamp ordinate, original response/gain/exposure reduction
and shot relay throughput are unresolved. No unique absolute experimental
correction factor is determined by the dissertation alone. Preserve calibration
uncertainty using the existing [radiometric response](../mhdlab/radiometric_response.py)
and [coverage regression](surface_coverage_inference.md), rather than assigning
a convenient factor or reinterpreting a derived display as counts.

## Native readback control

```powershell
python tools/probe_raikou_observation.py --executable "C:/path/to/raikou-native.exe" --source-revision 1997f65374ee9dbdaf2f2e6688af7e28cce008fe
```

The probe uses two manufactured decaying levels, checks an actual native
population response, and passes the explicit state IDs into the observation
projector. Its synthetic wavelengths, rates, path and gate are interface inputs,
not a physical H2 spectrum. Saved requests, stdout/stderr, projected outputs,
binary identity and artifact hashes permit readback. A successful probe does
not accept a surface-source receiver, establish radiative energy closure, run
FLASH, correct experimental calibration or infer neutral density/temperature.

`--saved-native-dir` reuses hash-checked requests and outputs from a completed
probe without another native evolution. `--clock-origin-s` explicitly selects
the manufactured photon-clock origin, including negative before-breakdown
coordinates. The projector never guesses or shifts an experimental clock.

The [registration and verification record](../supporting_documents/raikou_observation/registration.json)
binds the external reference and current source hashes. The initial 50 focused
checks passed; after signed-clock support was added, only four clock/partial-gate
cases and ten affected domain guards were run. The current module contains 53
tests. The real native control passed once; a signed-clock readback used its
saved output with zero further native evolution. Unchanged project suites were
not replayed. The exact dissertation PDF was also inspected locally before the
earlier session was supplied; this established no new absolute calibration
factor. Prior completed source/calibration results remain preserved.

The local Python installation has an editable mapping to the retired checkout.
Test commands remove that single mapping only from their own process. The probe
loads this worktree's projector source explicitly, so its evidence cannot be
silently generated by the stale installation. No installed package was changed.

## Remaining owner handoff

The external model must supply the canonical H2O/fragment/ion state map and
rates, including H2 Lyman/Werner excitation, branching/cascades and losses, plus
the evolving spatial populations and absorption needed by the selected case.
Keep translational gas temperature, rotational/vibrational populations and the
electron EEDF distinct. A rotational population fit is not automatically a
translational gas-temperature measurement.

Use the same observation operator, retained bin support and calibration nuisance
model for width-dependent heating and the matched uniform control. Feed those
qualified predictions into a subsequent neutral-state fit and the existing
coverage regression. The
[Yacora comparison plan](external_plasma_yacora_comparison.md) remains separate:
common-data checks precede native-model comparisons, and visible Fulcher
benchmarks do not establish a VUV Lyman/Werner model.
