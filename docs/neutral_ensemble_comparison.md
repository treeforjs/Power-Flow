# Neutral forward-ensemble comparison

`mhdlab/neutral_spectral_inference.py` compares complete, externally produced
neutral forward cases with bounded curve readouts. It does not solve chemistry,
transport, heating or excitation. It does not turn an upper-state population
into a ground-neutral density by assuming a closure.

Run the manufactured example from the repository root:

```powershell
python tools/compare_neutral_ensemble.py --config configs/neutral_ensemble_manufactured.json --output-dir runs/new_neutral_comparison --plot
```

The output directory must be new. The command saves the exact input, all
case/scale scores in CSV, a JSON report, input/source/output SHA256 hashes and an
optional figure. It explicitly loads the implementation from its own checkout,
so an editable installation pointing at another checkout cannot select an old
implementation. The fixture is a numerical recoverability example, not a
dissertation fit or physical spectrum.

## Input contract

The JSON schema identifier is `mhdlab-neutral-ensemble-config-v1`. Its required
top-level fields are `schema`, `observation` and `cases`. Optional fields are
`calibration_scales` and a manufactured-only `manufactured_truth` record. Truth
records are saved for independent validation and never used by the comparison.
Unknown constructor fields and duplicate JSON keys are rejected.

An `EnsembleObservation` declares a categorical `reconstruction_id`, a list of
`ReadoutEnvelope` curves, status, mapping status, signed gate endpoints,
`time_registration`, `observation_contract`, `excitation_contract`, calibration
provenance and observation provenance. Mapping status is one of:

- `qualified`: the supplied mapping is declared qualified. Admission still
  requires independent review of its evidence.
- `scenario_assumption`: a named conditional reconstruction. A result remains
  conditional on that assumption.
- `unavailable`: the comparison rejects inferential use until a defensible
  published-to-forward mapping is supplied.

Each curve contains these fields:

| Field | Meaning |
| --- | --- |
| `curve_id` | Unique observable identifier |
| `coordinate`, `coordinate_units` | Strictly increasing retained abscissas and their literal units |
| `center`, `lower`, `upper` | Readout center and asymmetric bounds, with `lower <= center <= upper` |
| `support_segment_ids` | Nonnegative integer per retained point; contiguous runs explicitly separate missing intervals |
| `quantity`, `signal_units` | Literal published or reconstructed physical quantity and units |
| `scale_role` | `nominal_energy_radiance`, `nominal_photon_radiance` or `published_derived_unscaled` |
| `provenance` | Readout/reconstruction source and limitations |

Radiance roles require energy- or photon-radiance quantity identifiers:
`spectral_energy_radiance`, `integrated_energy_radiance`,
`spectral_photon_radiance` or `integrated_photon_radiance`, respectively. Signal
unit strings must exactly match each prediction. The interface performs no
implicit dimensional conversion. A derived transition-density ordinate cannot
be relabeled as radiance or a ground-neutral density.

Coordinate-readout bounds, source pixels, glyph masks and other reconstruction
metadata remain in the evidence manifest. They are not automatically converted
into stochastic errors. A forward producer must evaluate the declared
observable at each retained nominal coordinate under an explicit sampling
prescription. This interface neither shifts wavelengths within readout bounds
nor fills missing segments.

Each `NeutralForwardCase` supplies:

- A unique case ID, hypothesis and reconstruction ID; one
  `parent_water_areal_inventory_m2`; one `thermal_parameter_K`; and a shared
  `thermal_coordinate` and `thermal_coordinate_definition`.
- Complete `ProjectedCurve` values with exactly the same curve IDs,
  coordinates, support-segment IDs, quantities and units as the observation.
- The same gate, clock registration, observation prescription and excitation
  contract as the observation. Different electron assumptions require a
  separate declared comparison.
- A source receipt; population and transport provenance; `conservation_status`
  equal to `passed` and its provenance; optical-admission and numerical-
  qualification evidence; and overall provenance.
- Canonically keyed `neutral_columns_m2` and
  `weighted_neutral_densities_m3`, with the same neutral species IDs, plus a
  `neutral_diagnostic_definition` stating collection/gate weighting and density
  normalization.

Manufactured observations and physical cases cannot be mixed. Source,
conservation, opacity and numerical records are supplied attestations; matching
labels do not independently verify a forward model. The producer must regenerate
population, transport and energy histories for each physical inventory/thermal
case. This layer never scales a solved population to manufacture a new neutral
temperature or density.

Physical clock, observation, excitation, source, transport, population, optical
and numerical-evidence fields reject exact missing-evidence sentinels such as
`unknown`, `unavailable`, `pending`, `TBD` and `N/A`, including case/spacing and
simple punctuation variants. A conditional reconstruction does not waive clock
registration. This check establishes only that an attestation is present: it
cannot authenticate a free-text claim or an invented numeric gate. Physical
dissertation admission additionally requires binding the observation to the
evidence manifest, literal recovered data and concrete mapping/time/collection
evidence. Descriptive manufactured records remain numerical fixtures.

## Score and calibration scenarios

For prediction `p`, readout center `c`, lower bound `l` and upper bound `u`, the
point score is `abs(p-c)/(u-c)` above the center and `abs(p-c)/(c-l)` below it.
The case score is the maximum over every retained point and every curve. A score
of at most one means the prediction is inside all readout envelopes. Repeating
correlated figure samples does not add statistical weight. A zero half-width
admits exact equality; displacement toward that side is unbounded. Such scores
are saved as JSON `null` with an explicit `unbounded_departure` flag, not as NaN
or an invented tolerance.

By default, comparisons containing nominal radiance use the agreed 41-point
scan `10**linspace(-2, 2, 41)`. `C_abs` multiplies the nominal **observation's**
center and bounds to obtain scenario-corrected radiance. Forward spectra and
their populations remain unchanged. This is an exploratory sensitivity range,
not a calibration confidence interval, distribution or unit correction.

Only radiance-role curves receive this multiplier. Derived published displays
remain unscaled. A derived-display-only comparison admits only `C_abs=1` and
cannot acquire an absolute-radiance interpretation from this scan. Categorical
unit/reduction scenarios are separate input reconstructions, not additional
continuous fit parameters.

## Interpretation and validation

The report retains every case/scale row, all tied supplied minima, all compatible
case IDs, and discrete compatible inventory/thermal values. It reports supplied
neutral column/density diagnostics with their definitions. It neither interpolates
between forward cases nor extrapolates outside the finite ensemble. Compatible
values are conditional alternatives, not confidence intervals. If no case is
compatible, the report says so; it does not synthesize an off-grid best fit.

The declared thermal coordinate must remain explicit. A source translational
temperature, gate-weighted translational temperature, rotational distribution
descriptor and vibrational descriptor are distinct quantities. A good spectral
match does not establish their equality. Surface coverage can subsequently be
derived from the parent inventory through an independently sourced monolayer
definition; density and coverage are not separate amplitude parameters here.

Tests cover asymmetric envelopes, calibration direction, exact constraints,
unit/clock/gate/support rejection, retained temperature/inventory alternatives,
ground-neutral versus excited-state degeneracy, and the matched-source null
through the resolved-emission projector. Equal total upper-state histories with
equal full-width collection produce equal spectra under that fixed linear
observation. Hypothesis discrimination therefore needs justified changes in
collection, transport or excitation rather than an amplitude artifact.
