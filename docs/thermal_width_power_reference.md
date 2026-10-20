# Matched-input width heating and finite release

This manufactured comparison distributes the same total prescribed power
uniformly or unevenly across three independent stationary surface patches.
Each patch uses the qualified
[temperature-dependent forced-release reference](thermal_forced_feedback_reference.md).
The initial parent coverage, temperature, areal heat capacity, state identity,
kinetic parameters and energy reference are common to both hypotheses. Only
the spatial assignment of prescribed power changes.

This isolates a nonlinear power-to-release comparison. It supplies no lateral
thermal conduction, field-derived heating, physical surface parameters or
spectral observation model. All geometry and coefficients are manufactured.

## Geometry and exact input matching

The declared transverse edges are −1/512, −1/1024, +1/1024 and +1/512 m,
with common length 1/128 m. The three areas are 1/131072, 1/65536 and
1/131072 m², giving total area 1/32768 m² and exact weights 1/4, 1/2, 1/4.
This constructed geometry is not associated with the experimental foil.

| Hypothesis | Left patch | Middle patch | Right patch |
| --- | ---: | ---: | ---: |
| Uniform power, W m⁻² | 60,000,000 | 60,000,000 | 60,000,000 |
| Nonuniform power, W m⁻² | 90,000,000 | 30,000,000 | 90,000,000 |

For patch area `A_i` and prescribed areal power `P_i`, both satisfy

```text
P_total = sum_i A_i*P_i = 234375/128 W
```

They share the same finite duration, nominally 200 ns. Exact rational
arithmetic on the represented coordinates, areas, powers and duration checks
the common power and input energy. No fitted normalization is introduced.
The displayed common energy is approximately 0.000366211 J; its exact
represented-input value remains in the saved metadata.

## Histories, reuse and area

There are only three distinct forcing histories: uniform, edge and middle
power. Each is evolved once at unit reference area for 64 steps, retaining
65 samples. The uniform history supplies all three uniform patches; the edge
history supplies both nonuniform edge patches. Reuse retains separate patch
geometry and source identity without repeating evolution.

For any per-area number or energy quantity `q_i`, patch and total quantities are

```text
Q_i = A_i*q_i
Q_total = sum_i Q_i
```

The adapter applies the actual patch area to areal fields once. It does not
rescale the scalar history's already extensive fields. Geometry, parent/state
parameters, energy reference, logical time grid and history associations are
checked before assembling the comparison.
Temperature remains an intensive state variable; it is not multiplied by area.

The saved data retain each patch's coverage, temperature, instantaneous parent
flux, emitted number and cumulative input, removal, internal, translation and
gross wall-debit channels. Signed thermal change, adsorbate-reference changes
and numerical residuals remain separate. Instantaneous flux is evaluated once
per unique history sample and reused with that history. The same call saves
all right-hand-side fields, including the temperature derivative in K s⁻¹;
retaining these fields requires no additional evolution or flux evaluation.

Patch temperatures and their time dependence remain explicit. An area- or
number-weighted summary temperature does not replace the full source mixture
or define a single outgoing velocity distribution.

## Comparison and interpretation

Signed differences use **nonuniform minus uniform** on the same finite time
grid. Equal total input does not impose equal released parent number or equal
outgoing translation: the local activation barrier, inventory and wall
temperature evolve together. The sign and magnitude of any contrast must
come from the saved calculation and independent reference.

Each hypothesis retains the separate wall and gas balances. Because the
initial inventory and thermal state match, differences in final wall energy,
adsorbate-reference change and outgoing energy must also account for the
equal input. A small balance residual alone cannot verify a release contrast.

Independent checks use fresh transformed-IVP and positive-channel references
for the three new forcing cases, then aggregate patch quantities independently
with exact geometry. Repeated patch assignments reuse the same independent
reference. Difference checks use resolved contrast scales separately from
gross inventory or input-energy scales. Initial structural zeros and possible
near-zero differences remain explicit.

This comparison qualifies a manufactured upstream source calculation. It
does not establish actual foil heating, monolayer coverage, H₂O fragmentation,
Hall transport, absolute VUV calibration, spectral discrimination, a gap
solution or owner acceptance.

A resolved source contrast still needs a state-resolved transport/emission
prediction and a justified observation map before the recorded spectra can
test the two hypotheses. Source counts and translation energy alone do not
supply that discrimination.

## Interface and reproduction

`assemble_width_power_reference` accepts three supplied typed
`ForcedFeedbackHistory` records, exact geometric support, explicit reuse maps,
canonical history hashes and manufactured provenance. It performs admission,
cached right-hand-side evaluation and area assembly; it calls no evolution
routine. Its result is recursively immutable, and `to_dict()` produces a
detached JSON-compatible tree.

The command below creates the three distinct scalar histories and their saved
assembly. Use a fresh output directory; an existing destination is rejected.

```sh
.venv/bin/python tools/run_thermal_width_power_reference.py \
  --config configs/thermal_width_power_manufactured.json \
  --output-dir /tmp/mitl-width-power-new --plot
```

`width_assembly.json` preserves complete histories and exact extensive
arithmetic. Three scalar CSVs retain cached derivatives and state identity;
`patch_samples.csv`, `aggregate_samples.csv` and `width_comparison.csv` expose
the corresponding assembled records. Exact rational values and signed
projection diagnostics accompany projected floats. File hashes, dependency
hashes, input identity, reuse and scope remain in `report.json`.

There are 195 unique scalar samples, 390 patch rows, 130 aggregate rows and
65 comparison rows. Patch copies are spatial assignments of the same saved
source history, not additional trajectories. This bookkeeping is separate
from independent reference evaluations and synthetic admission tests.

## Saved manufactured comparison

At the final nominal 200 ns sample, both hypotheses display cumulative input
0.0003662109375 J. The uniform patch temperature is 1036.838925 K;
the nonuniform edge and middle temperatures are 1143.302031 K and
930.237778 K. These are constructed-case temperatures.

| Extensive final quantity | Uniform | Nonuniform | Nonuniform minus uniform |
| --- | ---: | ---: | ---: |
| Released parents | 3.1290170084603115 × 10¹⁵ | 3.120527579286983 × 10¹⁵ | −8.4894291733285 × 10¹² |
| Expected outgoing translation, J | 8.221273856823092 × 10⁻⁵ | 8.252207951974316 × 10⁻⁵ | +3.093409515122285 × 10⁻⁷ |
| Gross wall debit, J | 1.574113076433863 × 10⁻⁴ | 1.5751662511903592 × 10⁻⁴ | +1.0531747564961036 × 10⁻⁷ |
| Wall thermal change, J | 2.0879962985661355 × 10⁻⁴ | 2.0869431238096467 × 10⁻⁴ | −1.0531747564889208 × 10⁻⁷ |

For these inputs, uneven heating releases fewer parents while their total
outgoing translation is larger. Local hot-edge release alone does not fix the
sign of the total parent-number difference.

The 34 new admission/assembly tests passed using explicitly synthetic histories
and mocked right-hand-side records. The three actual histories were then
generated in one separate batch. A fresh independent transformed-IVP reference
used 80/120 decimal digits, Taylor orders 18/26, two/four subdivisions and
positive Gauss–Legendre quadrature orders 12/18. It checked all saved scalar
states, eight cached rate fields, patch assignments, aggregate records and
signed differences. Old case grids and passing suites were not replayed.

Across the sampled trajectory, the largest error divided by each quantity's
own maximum resolved width contrast was 7.225 × 10⁻⁸. For parent release it
was 2.358 × 10⁻⁸ and for translation 2.708 × 10⁻⁸, all below the predeclared
2 × 10⁻⁵ allowance. The largest reference-ladder change on those contrast
scales was 6.707 × 10⁻³⁴. Initial structural zeros remain exact; pointwise
near-zero contrasts are qualified separately. These finite comparisons are
numerical error estimates, not certified bounds or physical uncertainty.

The saved numerical run used the pre-correction plotting source. A subsequent
plot-only change removes unsupported zero-baseline closure at outer patch
edges. Its original numerical report, data and figures remain immutable;
corrected figures and the exact plotting-only source change are saved
separately. No scalar evolution or rate evaluation was repeated for that fix.
