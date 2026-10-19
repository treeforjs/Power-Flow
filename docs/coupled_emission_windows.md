# Coupled emission windows from saved surface histories

This manufactured boundary interface turns accepted constant-power surface
histories into finite emission windows. It preserves the joint emission-time,
parent-state and surface-temperature measure needed by a later transport
handoff. It supplies no velocity sampler, gap evolution or spectral model.

The inputs are the three saved histories from the
[matched-input width comparison](thermal_width_power_reference.md). Their
original numerical report, raw files, canonical history identities and actual
independent acceptance are explicit inputs. The exporter evaluates no surface
ODE, source right-hand side or width assembly. Plotting-only source changes
remain separate from the source version that generated those histories.

## Positive interpolation on finite support

Let `M` denote cumulative emitted parents per unit area. In a saved cell of
exact duration `h`, define local coordinate `u=(t-t0)/h`. The cubic Bernstein
controls for emitted number are

```text
[M0, M0 + h*Gamma0/3, M1 - h*Gamma1/3, M1]
```

The temperature controls are

```text
[T0, T0 + h*Tdot0/3, T1 - h*Tdot1/3, T1]
```

The endpoint values and derivatives come from the same saved source identities.
Nondecreasing emission controls provide a conservative sufficient condition
for nonnegative emission throughout the cell. Temperature controls must lie
inside the supplied positive temperature domain. Unsupported controls reject;
there is no slope limiter, clipping, constant-temperature fallback or source
refinement during export.

A requested subwindow is intersected with every supported source cell. Exact
rational de Casteljau subdivision retains the restricted temperature and
emission controls and their local time coordinates. The restricted derivative
already includes the subinterval scale. Its integrals receive no second time
factor. No endpoint padding, extrapolation or held tail is supplied.

The interpolated measure uses stored `M`. Stored coverage `N`, its independent
inventory residual and the original source channels remain attached. The
interface does not redefine `N=N_initial-M` or claim that the interpolant
solves the nonlinear release law between source samples.

## Number, joint moments and energy

Each window uses a declared local time anchor `a`, with `tau=t-a`. Direct
polynomial integrals against the same positive measure `dM` retain

```text
integral [1, tau, tau^2, T, T^2, tau*T] dM
```

The complete restricted coefficients also remain available. Means, variances
and time–temperature covariance are statistics of this retained joint measure.
A mean temperature does not replace its components or define a single thermal
velocity distribution. Centered time statistics keep the exact anchor separate
from rounded absolute-time displays.

Expected outgoing translation is computed directly:

```text
X = 2*kB*integral T dM
B = D*DeltaM
I_out = I*DeltaM
W = B + I_out + X
```

Here `D` and `I` are the separately supplied removal and internal energies;
the kinetic activation barrier is not debited. The exact SI Boltzmann constant
and elementary charge define the unit conversions. These scalar manufactured
inputs do not provide molecular-state branching, fragment recoil, particle
masses or rotational/vibrational energy pathways.

Exactly zero emitted mass is explicit. Its normalized statistics are
unavailable, with a reason, rather than zero or a pseudocount. Nonzero exact
quantities whose float projection underflows or overflows reject. Exact
rational values and signed projection diagnostics accompany finite displays.

## Mixtures and comparisons

There are three declared windows for each of three distinct saved sources:
complete support, an early knot-aligned interior window and an off-knot
interior window. Nine unique source/window results supply 18 patch references,
six hypothesis/window aggregates and three nonuniform-minus-uniform contrasts.
The patch references preserve actual manufactured areas, state identity,
source hashes and the full time/temperature pieces.

Patch area multiplies each areal unnormalized quantity once. Aggregates sum
those exact extensive moments before normalizing by total emitted number.
They do not area-average patch means or variances. Local time anchors must
agree before moments are combined. Signed extensive contrasts subtract the
unprojected aggregate values; differences between separately normalized
statistics have a distinct meaning.

Equal total prescribed input across the heating hypotheses does not impose
equal nonlinear emission or energy. It does not establish the separate fixed
linear full-width response equality used in owner response controls.

## Three distinct error layers

Exact polynomial algebra qualifies the declared interpolation and its finite
representation. It does not establish continuous kinetic accuracy.

At original-knot endpoints, direct polynomial translation can be compared
with the difference of saved cumulative translation. Source inventory,
internal allocation, input accumulation and wall-channel residuals remain
separate from that interpolation discrepancy. For off-knot endpoints, saved
cumulative-channel comparisons are unavailable unless an additional
representation has been explicitly declared.

The polynomial energy defect is retained as

```text
R_poly = C_A*(T(b)-T(a)) + (D+I)*DeltaM + X - P*(b-a)
```

It is not forced to zero. Knot-aligned comparisons retain the separate source
energy defect and the contribution of each source allocation/input discrepancy.
No energy channel is inferred as a remainder to repair the polynomial measure.

Finally, new independent moment quadrature and off-knot queries compare with
the underlying coupled IVP. Accepted independent knot states seed those new
calculations; unchanged full-state grids are not repeated. Numerical reference
stability, source integration error and interpolation error are qualified
separately. Finite numerical comparisons do not provide certified enclosures,
physical uncertainty or spectroscopic significance.

Only constant-power histories are supported. Scheduled power jumps require
one-sided temperature derivatives and a separate admission contract. No
physical foil, H2O coverage, trigger/gate registration, VUV calibration,
state-resolved chemistry, Hall transport, FLASH run or owner acceptance is
established by this manufactured interface.

## Interface and input identity

`admit_saved_width_source` consumes ten explicitly named byte inputs and their
expected SHA-256 values: the source report, assembly, input configuration,
reuse map, three raw histories, independent completion/report and scientific
review. It verifies byte identity before parsing and checking the associations.
File location, copying and path safety belong to the command-line loader.

`CoupledEmissionCell` and `WindowQuery` declare exact finite support and
manufactured provenance. `integrate_emission_window` computes a single supplied
polynomial measure. `build_coupled_emission_windows` calculates the nine unique
source/window pairs and assembles their patch references and aggregates. The
result is immutable; `to_dict()` returns detached JSON data.

The manufactured query configuration fixes fractions of the exact source
duration: `[0,1]`, `[1/16,5/16]` and `[23/128,71/128]`. The final query has
off-knot bounds. These fractions refer to source emission time; they are not
registered experimental gates.

```sh
.venv/bin/python tools/export_coupled_emission_windows.py \
  --config configs/coupled_emission_windows_manufactured.json \
  --source-dir /path/to/accepted-width-pass \
  --output-dir /tmp/mitl-coupled-windows-new
```

The destination must be new. `coupled_emission_windows.json` retains the complete
joint measures and comparisons. Four CSVs expose source windows, patch windows,
aggregates and width contrasts; `polynomial_pieces.csv` retains every restricted
piece. Exact values, projections, availability and reasons remain explicit,
including literal `null` for absent CSV quantities. `report.json` binds code,
configuration, copied inputs and result hashes. The copied `source-inputs/`
directory preserves the required source layout and original bytes for review;
these are reused manufactured inputs, not new experimental evidence.

## Saved numerical qualification

All 50 new tests passed once with warnings treated as errors. Their exact
polynomial fixtures qualify the new mathematics and admission guards; they
are not surface-kinetic trajectories. The sole separate saved-data export
produced nine source windows, 18 patch references, six aggregates, three
contrasts and 315 polynomial pieces, with zero source evolution, source-rate
or old assembly calls.

Fresh independent checks formed polynomial targets in an ordinary power basis,
then performed new positive moment quadrature and off-knot queries from the
accepted 80/120-digit independent knot states. The two levels used Taylor
orders 18/26, two/four local subdivisions and positive quadrature orders 12/18.
There were 384 seeded cells, 1152 local substeps, 18,432 node-flux evaluations
and 12 new off-knot endpoint state pairs. No old full-state grid was repeated.

The new interpolation allowance was fixed at 2 × 10⁻⁴ before calculation.
The largest observed scaled window error was 1.656 × 10⁻⁷, the largest gross
mixture error was 2.980 × 10⁻⁸, and the largest error on a width contrast's
own resolved scale was 9.710 × 10⁻⁸. The largest width-contrast reference-ladder
change on those scales was 8.715 × 10⁻³⁴. These are finite estimates with
explicit scales, separate from the earlier scalar-history qualification.

The exported extensive contrasts, **nonuniform minus uniform**, are:

| Source-duration window | Parent-number contrast | Translation contrast, J |
| --- | ---: | ---: |
| `[0,1]` | −8.48943 × 10¹² | +3.09341 × 10⁻⁷ |
| `[1/16,5/16]` | +2.87616 × 10¹¹ | +4.67155 × 10⁻⁸ |
| `[23/128,71/128]` | −3.64858 × 10¹¹ | +1.38115 × 10⁻⁷ |

The independent reference resolves these signs numerically. The early-window
parent contrast has a different sign from the complete-support contrast in
this constructed case. This is a source-emission result; no experimental gate
or spectral discrimination is inferred.

Over complete support, the exported emitted-number weighted mean temperatures
are 951.520 K and 957.699 K for the uniform and nonuniform hypotheses. Their
temperature variances are 1572.38 K² and 5055.68 K². The complete coefficients,
time statistics and covariance remain attached to these summaries.

Six source/window comparisons have both bounds at saved knots. The three
off-knot comparisons explicitly lack stored cumulative-channel differences.
Independent readback checked 351 CSV rows, 17,781 cells, 7875 piece
coordinate/control associations, 756 quantity/statistic records and 1078
source-defect records. No interpolation, inventory or energy residual was
removed to achieve these checks.
