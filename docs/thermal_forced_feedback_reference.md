# Prescribed heating with temperature-dependent release

This separate manufactured reference evolves one stationary wall patch and a
finite parent inventory under a finite prescribed power schedule. A strictly
positive base activation energy makes the release rate depend on temperature;
the optional coverage dependence uses a fixed declared reference inventory.
All material, kinetic and energy coefficients are supplied manufactured data.
They do not constitute a water/OH parameter set or a physical foil closure.

The [zero-power reference](thermal_feedback_reference.md) and
[zero-barrier prescribed-power reference](thermal_power_reference.md) retain
their separate evidence. The present branch requires its own forced coupled
reference; neither earlier trajectory verifies temperature-dependent forced
kinetics.

## Equations and energy ownership

With parent coverage `N` in m⁻² and temperature `T` in K, the declared law is

```text
E_eff(N) = E_a0*(1 - alpha*N/N_ref)
Gamma = nu*N*exp[-E_eff(N)/(k_B*T)]
dN/dt = -Gamma
C_A*dT/dt = P(t) - Gamma*(D + I + 2*k_B*T)
```

`E_a0` is strictly positive and the supported effective barrier is nonnegative.
`N_ref` remains fixed as the inventory depletes and across power discontinuities.
It is not replaced with the instantaneous or segment-start coverage. The
coefficients and their supported temperature/coverage domain are explicit.

`C_A` is a supplied constant areal heat capacity in J m⁻² K⁻¹. `P` is areal
power in W m⁻². The constant removal cost `D` and outgoing internal energy `I`
are declared per parent in a common energy reference. Input eV values convert
to joules using the exact elementary charge before entering these equations.
The kinetic barrier is distinct from either energy debit. A positive barrier
does not require a nonzero `D` or `I` in a manufactured test.

The expected outgoing thermal-flux translation is `2*k_B*T` per emitted
parent. It does not specify fragment velocities, dissociation allocation or a
molecular internal-state distribution. State-resolved chemistry and gap
transport remain separate workstreams.

The numerical state retains independent cumulative channels:

```text
dM/dt = Gamma
dE_input/dt = P(t)
dE_binding/dt = D*Gamma
dE_internal/dt = I*Gamma
dE_translation/dt = 2*k_B*T*Gamma
dW/dt = (D + I + 2*k_B*T)*Gamma
```

The adsorbate energy reference is `Uads=-D*N`. Separate signed balances include

```text
N + M - N0
C_A*(T-T0) - E_input + W
C_A*(T-T0) + D*(N0-N) + E_internal + E_translation - E_input
```

Input energy, gross wall debit and individual gas/removal channels remain
separate even when the signed net wall change is small or zero. Numerical
input is also compared with directly evaluated segment support. Patch area
converts each areal quantity to an extensive quantity once. Residual closure
alone cannot verify a trajectory or its allocation.

## Finite schedules and represented states

The caller supplies complete constant or piecewise-constant nonnegative power
over the requested interval, with explicit provenance and a positive step
count per segment. No gap, overlap, extrapolation or implicit final-value hold
is supplied. Command input checks exact raw JSON decimal adjacency before
projecting coordinates to binary64.

Each Runge–Kutta step remains in one segment. Its last stage uses that
segment's power; the next segment begins with the new power and the same
carried physical state. Logical times use exact ratios of represented
endpoints and step counts. Stored floating-point times are display stamps,
with segment/step indices retaining the logical support.

The parent rate, release amount, energy channels and signed thermal increment
must each fit the supported representation and domain. Unsupported positive
transfers or state changes are rejected. No inventory clipping, imposed
thermal floor, hidden rate cap or automatic step-size repair is applied.
Qualified representation does not supply an accuracy certificate.

Arrhenius suppression is evaluated in a local 50-digit Decimal context before
population and time scaling. Its computed decimal value is retained as a
rational factor through the stage products. This avoids losing a small
positive suppression factor merely because that factor alone would round to
zero in binary64. The exponential remains a numerical approximation.

This branch explicitly supports dimensionless barriers
`0 <= E_eff/(k_B*T) <= 4096`. An exactly zero effective barrier uses factor one;
a positive barrier whose suppression cannot be distinguished from one is
rejected. The barrier bound and representation checks define a numerical
domain, not a physical upper limit or an adopted material law. Instantaneous
derivatives and finite-step transfers are qualified in their own returned
representations.

## Independent transformed reference

The independent reference uses different state variables and high-precision
integration. At a segment start let `N_s` be its own carried inventory, and
define

```text
s = log(N_s/N)
A = (D+I)/(2*k_B)
b = 2*k_B*(N_s-N)/C_A
Z = (T+A)*exp(b)
```

On a segment with constant power,

```text
ds/dt = nu*exp[-E_eff(N)/(k_B*T)]
dZ/dt = (P/C_A)*exp(b)
N = N_s*exp(-s)
T = Z*exp(-b)-A
```

Thus `Z` is constant on an off segment. At a power jump, physical `N`, `T`
and all energy channels carry forward; only the numerical segment reference
is recentered. `N_ref` is never reset. Empty inventory has a separate exact
heating-only branch because its log-depletion coordinate is undefined.

Translation and gross wall debit are separately integrated from their positive
physical integrands along the independent trajectory. Neither is recovered
from a thermal or total-energy residual. Input uses exact logical segment
support; removal/internal costs use independently computed release. Each
channel and the resolved temperature excursion have their own error scale.
Large gross input cannot hide an incorrectly resolved small thermal change.

Precision, method-order and grid changes are numerical estimates. They are
not interval certificates, physical uncertainty or validation of an actual
surface, material phase, detector response or owner implementation.

## Interface and saved validation

The [module](../mhdlab/thermal_forced_feedback_reference.py) exposes
`forced_feedback_rhs` and `evolve_forced_feedback`. It uses the existing
`FeedbackProblem`, reservoir and power-schedule records as parameter data;
it does not call the earlier numerical evolutions.

Run the [manufactured configuration](../configs/thermal_forced_feedback_manufactured.json)
into a new output directory:

```sh
python tools/run_thermal_forced_feedback_reference.py --config configs/thermal_forced_feedback_manufactured.json --output-dir runs/thermal_forced_feedback_reference --plot
```

Four new cases cover constant-barrier heating, coverage-dependent cooling,
positive kinetic barrier with zero removal/internal costs, and a three-segment
on/off/change-power schedule. Seven saved trajectories contain 463 samples.
The first three cases use 32 steps each; the scheduled case uses 8, 16, 32 and
64 steps per segment. Source/input/artifact hashes, declared units, logical
time indices, separate gross channels, residuals and area projections remain
associated with the saved CSV/JSON and figures.

The 32 new core cases were qualified through 33 executions: 31 initial passes,
one failed redundant exact-temperature assertion, and one focused corrected
case. That correction used the same predeclared 5 × 10⁻¹³ K tolerance already
applied to every analytic empty-inventory sample. The production core did not
change and the 31 passing cases were not repeated. Controls include the
small-suppression/large-frequency rate, distinct activation/removal costs,
fixed reference coverage, finite power jumps and representation guards.

Parser evidence retains its own source versions: 17 initial executions and
five focused checks after a preproduction numeric-range guard correction.
The initial parser set was not rerun wholesale on the final command. One
positive four-case batch produced the saved trajectories; verification reads
those files without repeating the batch or earlier suites.

The fresh independent reference uses 80/120 decimal digits, Taylor orders
18/26, two/four subdivisions of each union observation interval and positive
Gauss–Legendre quadrature with 12/18 nodes. It checked all 463 saved rows,
7,871 CSV/JSON numeric cells, 6,482 area projections and 463 logical time
stamps against 292 reused reference samples. Eight completed case/precision
checkpoints preserve the reference work.

One numerical reference pass satisfied all declared budgets. At 64 steps per
segment, maximum scheduled-case errors over the saved rows, divided by their
own reference budgets, are approximately:

| Quantity | Maximum absolute error / quantity budget |
| --- | --- |
| Temperature | 1.63 × 10⁻¹¹ |
| Emitted number | 5.16 × 10⁻¹² |
| Outgoing translation | 5.09 × 10⁻¹² |
| Gross wall debit | 5.13 × 10⁻¹² |

The temperature budget is the resolved 28.787538 K excursion, rather than
the initial 900 K or an input-energy denominator. Across all grids the largest
normalized row error is approximately 6.70 × 10⁻⁸, in the coarsest scheduled
thermal result. Eighteen finite self-order estimates span 4.00268–4.00764.
The largest observed precision/order change is about 4.59 × 10⁻³⁰ of its
declared scale; this remains an estimated reference-stability check.

In the manufactured scheduled case, 6 J m⁻² of input accompanies approximately
4.560623 J m⁻² of gross wall debit and 2.333912 J m⁻² of outgoing translation.
The final signed wall change is +1.439377 J m⁻², giving about 928.787538 K.
Its independently calculated release is approximately 9.26536 × 10¹⁹ parents
m⁻². Release continues through the off segment while the wall cools; the
emitted population therefore retains a time-varying thermal-energy history.

These results qualify the declared manufactured coupled equations and their
numerical range. They install no physical rate set, monolayer inference,
spectral calibration, lateral thermal conduction, MHD heating history,
source-history feedback or accepted host/gap transaction.
