# Prescribed power and desorption energy

This manufactured reference adds a finite prescribed heating-power schedule to
one stationary patch with a finite parent inventory and supplied constant areal
heat capacity. Its first qualified branch requires **zero activation barrier**.
Number release is then independent of temperature: this branch verifies input
power, thermal evolution and outgoing energy allocation. It supplies no
temperature-dependent water-release prediction.

The [zero-power reference](thermal_feedback_reference.md) and its saved evidence
remain separate. Neither reference installs a source-history, material, FLASH
or host coupling. All coefficients and schedules are explicitly manufactured.

## Equations and channels

For zero activation barrier, the supplied first-order law and wall equation are

```text
Gamma = nu*N
dN/dt = -Gamma
C_A*dT/dt = P(t) - Gamma*(D + I + 2*k_B*T)
```

`P` is areal power in W m⁻² and `C_A` is areal heat capacity in J m⁻² K⁻¹.
The constant removal cost `D` and fixed outgoing internal energy `I` share
an explicitly declared energy reference; neither is inferred from the kinetic
barrier. The term `2*k_B*T` is the expected outgoing thermal-flux translation
per parent. It does not specify fragment energies or an internal-state table.

The cumulative channels remain separate:

```text
dM/dt = Gamma
dE_input/dt = P(t)
dE_binding/dt = D*Gamma
dE_internal/dt = I*Gamma
dE_translation/dt = 2*k_B*T*Gamma
dW/dt = Gamma*(D + I + 2*k_B*T)
```

`W` is the gross wall debit. Its value and gross input energy must remain
positive when the net thermal change vanishes. With `Uads=-D*N`, independent
signed balances include

```text
N + M - N0
C_A*(T-T0) - E_input + W
C_A*(T-T0) + D*(N0-N) + E_internal + E_translation - E_input
```

Binding/internal allocation, the wall-channel sum and bound-reference changes
are retained separately. Patch area converts each areal number, energy and
residual to an extensive quantity once. A total-energy residual alone cannot
verify either trajectory or allocation.

## Finite schedules and numerical evolution

The caller supplies complete constant or piecewise-constant nonnegative power
on the requested interval. Segments must have exact adjacency in the declared
time representation, without gaps, overlaps or an implicit final-value hold.
Each segment supplies its own positive step count and provenance.

Runge–Kutta stages remain inside one power segment. The last stage before a
jump uses the old segment's power; the next segment's first stage uses the new
power. The boundary has one shared state. Splitting at the jump prevents a
single step from averaging incompatible one-sided derivatives.

Logical step times are exact ratios of the canonical segment coordinates and
step counts. The stored `time_s` value is a floating-point display projection;
segment/step indices retain the logical coordinate used for input-energy
support and reference comparisons.

Inventory and temperature evolve as separate numerical states, alongside the
gross energy channels. Exact rational stage products retain independent
positive heating and release terms before projection. Signed thermal increments
may warm, cool or vanish. A nonzero increment that cannot be represented is
rejected rather than rounded deliberately to a balanced state.

Supported temperature/coverage domains and representation limits are explicit.
No temperature floor, inventory clipping, automatic time-step repair or schedule
extrapolation is supplied. Passing those checks is not an accuracy certificate;
temporal refinement and an independent reference are still required.

## Independent zero-barrier reference

The number history is `N(t)=N0*exp(-nu*t)`. Define

```text
A = (D+I)/(2*k_B)
B(t) = 2*k_B*(N0-N(t))/C_A
T(t) = exp[-B(t)]*(T0+A)
       + integral_0^t exp[B(s)-B(t)]*P(s)/C_A ds - A
```

The temperature integral is split at every power discontinuity. Its scaled
exponential kernel is at most one for earlier times. For empty inventory,
`T(t)=T0+E_input(t)/C_A`, with exactly zero release and gas/removal energy.

Input energy is also evaluated directly from segment power times elapsed
support. Translation is independently integrated as
`integral 2*k_B*T_reference(t)*nu*N0*exp(-nu*t) dt`; it is not recovered by
subtracting the thermal residual being tested. Separate number, temperature,
channel and area comparisons can therefore detect incorrect allocation even
when total-energy cancellation looks small.

Nonzero activation barriers are rejected by this first branch. Their rates
depend on temperature, so this exponential/integrating-factor reference would
not verify the forced coupled dynamics. They require a new independent
coupled-IVP reference before qualification.

## Interface and verification

The [module](../mhdlab/thermal_power_reference.py) exposes `PowerSegment`,
`PowerSchedule`, `power_rhs` and `evolve_power`, using an explicitly supplied
`FeedbackProblem` as parameter data. It does not call the earlier reference's
numerical evolution. Library time coordinates are canonical binary64 values;
the command additionally checks exact raw JSON decimal adjacency and complete
support before projection. A gap, overlap or positive duration cannot disappear
silently through rounding.

Run the [manufactured configuration](../configs/thermal_power_manufactured.json)
into a new directory:

```sh
python tools/run_thermal_power_reference.py --config configs/thermal_power_manufactured.json --output-dir runs/thermal_power_reference --plot
```

Four cases cover empty-inventory heating, nonempty net cooling, stronger net
heating and a three-segment on/off/change-power schedule. Seven saved trajectories
include that schedule at 8, 16, 32 and 64 steps **per segment**. The command
retains CSV/JSON histories, code/input/artifact hashes, units, segment/step
indices, separate gross channels and optional PNG/PDF figures. Existing output
destinations are refused; incomplete staging evidence is preserved on failure.

All 37 new focused tests passed on the first run, including an exact
gross-positive/net-zero instantaneous thermal balance. Thirteen malformed
configurations were rejected, including raw decimal gaps/overlaps smaller than
a binary64 spacing and nonzero input values that would convert to zero.

An independent Decimal implementation reads the saved files without importing
or calling the tested core. It evaluates exact exponential inventory,
integrating-factor temperature, direct segment input and positive nested
translation quadrature. It checked 463 saved rows, 7,871 CSV/JSON scalar cells,
6,482 area projections and every logical-time stamp against 292 reused
reference samples. Separate temperature-excursion, emitted-number and energy
channel budgets prevent a large gross input from hiding a small thermal change.

For the schedule at 64 steps per segment, maximum errors over its saved rows,
normalized to the corresponding declared budgets, are approximately:

| Quantity | Maximum absolute error / quantity budget |
| --- | --- |
| Temperature | 6.16 × 10⁻¹¹ |
| Emitted number | 9.66 × 10⁻¹² |
| Outgoing translation | 7.44 × 10⁻¹² |

The largest normalized row error across all grids is about 2.64 × 10⁻⁷ in
the coarsest schedule's temperature. Finite self-order estimates span
4.019–4.041. Changes between 80-digit/16-node and 120-digit/32-node independent
quadratures are retained as numerical estimates, not certified error bounds
or physical uncertainty. The reference uses the explicit binary64 parameter
values, exact logical step ratios and declared SI decimal constants.

In the manufactured scheduled case, 6 J m⁻² of input accompanies approximately
6.171716 J m⁻² of gross wall debit. The signed thermal change is only about
−0.171716 J m⁻², taking the wall from 900 K to approximately 896.565673 K.
Retaining only net thermal change would discard the separately verified gross
transfers and gas energy history. At zero activation barrier, changing power
changes temperature and outgoing translation while leaving the analytic
number history unchanged.

No actual foil coefficients, monolayer convention, molecular state data,
spectral calibration or experimental clock registration are inferred here.
