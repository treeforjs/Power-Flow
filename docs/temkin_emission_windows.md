# Partial source-emission windows

`mhdlab.temkin_emission_window.ReducedTemkinEmissionWindow` conditions the
existing frozen-temperature Temkin law on a positive source-emission window
inside its original interval. It does not evolve gap particles or select
photons. A particle emitted before a recorded photon gate can still contribute
photons during that gate after transport, excitation and radiation.

The low-level interface accepts explicit reduced-law coefficients:

```python
window = ReducedTemkinEmissionWindow(
    coupling, log_initial_rate_s, parent_duration_s, start_s, stop_s,
    log_mass_tolerance=1e-7,
)
```

The tolerance shown is an example supplied by the caller, not an experimental
uncertainty. All inputs are required. Coupling lies in `[0,512]`; the initial
log rate is the logarithm of its numerical SI rate in inverse seconds. Times
must satisfy `0 <= start_s < stop_s <= parent_duration_s`. Empty windows are
explicitly unsupported. Masked, boolean, nonfinite and malformed inputs are
rejected. This reduced interface contains no material, molecule, state, mass,
monolayer or calibrated inventory assumptions.

## Law and local coordinates

Let `s=-log(N/N0)`. For initial per-parent rate `r0` and coupling `c`,

```text
ds/dt = r0 exp[-c(1-exp(-s))]
I(s,c) = integral_0^s exp[c(1-exp(-v))] dv
t = I(s,c)/r0.
```

Starting depletion `s_a` is obtained from the original law without the
whole-interval quantile interface's large-depletion sentinel. At local delay
`tau=t-a`, write `s=s_a+delta`. The identical law has

```text
c_a = c exp(-s_a)
log(r_a) = log(r0) - c[-expm1(-s_a)]
I(delta,c_a) = r_a tau
log(window_fraction_of_N0) = -s_a + log[-expm1(-delta_b)].
```

The original population is not restarted or duplicated. The parent-conditional
log probability divides by the original interval's mathematical released
fraction in logarithms. It does not divide rounded inventory differences.
Stable `expm1`/`log1p` forms retain narrow releases and late survival.

`local_duration_s` is the represented original endpoint difference. A
long-double subtraction measures its relative rounding error, which enters
the numerical allowance. Keep original `start_s` and local delays separately.
Adding them can erase valid conditional spread, and taking a logarithm of that
absolute sum can erase it even when individual sums remain distinct.

`cdf(delay_s)` and `log_cdf(delay_s)` use this local delay coordinate, including
both endpoints. The logarithmic CDF at zero is `-inf`. A sufficiently small
positive CDF can round to zero; its logarithm remains available. Open-probability
`log_quantile(u)` returns log local delay; `quantile(u)` requires a positive
float64 delay whose realized log projection error meets the supplied tolerance.
It can reject a coarse subnormal projection while its log quantile remains
qualified. Projection is assessed for the actual requested probability; one
subnormal median can land accurately on the floating-point grid even when
other quantiles do not. Returned arrays are immutable.

## Estimated numerical qualification

This is a bounded float64 implementation using NumPy/SciPy. It does not use a
high-precision fallback. The numerical envelopes assume that sixteen times
SciPy's quadrature error estimate, plus explicit rounding allowances, covers
the integral error. This is an **estimated qualification, not a certified error
bound**. A small estimate alone does not constitute physical validation.

For depletion up to 32, endpoint-scaled positive quadrature evaluates the
separated integral without overflow. For larger depletion the code uses

```text
I(s,c) = exp(c) [s - Jinf(c) + T(s,c)]
Jinf(c) = integral_0^1 [-expm1(-c*t)]/t dt
0 <= T(s,c) <= c exp(-s).
```

`Jinf` is evaluated by positive quadrature, with a bounded alternating series
at small coupling. The omitted positive tail contributes its explicit
`log1p(tail/base)` allowance. This avoids applying ordinary quadrature to a
large-depletion boundary layer too narrow for its nodes to resolve.

Monotonic inversion uses `s <= I(s,c) <= s exp(c)` and
`d log(I)/d log(s) >= 1` to refine estimated brackets in log depletion. Input
log transforms, quadrature estimates, final output rounding and coherent
rebasing uncertainty enter the allowance. In particular,

```text
abs(d log(window mass)/d s_a) <= 1+c.
```

The constructor propagates an *absolute* starting-survival error through that
bound. It also propagates local-depletion and parent-normalization errors.
Large late starts can therefore fail despite a small relative time residual.
The requested `log_mass_tolerance` must cover both reported source-fraction and
parent-conditional log-probability estimates. Smaller requests are rejected;
the code does not enlarge the caller's tolerance.

For quantiles, normalization uncertainty is amplified by the actual conditional
depletion and time-map sensitivities. A probability too close to one can fail
the requested numerical qualification even when construction succeeded.
CDF queries also receive a conditional-log-error check. Their per-query
envelopes are conservative; construction does not promise support for every
representable open probability or local delay.

A rebased positive coupling can underflow. If so, its logarithm is retained,
the zero-coupling approximation is explicitly flagged, and its local log-mass
error is bounded by `c_a`. A positive representable upper allowance and the
finite log allowance are recorded. Rounded zero never denotes an exactly
constant mathematical law. Initial `c=0`, supplied explicitly, is the separate
exact constant-barrier case.

Positive depletion increments below float64 range retain their log coordinate;
the represented depletion can be zero with an explicit flag. Finite log
fractions remain usable when a fraction underflows but a physical count can
be recovered by adding `log(initial_number)` before exponentiation. Conversely,
finite log mass does not establish a positive representable launch budget.
Depletion magnitudes or error allowances above finite float64 range are
unsupported. The numerical domain is separate from the physical model domain.

Malformed inputs raise `ValueError`; a valid mathematical request outside the
estimated numerical qualification raises `TemkinWindowRangeError`, a
`ValueError` subclass. No partial inventory change or host action occurs.

## Metadata and accounting boundary

`numerical_diagnostics` returns a detached JSON-compatible mapping with the
qualification identity/assumptions, requested tolerance, estimated log source
and parent-probability errors, starting-survival error, local/parent log
depletion errors, approximation flags and root evidence. `manifest()` returns
the inputs, reduced coefficients, source logarithms and these diagnostics.
Its JSON representation uses `null` for exact-zero log coordinates such as
`log(c)` at initial `c=0`; the corresponding Python property is `-inf`.
Changing an exported dictionary does not change the immutable window.

The separately qualified physical adapter must retain the actual reservoir,
initial density, temperature, state/pool, patch area, status and provenance.
It must distinguish mathematical log fractions from represented partial
counts, energy diagnostics and balance residuals. An underflowed areal amount
can coexist with a representable extensive amount; multiplying a rounded zero
fraction by a large initial population is unsuitable. Supplied binding/internal
costs and mean outgoing `2*kB*T` remain distinct from the kinetic rate barrier.
This reduced window accepts no conservative host budget or energy transaction.

Verification uses the saved constant-barrier and nonlinear high-precision
references as static inputs, without rerunning their oracle suites. Tests
exercise ordinary and stiff coupling, uncapped late survival, underflow/count
rescue, local timing, projection failures and explicit numerical rejection.
Source timing, photon timing, spatial transport, chemistry and recorded spectra
remain separate modelling responsibilities.

## Physical association and CLI inspection

`mhdlab.temkin_source_window.FrozenTemkinSourceWindow` supplies the physical
association described above:

```python
source_window = FrozenTemkinSourceWindow(
    reservoir, initial_coverage_m2, temperature_k,
    parent_duration_s, start_s, stop_s, patch_area_m2,
    log_mass_tolerance=..., status=..., provenance=...,
)
```

Its initial density, temperature and actual patch area are strictly positive.
It derives the reduced coefficients from that same supplied reservoir and
retains its independent kinetic and energy provenance. Half the supplied
allowance is assigned to the reduced time law and half to represented amounts.
The physical manifest keeps density, actual particle number and each energy
quantity's log value, represented value, qualification reason and rounding
allowance separately. An exactly zero supplied energy cost is distinguished
from a positive energy that underflows. Initial parent density and number must
be representable; partial quantities can remain diagnostic with explicit
unqualified flags.

`all_budget_values_qualified` reports representation status.
`represented_budget()` refuses any unqualified field, retains signed inventory
and energy residuals and always reports `host_budget_accepted: false`. It does
not remove particles, cool a wall, repair residuals or create launch packets.
Expected outgoing translation remains `2*kB*T`, not a sampled particle energy.

Inspect the manufactured component in a new output directory:

```sh
python tools/inspect_temkin_emission_window.py --demo --output-dir source_window --plot --require-represented-budget
```

The default explicit input is
`configs/temkin_window_manufactured.json`. For a canonical planar source, use
`--source source.h5 --config selection.json`; the configuration selects its
original interval, reservoir, current and width indices. Window start/stop
remain offsets inside that original source interval. Supply
`--source-run-summary report.json` to retain a supported source-run report's
exact source hash, parameter status and provenance association. Without that
association the external parameter status must be `unvalidated`.

The export retains exact configuration bytes, selected source and associated
report when supplied, `quantiles.csv`, `report.json` and an optional
`source_window.png`. Its report contains physical/reduced manifests, numerical
and representation qualifications, source geometry/timing context and
input/code/output hashes. Quantiles with qualified finite logs but unsupported
linear projections keep their log delays and an explicit linear-projection
flag. `--require-represented-budget` requests represented fields; it grants no
host transaction acceptance. No instrument photon gate or spectrum is inferred
from the selected source-emission offsets.
