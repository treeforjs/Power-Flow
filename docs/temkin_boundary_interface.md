# Coverage-dependent thermal boundary interface

The boundary modules turn a qualified Temkin source into weighted thermal
launch packets. They create no gap trajectories, collisions, molecular
transitions or spectra. Heating, source parameters and emitted molecular-state
populations must be supplied with their own evidence.

## Conditional emission times

For a reservoir at frozen temperature, let Nₐ be its initial coverage, N(t) its
remaining coverage, and ΔN = Nₐ − N(Δt) its interval release. Conditional on a
particle having left during the interval, its emission-time CDF is

\[
F(t)=\frac{N_a-N(t)}{N_a-N(\Delta t)}.
\]

For a supplied open probability u, target the population Nₐ − uΔN, then evaluate
its elapsed time using the same separated nonlinear law as the release kernel.
With s = log(Nₐ/N), c = Ed α Nₐ/(kBT Nref), and initial per-particle rate r₀,

\[
t=\frac{1}{r_0}\int_0^s \exp\{c[1-\exp(-v)]\}\,dv.
\]

`FrozenTemkinEmission(reservoir, initial_coverage_m2, temperature_k, dt_s)`
provides `cdf`, `quantile` and `log_quantile`. CDF queries include interval
endpoints; quantile probabilities must lie strictly inside (0,1). It caches a
validated release trial and exposes distinct kinetic and time identities.
The core's c ≤ 512 numerical qualification applies unchanged.

The normalization uses a log released fraction from the mathematical law,
rather than dividing rounded extensive inventories. This preserves positive
tiny releases and meaningful conditioning for subnormal budgets. A stiff,
fully depleted interval can still have finite early emission times; quantiles
are not computed by multiplying a normalized time ratio that already underflowed.
The near-unity-probability calculation retains small survivors through a
complementary form. A certified tiny-hazard limit is conditionally uniform;
general coverage-dependent emission does not use an initial-rate exponential.

The exported log tolerance bounds numerical interval projections. Zero-emission
queries fail explicitly. `quantile` also fails when a positive offset underflows
or its float64 quantization exceeds the timing tolerance; `log_quantile`
preserves the finite logarithm. Keep local offsets separately from the absolute
time origin. An absolute timestamp may lose a small valid offset.

Independent checks compare the CDF to an ODE and quantiles to a separate
population-to-time integration. They cover constant-barrier limits, stiff
intervals, extreme probabilities, reference normalization and timing precision.
These checks qualify the specified frozen-temperature model, not the original
FLASH snapshot cadence or actual wall-temperature history.

## Weighted launches and paired spatial controls

`sample_temkin_emission_packets` takes a distinct `TemkinReleaseHistory`, exact
state masses, one interval, and an n×6 array of supplied open uniform variates.
No random generator is hidden inside the kernel. Its variates select component
refinement, conditional time, width position, normal speed and two tangential
Gaussian coordinates. Each positive reservoir/width component gets at least one
packet; remaining packets refine the component distribution. A component's
physical count is divided by its packet count, preserving rare components even
when their relative probabilities are extremely small. These weights are
unequal and require explicit downstream support.

Normal speed is drawn from the outward flux-weighted Rayleigh law. Both
tangential components are Gaussian at that component's supplied local
temperature and exact-state mass. The expected outgoing translational energy
is 2kBT per particle. Assigned internal energy and the separately supplied
desorption cost are retained without charging the rate barrier as heat.

Resolved and uniform hypotheses share exactly the same time, state, mass,
velocity, momentum and weight arrays. Only their width positions differ.
The uniform control therefore retains the original temperature mixture and
emission-time correlations. Rethermalizing it at one mean temperature would
change the experiment being tested.

Each component retains expected and actual weighted translational energy and
their sampling residual. Finite velocity samples do not exactly reproduce the
analytic mean. A coupled host must balance its actual applied launch energy
against the wall, or use a separately verified moment-conserving representation.
The sampling defect is not another physical heat term. This reference does not
cool a prescribed wall or rescale sampled velocities to hide the difference.

```sh
.venv/bin/python tools/run_temkin_release.py --demo \
  --output-dir runs/temkin_source
.venv/bin/python tools/sample_temkin_boundary.py \
  --source runs/temkin_source/temkin_release.h5 \
  --source-run-summary runs/temkin_source/summary.json \
  --masses configs/manufactured_boundary_masses.json \
  --interval 20 --packets 4096 --output-dir runs/temkin_packets
```

These are manufactured inputs for generic state M, not an H₂O mechanism. The
export contains safe NumPy arrays (`allow_pickle=False`), SI units, exact input
and implementation hashes, source and mass provenance and the new packet schema.
The optional source-run summary must hash-match the exact release file before
its manufactured/sourced/unvalidated parameter label is retained. Without it,
that label stays unreported; it is not inferred from a mass or heating label.
Existing output paths are refused.

## Resolved planar release and 3D geometry

`evolve_temkin_planar_release` retains temperatures and independent reservoir
coverages on every current/width patch. It applies local nonlinear depletion
before area integration. Actual patch area is Δz Δx; there is no second current
length multiplier. Nested extraction geometry, material sampling depth and
explicit boundary-proxy references remain attached to the source.

```sh
.venv/bin/python tools/run_temkin_planar.py --demo --plot \
  --output-dir runs/temkin_planar
```

The default uniform control redistributes only width, separately for every
time/reservoir/current patch. A full-face control is a distinct optional
hypothesis. Each pool uses its own local coverage; this does not infer shared
site competition or automatically preserve a physical layer when it is split
among molecular states. The molecular owner must supply that physical mapping
without duplicating inventory.

The planar schema declares the independent-pool model and numerical
qualifications. Uniform controls use exponent-scaled weighted means; positive
uniform moment densities below float64's normal range are rejected explicitly,
because severe subnormal rounding can change number or energy ratios. Exact
zero remains supported. The bound is a numerical SI representation limit, not
a physical desorption cutoff. The HDF loader preflights local storage, shapes,
units and memory limits before materializing datasets.

`sample_temkin_planar_emission_packets` uses the same qualified time and thermal
laws with an n×7 variate array. It allocates at least one packet per positive
reservoir/current/width component, uses each row's actual Δz once and samples z
inside that patch. Resolved and uniform-width 3D positions retain the same y
face, z coordinate and all other packet properties. Explicit normal coordinate
and sign are checked against raw extraction geometry when available. Current
patch temperatures are a stated local finite-volume approximation, not a claim
of uniform heating across the entire current extent.

The independently runnable example is
[`sample_temkin_planar_boundary.py`](../examples/sample_temkin_planar_boundary.py).
All these interfaces are owner handoffs. No Quark source, build or runner is
needed for their independent verification.
