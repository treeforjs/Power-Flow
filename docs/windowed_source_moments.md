# Deterministic moments of window-conditioned boundary sources

The local-time and joint interfaces describe the same original source-emission
selections as the windowed packet interface. They retain initial population,
state, frozen temperature, current/width patch and actual area. Source-emission
selections do not define photon recording gates or recover diagnostic timing.

`mhdlab.window_time_moments.conditional_window_time_moments(window,
quadrature_order=...)` integrates the `ReducedTemkinEmissionWindow` local-delay
law at open probability-space Gauss-Legendre nodes. Log quantiles are anchored in
actual local seconds. A centered scaled calculation forms variance without
subtracting two large raw moments. Neither the original interval origin nor the
window start is added to these conditional local delays.

The immutable `WindowTimeMoments` record retains the local mean, raw second
moment and variance, their logarithms, SI underflow flags and quadrature order.
Positive subnormal SI values are rounded projections; the logs retain the
estimates before that projection. Unsupported log quantiles or unresolved
relative quantile spread cause explicit rejection. A source-window constructor
can succeed while a requested quadrature node fails its separate numerical
qualification.

`refine_window_time_moments` computes orders `n` and `2n`, with
refined-minus-coarse changes. `compare_window_time_moments` compares supplied
estimates with estimate-minus-reference signs. Changes use log moments even
when both corresponding SI projections underflow. The caller must associate
compared records with the intended common source law and selection. Finite-order
changes are resolution indicators, rather than certified error bounds or
physical uncertainty. A small change between two orders alone is not proof of
accuracy.

`mhdlab.windowed_joint_source_moments.windowed_joint_source_moments(source,
quadrature_order=..., include_component_records=True)` combines original
conditioned components. It returns resolved and matched-width summaries,
conditionals by selection and reservoir, and optional component records with
their local timing and exact source associations. Exact rational arithmetic on
represented component counts, coordinates and moments prevents a tiny displayed
probability or a large coordinate origin from silently deleting a component.
It does not make the underlying quadrature or physical parameters exact.

Joint time means retain a fixed anchor chosen from the dominant-count
component's original interval origin and window start. The reported mean offset
uses that anchor at every quadrature order. Original origins, window starts and
local means remain separate before centered covariances are formed. Best-effort
global mean time has an independent precision flag. Spatial moments integrate
the supplied uniform launch rectangles analytically. The width-square reference
is the center of the full original face for every conditional population.

The matched-width summary changes only launch width position. It retains the
joint time/current/state/temperature/energy mixture. Fixed responses independent
of launch width therefore preserve their full-width means. Width-dependent
covariances and apertures can distinguish the launch hypotheses. This control
does not independently deplete a uniform-temperature surface and does not
establish distinguishability of recorded spectra.

Conditional mean outgoing translational energy is `2*kB*T`. Two different
variances are reported:

```text
Var(component conditional mean energy) = Var(2*kB*T)
Var(full emitted kinetic energy)       = Var(2*kB*T) + E[2*(kB*T)^2]
```

The second term is thermal scatter inside each flux-weighted Maxwellian, whose
kinetic energy has a gamma distribution with shape two. At a fixed component,
local time, launch position and thermal kinetic energy are independent. Their
mixture can have nonzero cross-covariances because components carry different
times, temperatures and geometries. Supplied internal and wall removal costs
remain separate source ledger quantities, rather than kinetic activation
barriers or additional outgoing translational energy.

The joint interface admits individually qualified physical counts and requires
supported local moment projections for its normalized statistics. It rejects
unresolved positive local timing projections and nonzero final statistics that
cannot satisfy their representation allowance. It does not substitute zero for
an unsupported positive component. Empty source populations have zero number
and null normalized statistics. The primitive log-moment interface remains
useful for diagnosing cases outside the joint SI representation domain.

Inspect a saved packet export into a new directory:

```sh
python tools/inspect_windowed_source_moments.py --packet-export saved_packets --orders 8 16 32 --output-dir window_moments
```

The command checks the saved report, exact input files, producer code pins,
packet array schema/units/shapes, status/provenance and canonical source
selection association. It performs
canonical read/validation and new moment calculations. It does not rerun the
original source evolution or regenerate packets. The original export files are
copied with explicit prefixes, each requested moment report is saved, and
`report.json` binds input/code/output hashes and finite-order changes. Existing
paths, including broken symlinks, are refused. Staging handles ordinary
single-writer failures; it is not a concurrent filesystem transaction protocol.
Structural admission and hash association do not independently recompute each
saved sample's timing, velocity or realized energy. Those checks can use the
separately recorded packet realization and independent readback.

Sampling qualification must respect the actual fixed component allocation. For
represented component fraction `alpha_c` and packet count `n_c`, the weighted
mean's conditional sampling covariance is

```text
Cov(mean estimate) = sum_c alpha_c^2 * Cov_c / n_c.
```

Between-component population spread is retained in the source distribution; it
is not extra random allocation variance once the component counts are
conditioned on their realized values. A plug-in weighted covariance has
expectation equal to population covariance minus this conditional covariance of
the estimated means. An empirical within-component correction needs at least
two packets in each relevant component. A singleton component instead needs an
independently supplied conditional covariance or an explicit unsupported flag.
A single global weighted Bessel factor generally does not correct distinct
source components. Uncertainty of a covariance estimate needs the appropriate
mixed higher moments or justified conditional-independence formulas; second
moments alone do not supply every such uncertainty.

These numerical comparisons concern source injection. They accept no host
transaction, wall feedback, kinetic gap evolution, charged trajectory,
collision/radiation calculation or molecular state data. Manufactured examples
retain their labels. Physical water coverage, experimental timing, absolute
VUV normalization and physical spectral predictions remain separate inputs.
