# Joint source moments and matched width controls

The source comparison needs correlations as well as integrated release. A
symmetrically heated foil can have `Cov(time,x)=0` even when its edges release
earlier than its center. The joint interface therefore includes
`Cov(time,(x-face_center)^2)` and
`Cov(expected_translation,(x-face_center)^2)`, as well as time/current and
time/thermal-energy covariances. They describe the upstream source, not a
recorded spectrum or transported gas distribution.

```sh
.venv/bin/python tools/inspect_joint_source_moments.py --demo \
  --output-dir /tmp/joint-source --quadrature-order 8 --plot
```

The example uses manufactured temperatures, independent surface populations,
rate parameters and bound-removal costs. To inspect a qualified source, use
`--source path/to/temkin_planar_release.h5`. Its parameter status is unreported
unless `--source-run-summary` supplies a hash-associated source-run record.
Renaming a file does not invalidate an association with its exact bytes.
Source provenance and any supplied reservoir definitions must agree too.

`joint_source_moments` integrates actual component counts over physical patch
areas. It uses conditional emission-time means and centered variances from the
identified coverage-dependent law. Uniform positions within each patch give
exact within-patch position moments, including fourth moments needed for the
variance of squared width. Width always refers to the full face center.
Original width bins and integration intervals retain conditional statistics;
empty populations have zero count and null normalized statistics.

Interval starts and local emission offsets remain separate until centering.
Selected-population anchors and centered spatial polynomials retain small local
correlations beside large common time/position offsets. Absolute reported means
still have float64 coordinate precision; the face-center endpoint parts are
retained in the reference metadata. Nonzero statistic products that round to
zero fail explicitly rather than creating a structural zero. The timing
kernel's separately flagged underflow/log moments remain attached. These are
numerical range qualifications, not physical release or energy cutoffs.

The uniform-width control redistributes those same component populations over
width. It preserves current position, molecular state, time and local thermal
mixture. Its width correlations disappear without replacing the temperature
mixture by one mean temperature. The hypotheses are alternative sources, not
two populations to inject together. The molecular owner must map physical
adsorbate inventory to states without duplicating coverage among independent
pools.

The thermal-energy variable here is each component's expected outgoing
translation, `2*kB*T`. Its variance is the variance of these local conditional
means. It excludes the kinetic-energy spread of the flux Maxwellian within one
component. It is also distinct from sampled packet energy, molecular internal
energy, bound-removal cost and the kinetic activation barrier. Conditional
mean timing and energy covariances therefore do not add a new wall debit.

The output contains full joint metadata in `report.json`, conditional CSVs,
the exact source HDF5, input/implementation/output hashes and an optional
figure. The demo also retains its canonical surface and input configuration.
The factory's numerical and source-temperature qualifications remain in force.
Quadrature order needs refinement independently of source space/time grids.
No arbitrary diagnostic-gate clipping is performed: a gate crossing an original
interval needs a separately qualified conditional truncation operation.

These correlations help choose source observables that retain a heating
signature. Establishing whether measured spectra distinguish the hypotheses
still requires the separately owned gap transport, chemistry, radiation,
diagnostic response and uncertainty analysis. The interface runs no FLASH or
production host and applies no wall feedback or motion.

The [identifiability control](source_identifiability_controls.md) states a useful
limit: a fixed linear response independent of launch width, with full-width
weighting and appropriate boundary conditions, gives identical observations
for these matched sources. Width-dependent geometry/transport/weighting or a
source-dependent nonlinear response can break that limit. A latent source
correlation alone does not demonstrate spectral discrimination.
