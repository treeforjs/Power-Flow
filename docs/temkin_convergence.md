# Refining the resolved thermal source

`tools/verify_temkin_planar_convergence.py` studies the manufactured temperature
field

```
T(x,z,t) = 300 + (t/duration)^2 * [250 + 1200*(2*x/W)^2*(0.6+0.8*z/L)] K.
```

The supplied snapshot values are exact spatial bin averages of this field.
Temperatures enter local nonlinear depletion before any area integration.
Width, current direction, original snapshot cadence and frozen-midpoint
subdivisions are refined in separate families. All other resolutions stay
fixed in each family. Refining subdivisions cannot remove the error from
linearly interpolating a coarse set of temperature snapshots.
The default snapshot family holds subdivisions **per snapshot interval** fixed,
so its actual frozen-temperature step also becomes shorter as snapshots are
added. That family's changes include both snapshot interpolation and the
associated integration changes. To isolate snapshot cadence further, use a
common fine integration grid: for example, 8, 16, 32 and 64 snapshot intervals
with respectively 64, 32, 16 and 8 subdivisions give 512 integration intervals
in every case. The fixed-knot integration family separately measures the
subdivision approximation.

Run the manufactured study in a new output directory:

```sh
.venv/bin/python tools/verify_temkin_planar_convergence.py \
  --config configs/temkin_convergence_manufactured.json \
  --output-dir /tmp/temkin-convergence --plot
```

`--families width current snapshot integration event_time` selects the families.
The report records each unique source case once, full resolution tuples, code
and configuration hashes, canonical array hashes, inventory residuals and
finite-reference comparisons. The CSV contains the integrated observables.
The plot omits zero or unavailable relative differences instead of inventing
positive error floors.

The reference in each family is a **finite finer result**, with its resolution
explicitly recorded. Observed orders use three consecutive solutions:
`log(||Qh-Qh/r|| / ||Qh/r-Qh/r²||) / log(r)`. They do not divide errors to the
finite reference. Reversed differences, roundoff-level changes and
nondecreasing differences yield a null order and a warning. A reported order
alone does not demonstrate the continuum limit; inspect the refinement path
and the fixed resolutions too.

`source_observables` integrates actual patch areas, within-patch position
variance and centered mixture variance. The matched uniform-width control
retains the same current-direction, state, time and temperature mixtures.
Its width positions are redistributed over the face. Normal velocity means
use an outgoing flux Maxwellian: `E[vn]=sqrt(pi*kB*T/(2*m))` and
`E[v²]=4*kB*T/m`, averaged over the emitted population. Each state needs an
explicit exact-state mass. The translational energy debit is `2*kB*T` per
released particle. No velocities are rescaled to manufacture an energy match.

Optional whole-source event-time moments use the actual law-specific
conditional quantiles in every active component. They can be expensive; the
default bulk study disables them and runs a separately labelled frozen
component quadrature family. The component is explicitly recorded and is not
an estimate of the whole source. Open Gauss-Legendre probability nodes integrate
actual log-seconds. Centered variance avoids subtracting nearly equal raw
moments. Logs and underflow flags preserve positive moments that round to zero.
`log_moment_comparisons` also retains their relative changes and logarithms of
absolute changes. Raw SI zeros alone do not establish quadrature agreement.
Increasing quadrature order measures a change in estimates, not a certified
error bound or physical uncertainty.

The default manufactured example initially gives release-number orders about
2.18 across width, 1.93 along current, 1.97 for the combined snapshot family
and 2.01 for fixed-knot integration. Some smaller spatial/velocity corrections
are non-asymptotic on those coarse grids; the report retains null orders and
the associated warnings. These values characterize this explicit fixture and
finite refinement paths. They are not accuracy bounds for recovered foil
temperatures or a convergence claim for a coupled experiment model.

All example rates, populations, energy costs and masses are manufactured.
These studies verify prescribed-temperature source calculations. They do not
calibrate water desorption, integrate kinetic gap particles, predict spectra,
apply cooling to the wall, establish moving-wall coupling or run FLASH. A real
heating study needs independent spatial/snapshot convergence, material inputs
qualified across its temperature range and uncertainty from the experimental
reconstruction.
