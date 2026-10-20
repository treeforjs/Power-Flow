# Monolayers, physical inventory and a coverage-source grid

`mhdlab.coverage_inventory` connects an explicitly defined layer coordinate to
the existing prescribed-temperature Temkin source model. It constructs source
histories, not spectra. A downstream owner must supply transport, chemistry,
radiation and detector-bin predictions before using those predictions in
`mhdlab.coverage_regression`. An inventory sweep alone measures no surface
coverage and resolves no absolute calibration ambiguity.

For one physical parent, let `D` be the supplied parent count per square metre
in one monolayer, `L` the nominal number of layers, `q(z,x)` a supplied local
multiplier, and `f_r(z,x)` the fraction assigned to independent pool `r`:

```text
Nparent(L,z,x) = L D q(z,x)
Npool,r(L,z,x) = f_r(z,x) Nparent(L,z,x)
sum_r f_r(z,x) = 1
Nparent,total = sum_z,x Nparent(L,z,x) delta_z delta_x
```

`D`, its counting definition, provenance and status are mandatory. There is no
default water coverage, site density or molecular packing assumption. A kinetic
`reference_coverage_m2` is a separate supplied convention. It remains unchanged
when the initial number of layers changes.

The multiplier is literal. A scalar `1` gives uniform parent density; an array
gives one value per original current/width patch. The code does not divide by a
mean, normalize by area or alter this profile. Consequently `L` is the nominal
layer amplitude. The area-mean number of layers is
`L * sum(area*q)/sum(area)`, which equals `L` only when that weighted multiplier
mean is one. Manifests retain both quantities and the original multiplier.

Pool fractions can also be scalar or patch arrays. They must sum to one at
every patch, including patches with zero multiplier. A tolerance of
`64 * float64 epsilon` accommodates supplied floating-point fractions. Accepted
fractions are retained without normalization; the signed residual
`Nparent - sum(Npool)` and its extensive number residual are reported. Reservoir
names must be unique. Two distinct binding pools may have the same state label
if each receives its own explicit fraction. Assigning the complete parent
inventory separately to every state is rejected.

These fractions partition the *initial parent population*. They are independent
inventory hypotheses, not dissociation branching ratios, daughter yields or
coupled occupancy of shared sites. One assigned parent is counted once. A
species conversion with different stoichiometric counts requires the separately
owned chemistry interface; it cannot be represented by relabelling this
partition. This bounded implementation accepts one parent per sweep.

## Interface

The immutable input records are:

```python
MonolayerDefinition(areal_density_m2, definition, provenance, status)
InventoryPoolAllocation(reservoir, fraction, provenance, status)
ParentSurfaceInventory(
    parent_id, monolayer, allocations, patch_multiplier,
    profile_provenance=..., status=..., provenance=...,
)
```

Each allocation contains an existing, fully specified `TemkinReservoir`. No
rate, energy cost or state property is inferred. Explicit statuses include
`manufactured`, `reconstructed`, `sourced`, `unvalidated`, `owner_supplied`,
`observed` and `source_evaluated`; each component retains its own label and
provenance. A supplied label is a declaration, not independent validation.
Masked, boolean, nonfinite or negative densities/fractions/multipliers are
rejected. Arrays are immutable and retain their original scalar/patch meaning.

```python
plan = prepare_coverage_inventory(
    surface, parent, coverage_ml, subdivisions_per_interval=1,
)
sweep = build_coverage_source_sweep(
    surface, parent, coverage_ml,
    parameter_provenance="explicit source-model association",
    subdivisions_per_interval=1,
)
```

The surface is a validated `PlanarSurfaceHistory` with resolved current/width
patches and an accepted boundary-temperature interpretation. Supplied depth,
proxy and extraction metadata retain the existing source qualifications. The
coverage knots are finite, nonnegative and strictly increasing; a single knot
is allowed, and zero is supported as an empty-inventory control.

The plan holds parent densities `(coverage,current,width)`, pool densities
`(coverage,pool,current,width)`, actual-area parent/pool numbers, allocation
residuals and numerical-domain bounds. `sweep.histories` contains one distinct
`TemkinPlanarReleaseHistory` per supplied knot. `plan.manifest()` and
`sweep.manifest()` provide JSON-compatible associations, full fixed reservoir
records and number/energy ledgers. Direct plan and sweep construction also
validates their arrays and associations. A supplied history cannot replace the
temperature/time mixture, geometry, initial allocation or source metadata.

## Domain and number accounting

The complete inventory grid is checked before any source integration. Each
pool's largest initial density must fit its supplied `max_coverage_m2`. Its
existing reservoir definition also requires a nonnegative kinetic barrier
throughout its declared domain. The numerical Temkin qualification remains

```text
c = Ed/(kB*T) * alpha*N/Nref <= 512.
```

Preflight uses the maximum initial density across coverage knots at each patch
and that patch's minimum supplied temperature over time. Subsequent inventories
decrease and linearly interpolated temperatures stay within their input range,
so this conservatively bounds all later source steps. The bound may reject a
case whose actual midpoint steps would be supported; it does not claim that
such a rejected case is physically impossible. There is no barrier clipping,
rate floor, adjustment of `Nref`, or substitution of constant-rate kinetics.

Time subdivisions must remain distinct at the supplied time origin. Existing
per-history memory limits remain in force. The sweep additionally caps coverage
knots at 1,000 and the aggregate product
`coverage * intervals * pools * current_patches * width_patches` at 20,000,000.
These are computational limits, not physical restrictions or a guarantee that
every permitted allocation fits the available memory.

Density products and actual patch counts are formed with long-double
intermediates. A positive initial parent/pool patch count, density or reported
extensive total that rounds to zero in float64 is rejected. This is an explicit
representation limit, not a minimum physical molecule count. Counts may be
fractional expected populations. No extra extrusion length multiplies the
actual `delta_z * delta_x` area. Source manifests derive number and energy
totals directly from densities and areas before rounding, retain signed balance
residuals, and reject unrepresentable nonzero totals. Existing history and
matched-control numerical qualifications still apply separately.

Source evolution uses the unchanged pool records at every knot and local
temperature before any spatial averaging. Prescribed temperature makes the
wall-energy debit a diagnostic; it does not feed heating back into the surface.
The outgoing mean translational cost is the existing `2*kB*T`; kinetic barriers
and separately supplied binding/internal energy costs retain their distinct
roles. No mass, temperature accommodation or molecular state data is added.

Use each history's existing `matched_uniform_width()` to compare source width
positions while preserving every interval, reservoir/state and current-patch
number/energy mixture. Its full-face control is a different operation. These
controls redistribute the computed source; they do not independently evolve
a newly uniform surface or replace local temperatures with a face average.
The fixed linear-response identifiability control in
[`source_identifiability_controls.md`](source_identifiability_controls.md)
continues to apply before attributing a spectral difference to width structure.

Manufactured sweeps verify this interface and its conservation properties.
Physical layer inference still requires a justified parent-layer convention,
surface allocation, heating and desorption model, owner forward predictions,
and supported observation/calibration uncertainties.

## CLI export and source association

Build the manufactured example into a fresh output directory:

```sh
python tools/build_coverage_sources.py --demo --output-dir coverage_sources --plot
```

`--config configs/coverage_source_manufactured.json` makes the example's
explicit layer, profile, pool allocation and kinetic assumptions visible.
For supplied resolved heating, use `--surface boundary.h5 --config inputs.json`
and set the configuration's `geometry` to `null`; the surface file owns its
geometry. Supplied inputs must retain their own provenance and status. The CLI
requires a new output path, stages its output and applies a stricter interactive
work budget of 200,000 coverage/interval/pool/patch components.

The export contains `input_config.json`, `surface_input.h5`,
`source_0000.h5` and subsequent per-knot histories, `summary.csv`, `report.json`
and, when requested, `coverage_sources.png`. The report associates each source
filename/hash with its coverage knot, physical parent convention, literal
profile and allocations. It also records input, implementation and output
hashes. Keep the report and configuration with the source histories: the
canonical source HDF schema does not itself embed the new monolayer convention
or parent-allocation manifest. An isolated source filename or bare layer
coordinate is insufficient to reconstruct that association.

`summary.csv` has one row per coverage knot and pool. Its
`initial_parent_number` repeats the same physical parent total on each pool row;
do not sum that column across rows. `initial_pool_number` partitions that total,
subject to the reported floating-point allocation residual. Pool release,
remaining inventory and energy columns can be combined across the distinct
pools for the same knot with their stated meanings. Preserve the original
state/pool labels and the full source history when forming downstream controls.
