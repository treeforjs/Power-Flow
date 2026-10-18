# FLASH heating and finite-inventory surface sources

This workflow supplies boundary sources to separately owned gap and molecular
population solvers. It reads prescribed material temperatures, depletes finite
adsorbed inventories, records launch energy, and constructs a matched uniform
release control. It does not run FLASH or predict a spectrum.

## Execute the source workflow

From the repository root with the project dependencies installed:

```bash
python tools/run_surface_release.py --demo --plot
```

The demonstration uses an analytic hot-edge temperature history and two explicitly
manufactured binding reservoirs. Its parameters describe no fitted water surface.
Outputs are written to `runs/surface_release_reference/`:

- `input_surface.h5`: the validated temperature history and its location metadata.
- `surface_release.h5`: interval-integrated release, inventories, energy budgets,
  and matched spatially uniform number/energy moments.
- `integrated_source.csv`: interval and reservoir totals for the declared foil length.
- `summary.json`: parameter/source provenance, hashes, balances and contrasts.
- `overview.png`: temperature, release and depletion diagnostics when `--plot` is used.

For an externally exported boundary-temperature file, supply both inputs:

```bash
python tools/run_surface_release.py --surface /path/to/surface.h5 \
  --config /path/to/physical_reservoirs.json --subdivisions 8 --plot
```

The configuration structure is shown in
[`surface_release_manufactured.json`](../configs/surface_release_manufactured.json).
Supply each reservoir's state ID, attempt frequency, activation barrier, desorption
energy, assigned internal energy and initial coverage explicitly. Scalar initial
coverage is broadcast across width; an array can specify the resolved coverage.
Record a parameter provenance reference and classify its status. Do not relabel
the manufactured example as a physical water mechanism.

## Temperature association is part of the contract

The [surface history contract](flash_surface_contract.md) is now version 2.
It distinguishes a material boundary, a sampled material plane at positive depth,
and an explicitly justified boundary-temperature proxy. Sampling depth and its
justification survive HDF5 round trips. The source integrator rejects an unconverted
material-plane diagnostic.

If evidence supports approximating the boundary with a sampled plane, pass an
explicit `--boundary-proxy-reference` identifying the gradient, mesh-convergence
result or analytic fixture supporting that approximation. This records a model
decision and preserves the original depth; it does not extrapolate a boundary value.

A temperature averaged along a resolved current-direction interval also retains
that interval. Nonlinear release from such a temperature is rejected without an
explicit `--current-uniformity-reference`. Refining only width or time cannot
restore lost current-direction structure. Prefer resolved surface patches and
evaluate local depletion before integrating their released number and energy.

## Depletion law and time accuracy

For reservoir r in width bin i, each frozen-temperature step uses

`h = nu_r * dt * exp[-Ea_r / (kB T_i)]`,

`N_after = N_before * exp(-h)` and
`N_released = N_before * [-expm1(-h)]`.

The implementation evaluates scaled/logarithmic products to preserve representable
small releases and survivors. It imposes no temperature or rate floor. Release
cannot exceed its supplying reservoir. Distinct binding pools can emit the same
gas state without merging their inventories.

Temperature is linearly interpolated between supplied snapshots and sampled at
each integration interval's mathematical midpoint. Depletion is exact at that
frozen temperature. A changing-temperature history still requires temporal
refinement; `--subdivisions` controls this partition. Spatial bin averages must
also resolve the nonlinear source. A sum of independent first-order reservoirs
is the current source model; coverage-dependent barriers, re-adsorption, bulk
diffusion and surface chemistry require additional supplied models.

All inventories use float64. Releases smaller than one spacing of a large
remaining coverage can be represented even when its stored decrement rounds
away. The compensated `inventory_balance_residual_m2` records that limitation.
Do not claim exact cumulative inventory at scales below the stored precision.

## Energy and internal-state ownership

The activation barrier controls the escape rate. It is separate from the energy
cost of removing a particle from its adsorbed reservoir. The recorded wall debit
per released particle is

`desorption_energy + assigned_internal_energy + 2 kB T`.

The `2 kB T` translational mean belongs to an outgoing flux-weighted Maxwellian:
the normal velocity is Rayleigh distributed and the two tangential velocities
are Gaussian. It is not the `3/2 kB T` mean of an unconditioned gas population.
Fragments produced later must receive their reaction recoil from the gap solver.

The desorption energy is referenced to the specified zero-internal-energy gas
reference. Assigned internal energy is measured above that same reference.
Do not add a formation/binding energy again if it is already included in that
cost. These budgets are diagnostics under prescribed temperature; they do not
cool the FLASH material state until a conservative feedback interface is applied.

Reservoir state IDs and energies must match the molecular solver's registry.
The source does not assume a Boltzmann rovibrational population at the wall.
If one physical binding pool is split across initial molecular states, its state
fractions must sum to its original inventory; copying the full coverage into
every level would create particles.

## Uniform-source control

`matched_uniform()` preserves the width-integrated released number and each
energy budget separately for every time interval and reservoir. It redistributes
the source spatially after evaluating the nonuniform heating and depletion.
Applying the desorption law to a mean temperature gives a different experiment,
with a different release rate and energy distribution.

Moment preservation alone does not preserve an entire velocity distribution.
A paired boundary sampler can retain each packet's state, temperature, emission
time, velocity and weight while randomizing only its foil-width launch location.
This retains the original source-temperature mixture. A single Maxwellian at an
averaged temperature generally does not.

The runnable boundary export uses an explicitly supplied state-mass map:

```bash
python tools/sample_surface_boundary.py \
  --source runs/surface_release_reference/surface_release.h5 \
  --masses configs/manufactured_boundary_masses.json \
  --interval 80 --packets 4096 --seed 42
```

The sample mass file is manufactured. Supply molecular masses from the shared
physical state registry for a real case. The exporter writes a pickle-free NPZ
and a JSON contract with paired launch positions, states, velocities, times,
weights, source identity and moment diagnostics. The emission time within a
frozen-temperature interval follows its conditional truncated exponential.

Weighted stratification gives every positive reservoir/width component at least
one packet and conserves its released number through component-specific weights.
At least as many packets as positive components are required. These weights may
differ; an equal-weight collision operator needs a separately validated adapter.
Random thermal velocities conserve the expected `2 kB T` budget in expectation.
The export reports sampled versus expected thermal energy; it does not rescale
velocities to hide sampling error. The paired spatial hypotheses use exactly
the same velocity samples and therefore the same realized thermal energy.

The matched source is a controlled spatial-redistribution hypothesis, not an
independent local desorption/thermal solution. Its assigned uniform launch
locations do not establish a locally balanced wall-energy history. Spectral
distinguishability must be tested after both sources pass through the same gap,
radiative, line-of-sight and detector models.

## Boundary-source HDF5 layout

The output schema is `mhdlab-surface-release-v1`. Arrays `released_m2` and all
energy budgets have shape `(interval,reservoir,width)`. They are integrated over
the interval; divide by its duration for an average flux. `coverage_m2` includes
both initial and final states and has shape `(interval+1,reservoir,width)`.
`temperature_k` holds the frozen interval temperatures.

Reservoir records and source metadata are JSON root attributes. Every dataset
declares units. `matched_uniform/` contains the corresponding spatially uniform
number/energy moments. The loader validates shape, provenance classification,
inventory balance, energy convention and consistency of the derived uniform
moments. Width integration gives particles per current-direction metre or J/m;
multiply by the declared physical length only once.

The legacy reduced runner is unchanged. Use the standalone modules and source
contract for the new FLASH workflow.

## Manufactured convergence study

```bash
python tools/verify_surface_convergence.py --plot
```

The report separately refines the frozen-temperature integration and the width
partition, records all source/configuration hashes and compares reservoir number
and energy totals. Observed orders use three-level solution differences; a
finite finest reference is not treated as an exact solution. Warnings identify
non-asymptotic sequences rather than declaring convergence from a single small
error. Temporal subdivision converges interpolation of the supplied snapshots;
the FLASH output cadence itself also needs independent refinement.

## Resolve both directions on a physical surface

```bash
python tools/run_planar_release.py --demo --plot
```

The planar workflow stores temperature as `(time,current,width)` and coverage
as `(reservoir,current,width)`. Local depletion and energy exchange happen
before patch-area integration. Its default uniform-width control preserves
each current patch separately. A full-face uniform control is available only
as an explicit alternate hypothesis.

The raw-output bridge can export a chosen current-direction partition using
the same configured FLASH layout and physical face:

```python
from tools.export_flash_surface import load_configuration
from mhdlab.flash_planar import extract_planar_surface_history
from mhdlab.planar_release_io import write_planar_surface_history

configuration = load_configuration("/path/to/extraction.json")
surface, manifest = extract_planar_surface_history(
    configuration.snapshots, configuration.layout, configuration.face,
    z_edges_m=[0.0, 0.005, 0.010, 0.015, 0.020, 0.025],
    **configuration.source_metadata,
)
write_planar_surface_history("planar_material_plane.h5", surface)
```

The partition must match the configured physical extent; the numbers above
are an example. The bridge retains each patch's overlap reduction, raw-file
hashes and material sampling depth. It never inflates a previously averaged
temperature trace into resolved data. This reference exporter rescans inputs
for each current patch; a streaming optimization can follow native integration.

External release then requires the physical reservoir configuration and any
explicitly supported boundary approximation:

```bash
python tools/run_planar_release.py --surface planar_material_plane.h5 \
  --config /path/to/physical_reservoirs.json \
  --boundary-proxy-reference "reference to the verified boundary-gradient approximation"
```

Supply scalar initial coverage or an exact `(current,width)` array per reservoir.
Actual z-bin lengths determine areas; a separately configured foil length must
agree and is never applied as an extra multiplier. The planar source output
preserves local patch geometry and all interval/reservoir number and energy
budgets. Refinement in x, z and time remains required for physical predictions.
