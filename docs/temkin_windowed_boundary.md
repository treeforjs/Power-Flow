# Boundary launches from partial source-emission windows

The windowed source interface conditions emission within the original planar
source intervals. It preserves each independent reservoir, exact supplied state,
temperature, current/width patch and actual area. It does not define a photon
recording gate: particles emitted before a diagnostic gate can still contribute
later through transport, excitation and radiation.

`SourceWindowSelection(interval_index, start_s, stop_s)` identifies local offsets
inside one original source interval. Adjacent selections are allowed; overlapping
or repeated selections of an original interval are rejected before conditioning.
The explicit input order determines component order. Selections in different
intervals keep their different original inventories and frozen temperatures.
No remaining population is transferred between independent reservoirs or states.

```python
source = condition_temkin_planar_source(
    history, selections,
    log_mass_tolerance=1e-7,
    status="manufactured",
    provenance="Numerical source-emission test with explicit local selections",
)
packets = sample_temkin_windowed_emission_packets(
    source, state_masses_kg, uniforms,
    face_normal_coordinate_m=0.0, outward_normal_sign=1,
)
```

The source functions are in `mhdlab.temkin_windowed_source` and the sampler is in
`mhdlab.temkin_windowed_boundary`. The canonical history is revalidated against
its existing law, energy and geometry contracts. This validation does not evolve
a new temperature or surface history. Every positive original pool inventory
creates a conditioned component. Exact zero initial inventories have explicit
skip records. A positive component is retained even if a density underflows;
the sampler evaluates whether its physical number and energy can be represented
after incorporating the actual patch area. Finite logarithms alone do not
establish a launch budget. Unsupported components cause explicit rejection.

The reduced scalar law retains uncapped starting survival and computes partial
mass without subtracting nearly equal endpoint CDFs. Its numerical error
envelopes use estimated quadrature and floating-point allowances, rather than
certified interval bounds or experimental errors. A constructor can succeed
while a particular requested conditional quantile fails its separate numerical
qualification. The sampler requires a qualified positive local linear delay and
also retains the corresponding log delay.

The seven strictly open uniform columns select component refinement, conditional
local emission time, width position, outward normal-speed CDF, tangential
Gaussian radius, tangential angle and current position. One mandatory packet
represents each component. Remaining packets refine the component mixture, with
unequal particle weights preserving each component's represented inventory
within the reported floating-point residuals. The interface supplies all draws
explicitly and contains no random-number generator. The export command records
its explicit seed, generator, draws and runtime versions.

The paired comparison changes launch width position only. Every pair shares
the original component, local timing, current and normal position, exact state,
mass, temperature, particle weight, velocity and momentum. The redistributed
width positions sample the complete foil width. The comparison retains the full
temperature/state/velocity mixture; it does not rethermalize at an average wall
temperature. Fixed linear responses independent of launch width must recover
the paired full-width equality before a width-sensitive mechanism is tested.

Mean outgoing translational energy is `2*kB*T` for the flux-weighted Maxwellian.
Sampled translational energy is a separate finite-sample quantity and is never
forced to equal that expectation by rescaling velocities. Supplied internal and
wall removal energies remain attached to the exact reservoir, separately from
the kinetic activation barrier. Signed packet-number and sampled-energy
residuals are retained. A coupled host must specify its actual launch-energy and
inventory transaction policy. `host_budget_accepted` remains false throughout
this interface; no wall cooling, state branching or daughter yield is inferred.

Original interval origin, window start, local delay and log delay have distinct
meanings. A rounded absolute timestamp can merge with a boundary or lose a small
delay. Keep the separate coordinates and inspect the exported precision flags.
Position offsets similarly retain their local meaning if adding a large origin
rounds a coordinate onto a patch edge. Numerical source selection does not align
an experimental timing origin or recover a missing diagnostic trigger.

Export a canonical source into a new directory:

```sh
python tools/sample_temkin_source_windows.py --source source.h5 --config windows.json --source-run-summary source_report.json --output-dir window_packets
```

The configuration explicitly supplies selections, exact-state masses and their
separate status/provenance, numerical tolerance, stationary face geometry, packet
count and seed. A supported run summary must hash-associate the exact source,
source provenance and declared parameter status. Without that association,
external source parameters require status `unvalidated`. A declared status
preserves provenance and does not independently validate physical inputs.
The provided `configs/temkin_windowed_packets_manufactured.json` applies to the
previous manufactured coverage-source fixture; it supplies no physical water
or foil default.

The export contains exact configuration bytes, source and optional run summary,
uniform draws, `boundary_packets.npz` and `report.json`. NPZ files contain plain
numeric/string arrays and JSON metadata, and can be loaded with
`numpy.load(..., allow_pickle=False)`. Reports retain per-array units, component
associations, numerical qualifications, signed accounting residuals and
input/code/output hashes. Existing paths, including broken symlinks, are refused.
Staging guards ordinary single-writer failures; it is not a concurrent filesystem
transaction protocol.

Only upstream boundary injection is represented. No kinetic gap, charged-particle
trajectory, collision/radiation solve, molecular state data, thermal feedback or
actual FLASH run is supplied by this interface. Physical surface coverage and
absolute spectral normalization still require qualified experimental and
material inputs and separately supplied forward predictions.
