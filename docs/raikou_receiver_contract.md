# Source receipt, population history and collection contract

The v1 [validator](../mhdlab/raikou_handoff.py) checks a saved external receiver
bundle before its population histories enter the
[resolved-emission projector](../mhdlab/resolved_emission_observation.py).
It executes no source integration, transport, CR, FLASH or external runtime.
This is an owner handoff interface. The supplied manufactured control is not an
accepted physical surface source.

## Immutable bundle

The manifest schema is `mhdlab-raikou-handoff-v1`; the complete runnable example
is [manifest.json](../configs/raikou_handoff_manufactured/manifest.json).
Every numerical evidence artifact declares the same `manufactured` or
`owner_supplied` status. Each artifact's exact bytes have a SHA256 in the
manifest. The CLI resolves evidence paths inside the manifest's directory and
rejects paths that escape it. Runtime stdout may be arbitrary nonempty bytes;
the other artifacts are finite JSON without duplicate keys.

| Artifact | Required evidence |
| --- | --- |
| `source` | Positive joint boundary nodes, parent/reservoir/state/query/patch identities, patch area, outward normal and placement, source time and temperature, represented number weights and separate removal/internal/expected translational energies. |
| `receipt` | Bound source and runtime-output hashes, accepted transaction ID, duplicate/replay, rollback and restart evidence, per-component imported moments/covariance, one inventory debit per owning reservoir/patch, actual gas energy and stationary-wall/support ledgers, and patch-placement readback with receiving cell IDs. |
| `states` | Canonical species/state IDs, charge, electronic state, molecular vibrational mode counts/quantum numbers, rotor-appropriate rotational coordinates and nuclear-spin convention; neutral H distinguishes 2s and 2p. |
| `rates` | Bound state map, data revision, resolved transitions and an explicit disposition for excitation, radiative branching, cascades, quenching, predissociation, ionization, recombination and fragmentation. Missing families remain blockers. |
| `populations` | Receipt/state/rate/mesh/electron-history hashes, complete time-by-cell state densities in particles/m3, emitting-state velocity distribution or thermal closure, explicit source-to-photon clock registration and a contained photon gate. |
| `mesh` | Ordered cell IDs, fixed axis-aligned gas-cell bounds in metres and simulation provenance. |
| `collection` | Bound mesh, finite ray segments, already normalized dimensionless radiance-collection weights and optics provenance. |
| `electrons` | Shared source-time/cell coordinates, electron densities, mean energies, EEDF model and excitation-history provenance. |
| `opacity` | Bound population and collection hashes, all retained transitions, assessment method, maximum optical depth, chosen admission limit and optical model. |
| `runtime_output` | Exact receiver output associated with the receipt. The manufactured file explicitly states that no receiver was run. |

An owner-supplied receipt additionally declares the executed command, successful
exit, executable SHA256 and source revision. These identities remain distinct:
an executable hash and a source revision do not constitute a build attestation.
Hashes bind the supplied evidence but do not independently prove its origin,
physical rate completeness, host implementation or scientific validity. The
report always states `independent_runtime_execution_verified=false` because
this command validates saved bytes rather than executing the receiver.

`physical_receiver_receipt_present` means an owner-supplied receipt meeting
this contract is present; it is not a new runtime-verification claim.
Manufactured bundles always set it to false and remain blocked for physical
projection, even when their arithmetic passes.

## Source conservation and the existing quadrature

This interface retains the existing
[positive boundary quadrature](boundary_quadrature.md) and
[host transaction](boundary_transaction_requirements.md) distinctions. An
upstream export is transformed into the compact source view without rerunning
its source evolution. The owner retains that original export and its hash with
the source provenance. For a selected `patch_queries` entry, use these fields:

| Existing quadrature field | Source-view field |
| --- | --- |
| `summary.nodes[].number_weight` | `physical_particle_count` |
| `summary.nodes[].time_s`, `temperature_k` | Same names |
| `expected_translational_energy_j` | `expected_translational_energy_per_parent_j` |
| `desorption_energy_j` | `removal_energy_per_parent_j` |
| `internal_energy_j` | `internal_energy_per_parent_j` |
| Query, patch, source history and state identity | Group identity retained with an explicit reservoir owner and canonical state map |
| `local_time_anchor_exact_s` | Explicit represented `local_time_anchor_s`, with original exact record retained upstream |

The patch-view node weight already includes true area. Multiplying it by patch
area, current length, duration or another number density would duplicate the
source. The source-view moments are recomputed from these actual exported
fields, rather than ideal polynomial targets. The upstream rule-generation,
projection and IVP error qualifications remain separate. Exact-zero source
pieces remain catalogued upstream; active groups here contain positive nodes.
This source view supplies no velocity, molecular mass, fragment recoil or new
thermodynamic reference.

The validator checks each source group independently using the full
unnormalized moment tuple

```text
N, sum(w*tau), sum(w*tau^2), sum(w*T), sum(w*T^2), sum(w*tau*T),
sum(w*removal_energy), sum(w*internal_energy),
sum(w*expected_translational_energy)
```

The receiver must read back the same local time anchor and covariance. A rare
positive group cannot disappear into an aggregate residual. Parent counts use
zero absolute tolerance and a declared relative tolerance. Other quantities
have explicit absolute tolerances in their own physical units plus the
declared relative tolerance; there is no implicit unit-sized floor. The
manifest's tolerances are numerical acceptance choices, not experimental
uncertainty.

Resolved states sharing one physical parent pool have one reservoir/patch
inventory ledger. The total removed inventory and inserted parent number must
equal the sum over those state contributions. Inventory cannot become negative
or be debited once per state.

Source expected translation and actual applied gas translation remain separate.
For the supported stationary, externally held wall, the receipt closes:

```text
wall_total_change = wall_thermal_change + adsorption_reference_change
wall_total_change + actual_gas_translation + gas_internal = 0
wall_thermal_after - wall_thermal_before = wall_thermal_change
external_support_impulse = gas_added_momentum
```

The receipt must retain the actual gas energy; replacing it by `2*kB*T` does
not close a finite launch transaction. The supplied control deliberately has
different actual and expected translation. Moving/recoiling wall accounts and
source reactions that change the admitted parent inventory require a future
contract extension. This validator does not independently reconstruct mass,
nuclei, charge or applied kinetic energy from a receiver particle array. Those
host-budget qualifications remain part of the integrating owner's evidence.

## Gas geometry and signed clocks

`collected_path_weights` derives finite ray-cell intersections from the saved
gas mesh. It rejects cells overlapping in their interiors, zero-length rays,
nonfinite geometry and unnormalized collection weights. Weights must sum to one
within 1e-12; the adapter never normalizes them. A ray parallel to a cell face
uses `[lower,upper)` cell ownership, with the global outer maximum closed, so a
ray on an internal face has one owner. Finite segment clipping, reversed rays,
disjoint emitting support and subdivision of a cell preserve their expected
collected lengths.

The returned `path_weights_m` applies simulation path length and dimensionless
radiance collection once. It contains no detector etendue, foil thickness,
heating penetration, fitted slab length, gate duration or density. The existing
projector consumes these weights directly. Geometry must remain fixed over the
gate; moving cells or curved cells require a different geometric adapter.

Source patch placement is a separate check from the observation ray geometry.
An axis-aligned rectangular patch declares an origin and two oriented span
vectors. Their cross product must agree with its physical area and outward
normal. The receipt reads back that exact placement and identifies receiving
gas cells on the outward face. Their intersected face areas must cover the
whole patch without gaps. Missing placement or readback produces the physical
blocker `physical_source_placement_unverified`; a label and normal alone do not
establish where source material entered the mesh. This check binds declared
support to saved cells; it does not independently verify a receiver's insertion
code or subcell distribution.

`photon_time = source_time + photon_minus_source_time_s` is explicit and may be
negative before breakdown. The adapter does not guess experimental timing,
hold a tail or extrapolate outside the saved history. It rejects an offset so
large that distinct saved times collapse in finite arithmetic.

## Projection and effective Gaussian scope

Validation alone:

```powershell
python tools/validate_raikou_handoff.py configs/raikou_handoff_manufactured/manifest.json
```

Manufactured projection, requiring deliberate opt-in:

```powershell
python tools/validate_raikou_handoff.py configs/raikou_handoff_manufactured/manifest.json --instrument-config configs/raikou_handoff_manufactured/instrument.json --allow-manufactured-control --project-output runs/new_manufactured_projection.json
```

The output path must be new. Without `--project-output`, the projection and its
admission evidence are included in stdout. Both the validator and projector are
loaded from this checkout's exact sources, avoiding the stale editable Python
mapping. Their source hashes and the instrument-file hash accompany the output.
No global Python installation is changed.

The instrument/profile schema `mhdlab-resolved-gaussian-instrument-v1` declares
positive wavelength bins, one effective Gaussian FWHM or one per transition,
the population-history hash, effective-profile justification, velocity
treatment and profile-qualification reference. Its status must match the
bundle. A physical profile must be explicitly owner-qualified for this gate
and geometry. Physical projection rejects any admission blocker; a
manufactured control retains all its missing-physics blockers. Neither can
bypass a requirement for resolved radiative transfer.

The current projector uses supplied effective Gaussian profiles. It does not
fold the checked emitting-state velocities into Doppler line shapes, derive a
measured instrument response, perform Stark broadening, or project non-Gaussian
or time-varying profiles. Owners must establish the effective Gaussian
approximation externally for any physical application. A sample spacing of
0.03 nm/pixel is not that qualification. Translational temperatures in a
velocity closure are not automatically temperatures inferred from spectra;
internal molecular populations remain separate.

Optically thin admission requires a supplied assessment method, a maximum
optical depth no greater than both the owner's limit and **0.01**, and a
declared relative radiance error budget no larger than **1%**. The conservative
fixed-emissivity absorption attenuation proxy `1-exp(-maximum_optical_depth)`
must also fit that budget. A permissive owner limit such as 0.9 cannot bypass
these guards. This proxy does not qualify scattering, excitation feedback or
the owner's optical-depth calculation. Trapping, self-reversal, absorbing
foregrounds or a failed thin assessment require external resolved transfer.
No independent cell escape factors are multiplied into this projection.

The integrated control checks absolute line-energy units against independent
`A*h*c/(4*pi*lambda)` arithmetic, signed gate registration, geometry once,
normalized profiles and bin integration. It does not provide a physical H2
spectrum, calibrated dissertation radiance or neutral density/temperature fit.
