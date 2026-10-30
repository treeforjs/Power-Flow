# Native thermoelastic/source transaction v4

This package is an isolated derivative of the sealed heat/source v2 snapshot.
`vendor/pristine` preserves all25 sealed files byte-for-byte; the copied
`vendor/source-seal-v2.json` gives the owner paths, hashes and evidence scope.
`source-kernel` and `material-v2` preserve the exact kernel/material source bytes.
The source kernel is compiled once as `ss304-source-kernel-v4`, with its original
C++20 floating-point flags. No v3, live heat owner, or dependency is modified.

The working heat header retains sealed v2's19 cumulative entries, source/water
laws and native MFEM mass-lumped caloric/conduction solve. It adds a separate
signed mechanical cell-J assembly. The wrapper uses the actual sealed
`begin_bound_source`, `advance_bound_source` and `collected_source_heat` calls;
it does not implement a replacement source law or invent source samples.

One native owner holds one disconnected foil. Canonical grids use `(y,x)`, x
fastest; all reported energies/volumes include the invariant height exactly once.
Native cell energy is `Vcell/4 sum(h(Tcorner))`, the exact partition of the native
lumped nodal caloric state. The mechanical solid can have a separate 3D mesh;
conservative mapping and solid energy/work belong to the parent.

The v4 native accepted state includes thermal fields and canonical adsorbed
inventory, all19 heat/source cumulative terms, mechanical energy exchange,
producer/owner watermarks, full typed bound-layer states and surface origins,
particle status, actual neutral packets, canonical gas inventory, weighted source
current paths, and cumulative particle/electrode charge and electron number.
The packet and layer getters expose this same inventory. They are not additional
particles to add to another population. Gas transport and PIC transport are not
implemented in this worker.

## Reference and local energy law

The signed mechanical channel is independent of external wall/source terms:

`deltaH = Qmechanical + Qelectron_wall - Kwater - Uinternal_water - deltaUads`.

Joule and authored heat are exactly zero in this ABI. Certified electron energy
obeys the original source equation
`Qelectron_wall + deltaKe + deltaRe + deltaImage = Wfield`.
Field work, image energy, electron reference and gas/adsorbed reference remain
separate operands. The parent supplies the field/current owner and charges each
weighted path once. Zero-duration quantum bridge paths are retained.

H2O gas reference is zero. Adsorbed reference is the sealed negative reservoir
`-Eb0*q_e*(N-alpha*N*N/(2*N_ML))`; the stable positive removal debit is charged
once. Packet `adsorbed_reference_change_J` is a reservoir diagnostic and never a
transported gas reference. Gas packet kinetic energy and support impulse use the
actual weighted proper velocity. The fixed stationary wall support impulse stays
a diagnostic; no recoil or moving-solid work is silently added to mechanics.

The seven nodal arrays returned by `LEDGER` are deltaH, mechanical load,
conduction out, residual, external load, electron load and water load. They obey
`deltaH - mechanical - external + conduction = residual`, with
`external = electron + water`. Each channel is deposited directly, preserving tiny
source terms independently of larger mechanical loads. The combined solve retains
its declared floating-point residual bound. Original mechanical cell debits are
returned separately; shared nodes and conduction mean that a reaggregated cell
enthalpy change need not equal its original debit.

## Staging and admission

Every source evaluation uses the accepted interval-start patch temperature,
physical area, outward normal and caller-declared one-sided total E/B at that same
accepted stage. Preparing again after discard repeats the same source evaluation;
corrector temperature guesses do not silently redate source sampling.

All field/common/mesh/material/patch IDs, fields, material/planarity/electron-supply
admissions, source quadrature controls, layer policies, physical event weights,
open uniforms and water draws are explicit inputs. Event weights in component
fixtures are manufactured, not predicted electron counts/current. Source
admissions remain conditional model choices. This component does not certify
physical validity merely because a Boolean admission was supplied.
The native strong-magnetic-field fixture explicitly selects the source owner's
0.3 manufactured barrier criterion; it does not replace the owner's nominal0.2
selection or qualify that override for production.

Water-enabled patches must be disjoint on each physical face. Reviewed water
requires300..1300K while its inventory remains; the alternative298.2..3000K law
is an explicit analytical extension. Initial and trial accepted surface states
are checked against the selected range. Source uses the unchanged accepted-start
hazard and coverage law. The native material selects nominal1800 caloric support;
parent mechanics retains its independent298.2..773.15K and2%small-strain bounds.

A new Active layer must be advanced in the same PREP to the interval end; every
previously Active electron must also advance to that common end. Handed-off and
recollected records retain their event time. Active layers are excluded from
ordinary PIC by the parent. Continuation charges only incremental exchange;
birth packet ledgers remain metadata and are never charged again.

Collection accepts an existing handed-off electron, its unchanged physical
weight/work function and actual weighted kinetic energy matching that stored
particle. It tests boundary-event ownership without claiming an intervening PIC
trajectory or inventing field work. General PIC transport needs its own admitted
particle/work adapter. Magnetic recollection through the actual layer is a
separate native path and keeps its displaced material point.

## Internal C ABI and line protocol

`thermoelastic_v4.h` defines the complete ABI. The frozen non-PIC MFEM library is
retained unchanged, so the parent talks to the persistent executable
`ss304-thermoelastic-v4-server`. It inherits the parent's CPU/resource policy.
Every request is one line, strings use C++ `std::quoted`, outputs use17 significant
digits, and the response is one `OK` payload or one `ERR` line.

- `CREATE`: v3 geometry/string/nodal prefix in header order, then patch count,
  each patch's fields in header order, then electron state and parent strings.
- `PREP`: complete epoch in header order, then cell count and signed cell values,
  action count and complete action records, water draw count and draw records.
  Output is the27 budget fields in header order.
- `INFO trial`: trial0 accepted or trial1 candidate; output v3 ten-field prefix
  followed by patch/layer/gas/path counts. The SHA256 covers the complete joint
  state, and the candidate digest must equal the accepted digest after commit.
- `TEMP trial`, `CELL trial`, `VOL`, `COORD`: full canonical native arrays; COORD
  uses interleaved x,y coordinates.
- `LEDGER`: seven full nodal arrays in the order above; `MCELL`: original cell J.
- `LAYERS trial`, `PATHS trial`, `GAS trial`: rows of19,11,20 header-ordered values.
  Layer position is global; origin and outward normal permit independent image
  energy reconstruction. Recollected current kinetic energy is zero.
- `RESERVOIRS trial`: five values per patch: adsorbed Nhi,Nlo, gas Nhi,Nlo,
  negative adsorbed-reference diagnostic. Stable transfer debit is in the budget.
- `CUMULATIVE trial`: sealed19 entries, then mechanical exchange, particle charge,
  electrode charge, physical electron number. Total length23.
- `PREFLIGHT transaction base_digest`, `COMMIT transaction base_digest`,
  `DISCARD transaction`, `SAVE path`, `LOAD path`, `EXIT`.

Arrays have no redundant output counts; use INFO. The joint transaction ID is
`base_joint_sha256 + "/mechanical/" + decimal(base_revision+1)`.

Prepare is read-only for accepted state. Commit preflights the exact joint base,
then performs a noexcept native thermal swap and noexcept auxiliary-state swap;
all allocation precedes acceptance. The process acknowledgement can still be
lost. Parent orchestration must reconcile the exact candidate digest, or mark its
joint state uncertain and recover fresh owners from a common checkpoint. This is
a logical accepted-state barrier, not a distributed crash-atomic protocol.

Native SAVE publishes a checksummed joint payload by temporary-file rename and
rejects pending trials. LOAD requires a fresh identically configured owner,
checks signature, material/kernel/species identities, packet/inventory agreement
and active-layer clocks before acceptance. It also recomputes gas kinetic energy,
proper momentum and source charge/path-work consistency before accepting a
restart. The parent atomically publishes a
common checkpoint manifest including its mechanics and field/current state.
Native binary checkpoints are local frozen-runtime formats, not portable archives.

## Verification

The parent builds and runs centrally in the frozen runtime, CPU7, nice15,
one thread/build job, under the shared testing-core7 lock. Use the unchanged
MFEM_DIR `/home/treefor/projects/runtimes/mfem-4.10-pic-20261029/install-mpi/lib/cmake/mfem`.
Targets are the persistent server, `ss304-thermoelastic-v4-tests`, and
`ss304-sealed-heat-v2-regression`. The latter is the unchanged owner test source
compiled against the working derivative and the exact kernel.

These tests establish component interfaces and discrete ledger closure. They do
not establish source current statistics, moving-surface source validity, whole
field-driver conservation, mechanical support realism, foil-batch calibration,
physical H2O-plasma breakdown or production readiness. No GPU/production launch
occurs here.
