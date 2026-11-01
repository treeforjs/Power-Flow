# Native signed thermoelastic exchange, v3

This isolated package forks the sealed native heat component without editing its
live owner. `vendor/pristine/` preserves the exact v1 files and their material
snapshot; `vendor/pristine-manifest.json` identifies the bytes. The independent
working header keeps the native MFEM Q1 grid, mass-lumped nodal enthalpy,
backward-Euler insulated conduction and nonlinear caloric solve. It adds a
dedicated signed mechanical **cell joule** channel. It never puts internal
mechanical work into particle transfers, `authored_heat_J`, or Joule power.

`ss304_heat_material_v2.hpp` binds the copied owner v2 q1 nominal material and its
explicit reference-density convention. `material_v2/` contains unchanged owner
headers. The source seal is
`1a1845e04f1f6fda33338d9bd284f5e1ead9cbc4b7747d0c6355e32cb5661c2d`.
Its caloric support begins at 298.2 K and permits negative enthalpy relative to
300 K. This does not extend the separate solid-mechanics 773.15 K bound or
qualify the material for a particular foil batch.

## C ABI

`thermoelastic_v3.h` is the complete internal C ABI. The frozen MFEM archive is
non-PIC, so it is retained unchanged and linked into a persistent executable,
`ss304-thermoelastic-v3-server`; no dependency rebuild is performed. The Python
owner uses the line protocol below. One opaque handle owns one electrode. Constructor
input includes all canonical nodal temperatures, the full reference geometry,
explicit density, and mesh/material/source identities. There is no inferred
mirroring, support, gas law, or second adsorbed reservoir.

Canonical nodes are row-major `(ny+1,nx+1)`, x fastest; cells are `(ny,nx)`.
The height is included exactly once in all reported energies and volumes.
Native cell energy is `V_cell/4 * sum(h_vol(T_corner))`, the exact partition of
this component's lumped nodal enthalpy. It is not a Gauss integral of nonlinear
`h(T_Q1)`; that would silently change the conserved discrete state.

`prepare` requires the exact accepted revision, clock, native state SHA256,
mesh/material/source/EM identities, and interval-end stage time. The canonical
mechanical sequence is `base_revision+1` and the transaction ID is
`base_state_sha256 + "/mechanical/" + decimal(sequence)`. One cell array covers
the entire reference grid. Positive energy increases caloric energy; negative
energy decreases it. Uniform-cell Q1 load quadrature integrates each cell's J
to its four native corner nodes, with total partition-of-unity checks. Original
cell debits are retained separately from the nodal load.

The C ABI sets external Joule, authored and wall heat to exactly zero. Native
conduction still evolves nonuniform temperatures. Budget fields keep mechanical
exchange separate from all external channels. The local ledger is

`nodal_delta_H_J - nodal_mechanical_load_J + nodal_conduction_out_J = residual_J`.

Do not claim that a cell's reaggregated enthalpy change equals its original
mechanical cell debit: shared nodes and conduction redistribute energy.
The root corrector must combine this exact nodal law with its conservative
mechanical-to-thermal mapping and its independently computed mechanical work.

Prepare leaves accepted state untouched. Preflight/commit bind the prepared
transaction and base digest. Discard permits a corrected trial of the same
unaccepted interval; committed replay fails. Parent orchestration owns the
all-component acceptance barrier. Save refuses pending proposals. Restore
requires a fresh identically configured handle and verifies the v3 checkpoint
magic, payload SHA256, numerical/material/source signature, state dimensions
and replay watermarks. The frozen native binary checkpoint is not a portable
cross-platform interchange format.

## Persistent process protocol

Every request occupies one input line; every response is one `OK` payload line
or one `ERR` message line. Strings use C++ `std::quoted` encoding and numerical
output uses 17 significant digits. The process inherits the parent's CPU7/nice
and frozen environment. No launch or resource-policy override occurs inside it.

- `CREATE`: ABI config fields in header order, with quoted mesh/material/source/
  electrode strings, then nodal count and full temperature array.
- `INFO`: returns revision, time, node count, cell count, reference density,
  stored energy, cumulative mechanical exchange, minimum/maximum supported T,
  and state SHA256, in that order.
- `TEMP trial` and `CELL trial`: accepted (`0`) or prepared (`1`) full arrays.
- `COORD`: interleaved canonical x,y coordinates; `VOL`: nodal reference volumes.
- `PREP`: epoch fields in header order, then cell count and signed cell J array;
  returns the budget fields in header order.
- `LEDGER`: delta-H, mechanical-load, conduction-out, residual arrays, concatenated
  in that order; `MCELL`: original signed cell debit array.
- `PREFLIGHT transaction base_digest`, `COMMIT transaction base_digest`,
  `DISCARD transaction`, `SAVE path`, `LOAD path`, and `EXIT`.

All string arguments above are quoted. No array response has a redundant count;
counts come from `INFO`. Malformed or trailing request fields reject before a
mutation is attempted. `CREATE` requires a fresh server; `EXIT` releases its
native owner and any unaccepted proposal.

Native commit has no throwing cleanup after acceptance, but process/pipe
acknowledgement is not an atomic distributed commit. A lost `COMMIT` response has
an uncertain transport outcome: if the worker remains alive, reconcile `INFO`
and accepted `TEMP 0` against the prepared candidate's expected revision, time,
cumulative mechanical exchange and full nodal temperature array. The candidate
digest is not exposed before acceptance. Do not treat acknowledgement loss as
rollback or blindly replay the commit. If acceptance cannot be proved, parent
orchestration must stop preparing or saving that joint state and recover a fresh
set of owners from a common checkpoint. The acceptance barrier is logical; an
atomic common-checkpoint publication supplies recovery, not a claim that
preflight makes multiple processes crash-atomic.

## Build and tests

Use the existing frozen userspace and dependencies, CPU 7, nice 15, one build
job/thread and the shared `testing-core7.lock`. The parent runs these centrally.

```sh
cmake -S native/mfem_thermoelastic_v3 -B /isolated/build-directory \
  -DMFEM_DIR=/home/treefor/projects/runtimes/mfem-4.10-pic-20261029/install-mpi/lib/cmake/mfem \
  -DCMAKE_BUILD_TYPE=Release
cmake --build /isolated/build-directory -j1
/isolated/build-directory/ss304-thermoelastic-v3-tests /fresh/report-directory
```

The executable tests native nonuniform conduction, signed heating/cooling,
mechanical-cell and nodal ledgers, transactional identity/stage/rollback/replay,
checkpoint equivalence/corruption rejection, and cold/tiny/phase-crossing caloric
increments. These are native component checks; this package does not launch or
qualify a physical coupled water-plasma run, move the thermal grid, solve EM,
or invent mechanical fixture conditions.
