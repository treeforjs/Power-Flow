Native field-mesh integration — isolated electrostatic mode

Baseline: immutable supplied native-manufactured-source-v2.tar.gz, identified
in baseline-manifest.json. The sealed heating driver/native provider and its
benchmark header are separately copied/pinned. No active runtime/shared source
or reaction/neutral implementation was edited. Rebuild all native consumers;
GapPic's C++ class layout changed. Raw Particle/Options/Diagnostics sizes remain
80/168/200 bytes. The checkpoint format deliberately advances to PIC3.

GapPic(options) retains a uniform field mesh. GapPic(options,field_mesh) adds
PicFieldMeshOptions{margin_cells,growth_ratio,full_width_m,full_height_m}.
Both full extents zero selects nested growth; positive extents select fitted
axes. The kinetic grid/domain remains options.nx/ny/width/gap. Graded mode
requires finite foils, a 25 um-or-finer core, growth 1..1.15,<=300 margin cells
and<=1.01 M field nodes. The tested production-size case is 360 x 360+120-cell
margins,1.05 growth,361201 nodes,356.912 mm full width.

Transport insertion/escape, drift splitting, edge currents, continuity and
surface injection retain their original kinetic locator and function bodies.
Poisson deposition and field gathers use field_locate and actual physical
MFEM Jacobians. The original core coordinates and source weights are unchanged.
An explicit CoreDofMap connects old core DOFs to full field DOFs. Field-only
vacuum cannot receive a kinetic particle through the insertion API.

The actual heating driver is patched: snapshot scalar values use a typed
kinetic-cell sample, preserving the old exact(.5,.5) reference point; unit-bias
psi uses potential_value(x,y). Main phi and unit psi require identical field
topology. The native provider is unchanged: Fields::value receives PIC
coordinates and shifts internally; raw native-grid sampling uses the separate
PicNativeGridPoint conversion. No double shift or electric-domain translation.

This path is electrostatic plus prescribed native induction/B. It must never
be added to already-total Maxwell E/B. A Maxwell owner must bypass native
Poisson force E or explicitly own a compatible correction. This candidate
does not implement or qualify that Maxwell mode. Field extension also does
not change the kinetic escape box,25 mm axial closure, sheath model or source
geometry, and supplies no 3D/infinite-domain convergence claim.

Topology/checkpoints

The SHA256 identity includes actual axes, physical mesh vertices/connectivity,
Q1 element DOFs, essential owners, core DOF map, epsilon0, foil rectangles,
kinetic geometry/depth and finite-z bounds. Main/unit voltages are excluded
from geometry identity but still checked in checkpoint options as before.
Unversioned mesh mutation is detected before checkpoint save/load.

PIC3 binary order: uint64 magic 0x5241494b50494333; existing 168-byte options;
64 ASCII SHA256 bytes; existing 200-byte diagnostics; uint64 next_id; uint64
particle_count; existing 80-byte particles. Caller must construct the intended
field configuration before load; a changed digest rejects before import.
Default load rejects PIC2. Explicit import_legacy_kinetic_checkpoint(path,sha)
requires the exact legacy file SHA256 and records a field-model migration.
It never bypasses a mismatched PIC3 identity. A larger field does not silently
inherit the interpretation of an old kinetic checkpoint.

Reaction/arena integration boundary

The milestone reaction source is unchanged. Its 114 default-uniform checks
pass. Its existing load_bundle constructs GapPic(options), so graded bundle
restore fails closed until the arena owner threads PicFieldMeshOptions and
topology identity/configuration through restore and fingerprints. Same-owner
staging already keeps its existing immutable field mesh. No graded arena
restore or total-Maxwell qualification is claimed here.

CPU validation

The new field suite passes 30 checks: exact nine-array original 9 mm CSR/RHS
control; bitwise graded axes/source against the accepted nested 120 layout;
352652 captured phi/psi particle/snapshot samples bitwise equal; source charge
and source-work identities; manufactured scalar/gradient solve; unchanged
surface/escape/finite-z/edge-continuity paths; PIC3 corruption/layout/mutation
rejection and explicit legacy import. Actual captured fixture solves/advances
are zero. Manufactured CPU drift/solve tests are bounded and separate.
Peak field-test RSS is about 321 MiB. Original foundation and Q1 suites and
all 114 unchanged reaction tests pass. No numerical tolerance changed.

The unchanged archived reaction code needs explicit-std=c++17 with this frozen
compiler/vendor JSON combination. The first generic-language build failure
is preserved; CMake now applies the same C++17 setting as its supplied build
recipe, without changing reaction code. Native field/foundation code uses the
existing frozen MFEM 4.10/GNU 16.2.1 dependencies; OpenSSL comes from frozen root.

Build/qualification

Use CPU 7,nice 15,one thread and the exclusive testing-core 7.lock. The complete
commands are in the output qualification-commands-v3.json. Set the recorded
PIC_REACTION_NATIVE_SOURCE_SHA256 when configuring CMake, then build with
--parallel 1. Build outputs: heating-history-field-mesh, test_native_field_mesh,
test_native_foundation,test_q 1_scalar_gather,test_native_reactions.

Driver CLI retains its five inputs and optionally adds FIELD_MARGIN_CELLS
FIELD_GROWTH, then optional EXPLICIT_LEGACY_CHECKPOINT_SHA256 for a deliberate
PIC2 state import. The compiled driver uses this milestone's native CPU/MPI/
HYPRE backend. No CUDA backend/driver execution was performed. GPU backend
merging, per-step traces and saved-state replay remain explicit parent-owned
qualification steps. Do not deploy this isolated candidate to production.
