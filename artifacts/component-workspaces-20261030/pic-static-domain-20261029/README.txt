Static electrostatic-domain diagnostic

Status: prepared and CPU checked; no GPU solve or domain-sensitivity result.
This is an isolated field-only program, not a restart or a particle simulation.
No shared production source, checkpoint or field input is edited.

Bound input is the real46726-particle midpoint fixture at26.478443519271536ns.
Native types and HostSnapshot::load validate its raw arrays. IDs, particle order,
positions, proper velocities, q, m and weights remain byte-for-byte unchanged.
options.bin is read as a geometry record, never bypassed as a restart. Width9mm
uses the exact captured floating-point geometry value, not a recomputed9*.001.

Uniform25um cases:
  9mm:360x360 cells,130321 nodes
 13mm:520x520 cells,271441 nodes
 17mm:680x680 cells,463761 nodes
 25mm:1000x1000 cells,1002001 nodes
The hard ceiling is1010000 nodes. The graded exterior proposal is deferred;
the production uniform locator is never used on a nonuniform mesh.

Coordinate contract:
  X=x_captured-4.5mm, Y=y_captured-4.5mm.
  shift=(new_width-captured_width)/2.
  x_static=x_captured+shift, y_static=y_captured+shift; z is unchanged.
Both foil rectangles shift by the same amount in the new positive-coordinate
mesh. Physical foil width5mm, gap5mm, thickness200um and extrusion height25mm
remain unchanged. Native gathers ALWAYS use the original captured x/y minus
4.5mm. Static coordinates are never passed to the hard-coded native translation.

Equations and boundary contract:
MFEM H1 Q1 assembly uses exactly epsilon0=8.8541878188e-12 and point-charge
deposition q*w/depth as in the sealed solve_fields. Internal foil nodes are
Dirichlet0/V; all outer edges have homogeneous natural Neumann data. The initial
outer-y trial projection is preserved but not retained as an essential set.
For every domain the real unit Laplace psi and full Poisson phi are separate
RHSs of the SAME constrained matrix. The cuDSS-enabled program factors once,
solves twice, uses the unchanged original CudssPoisson header, FP64, exactly two
iterative refinements and the original residual guard:
  ||A x-b||2 <= max(1e-30,2e-13 ||b||2).
No tolerance or solver configuration is exposed as a command-line option.

Only static charge deposition is performed; no emission-source law is evaluated.
The field gather composes the actual phi/psi solution, actual saved native
A_z/rate pair and unchanged fraction/voltage. It does not invoke Boris or any
particle/source/heating/collision routine. FieldRecord output contains E[3],
B[3],phi,psi as8 doubles(64B), one record per original particle in original order.
Candidate phi/psi nodal arrays are raw FP64 vectors, not checkpoints.

CPU prerequisite:
static_domain_cpu has no CUDA linkage and no linear solve. It reconstructs the
original9mm topology/operator/charge RHS, loads the captured phi/psi and checks
their residuals and gathered E/B against the original force-loop references.
The topology/DOF map must match self_cells.bin bytes. An independent source route
uses the captured cell geometry, inverse local coordinates and Q1 basis. It is
not bitwise identical because the inverse-coordinate arithmetic differs from
the sealed normalized-coordinate locator; the measured discrepancy must lie
below1024 epsilon times total absolute deposited charge per unit depth. This
uses the existing ledger-roundoff factor, not a physical-source tolerance.
Raw source is copied BEFORE MFEM FormLinearSystem modifies its RHS. Both RHS
LinearForms remain alive while their MFEM b/x views are used.

The CPU prerequisite serializes full/constrained CSR and raw/transformed RHS
arrays for reviewable hashes. No CSR was stored by the original midpoint
capture; exact topology, literal sealed MFEM assembly and exact reproduction of
both saved residuals plus captured E/B provide the control evidence. Do not
describe this as a comparison against an independently stored captured CSR.

Acceptance gates:
  captured E:1e-6V/m +5e-12 abs(reference component)
  captured B:1e-13T +5e-12 abs(reference component)
  phi/psi nodal control:1e-10 global reference scale +256 epsilon scale,
    consistent with the existing checkpoint field regression policy.
  source sum:1e-11 times max(abs(direct charge),e).
  independent raw source:1024 epsilon abs(direct charge)/depth.
  original linear residual guard above for both RHSs.
For extended domains, reference field equality is NOT expected and is NOT a
gate; comparisons against9mm quantify sensitivity. No EEDF/trajectory result is
accepted or implied. Particle energies are fixed throughout this static study.

Parent execution boundary:
Only the parent may run static_domain_gpu with PIC_STATIC_DOMAIN_PARENT_GPU=1.
The parent request requires the9mm GPU solve to pass all captured-control gates
BEFORE any13/17/25mm case. Hash verification is required before/after every case.
Use shared CPU6/testing-core.lock and the separate testing-gpu.lock serially;
keep production frozen. No process-control logic is embedded in this package.
The parent must enforce the declared wall/RSS/GPU/output limits. This package
does not resume, freeze or signal production processes.

Resource contract proposed for parent review:
  one CPU rank, CPU6,nice15,one job; one geometry at a time.
  <=1010000 nodes;<=46726 particle field outputs;two RHS solves/one factor.
  per-case wall cap900s;host RSS cap4GiB;GPU used-memory increase cap8GiB.
  output cap64MiB/case;parent must stop the diagnostic case on a cap breach.
  resource limits are execution ceilings, not measured memory/time guarantees.
The largest case has~9.0million structural CSR entries. GPU factor fill/memory
is not measured yet; parent should inspect available memory before dispatch.

Build uses the existing MFEM4.10, GCC and pinned CUDA/cuDSS libraries only.
Compile with CPU6,nice15,the shared lease and --parallel1. input-bindings.json
binds source, binaries, headers, libraries and every actual fixture payload.
tools/verify_inputs.py checks these bindings without launching any computation.

See the parent request and final CPU receipt in:
<workspace>/outputs/pic-profile/em-field-consistency-audit/static-domain-preparation
