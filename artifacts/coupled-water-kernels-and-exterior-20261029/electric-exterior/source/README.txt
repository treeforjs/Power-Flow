Static graded electric exterior — CPU-qualified candidate, not integrated

Prefix: /home/treefor/projects/runtimes/pic-electric-exterior-20261029
Evidence: /home/treefor/sim-workspace/powerflow/outputs/pic-electric-exterior-20261029

The candidate appends graded vacuum quadrilaterals to the exact captured
360 x 360, nominal 25 um, 9 mm electric mesh. Particle and foil positions
never translate. The complete original fine region remains unchanged,
including the finite foil volumes and all 3618 Dirichlet nodes. The two
electrode potentials, point-charge source, 25 mm represented depth,
epsilon0, FP64 arithmetic, original residual tests and two cuDSS refinement
iterations remain unchanged. Native E/B still samples original particle
positions minus (4.5,4.5) mm; outer electric coordinates never enter that
mapping. No trajectory, source injection or current deposition is run.

The remote outer electric boundary remains the original natural homogeneous
Neumann condition. Moving this artificial boundary is a truncation study;
this is not a derived open-boundary condition. There is no Robin constant,
sheath smoothing, magnitude rescaling or thesis matching. Neither finite
extent nor a small successive difference proves the infinite-domain limit.

source/electric_coordinates.hpp provides non-interchangeable ParticlePoint,
ElectricPoint and NativePoint types, Layout/Axis, point deposition and scalar
value/gradient gathering. Core location preserves exactly the old normalized
x/width*nx arithmetic. Outside the core it searches physical nonuniform
axis intervals and uses the corresponding MFEM element transformation.
Each graded quadrilateral has its actual physical Jacobian. Shape weights
and gradients are not evaluated using a nominal 25 um exterior spacing.

The adapter is static-only. It is not a replacement for the evolving GapPic
driver: that driver still couples its Cartesian locator, current/continuity
grid, collision walls and absorbing boundaries to its original mesh. Those
paths require a separate coordinated integration and replay review. Enlarging
the electric field mesh does not enlarge the particle transport box, expand
the native B data, add water/plasma chemistry, or supply heat feedback.

CPU controls (cpu-controls.json, all passed)

* Uniform limit: full and constrained CSR indices/values and all three RHS
  arrays have exact matching SHA256 hashes against the prior CPU-qualified
  control. Captured phi/psi and all 46726 particle E/B values are unchanged.
* Every tested graded case: central cell coordinates and mapped charge RHS
  are bitwise equal to the capture. Transferring only the captured core nodes
  gives exactly zero E/B differences for every captured particle. Outer nodes
  in that test are zero test data, not an exterior solution; no residual or
  physical energy claim is made for those transferred vectors.
* Manufactured harmonic bilinear solution: max gradient error 9.50e-12 V/m,
  nodal error 2.51e-14 V. Signed point sources conserve charge and analytic
  source work, including points in graded exterior cells.
* Smooth nonzero-source Poisson test: exp(20x)cos(15y), rho=-175 epsilon0 phi,
  exact outer Dirichlet data. Under three successive graded refinements,
  max value errors are 1.201e-4,3.033e-5,7.620e-6 V and gradient errors
  0.3711,0.1866,0.09357 V/m. This verifies approximation behavior for that
  manufactured problem; it is not foil or physical plasma convergence.
* Zero-growth and exact uniform locator limits, boundary points, invalid
  margin rejection and nested shared-axis equality are checked explicitly.

The tested fixed-width cases are full widths 13/17/25/40/80/160 mm, using
32/40/48/64/80/96 margin cells on each side. At full width 80 mm, 120 and
160 margin cells separately refine the same outer extent. The largest tested
case has 463761 nodes and 348940 KiB peak RSS (about341 MiB).

For boundary-distance isolation, width argument0 selects a nested family:
first outer cell25 um, ratio1.05, margin count1..120. Counts40/60/80/96/110/120
have been tested; all overlapping physical axis coordinates are bitwise
identical. The exact extents are recorded from the actual generated axes,
not rounded labels. All cases stay below1.01M nodes. A 4 GiB address-space
limit bounded CPU tests. Actual GPU host/VRAM factorization memory is not
measured and must be bounded by the parent's existing resource monitor.

New charge/work audit and its preserved first failure

The initial extra Gauss test used1024 epsilon times the cancelling electrode
charge scale. It rejected the unchanged capture: 2.76708e-20 C defect versus
1.28096e-20 C. The original solver did not impose that diagnostic budget;
its Gauss guard requires finiteness. The failed uniform-control receipt is
preserved unchanged. It is never relabeled as a pass.

The replacement diagnostic is derived from the complete rounded full-operator
action and free residual, with finite-sum roundoff accounting; see
gauss-accounting.txt. The unchanged capture has a2.61781e-20 C full-operator
null-mode action sum and4.94577e-20 C/m free residual norm. The resulting
bound is5.12399e-19 C. The exact finite-vector bookkeeping identity closes
within1.11632e-22 C. This accounts for the original defect; it does not
establish charge-only accuracy under a driven-bias residual budget.

The point-source work identity is checked on its own actual source-work
scale:1.58328e-25 J error versus4.20615e-24 J budget. The driven field energy
identity is a separate check and is not evidence for resolving a tiny
space-charge field. Original source, charge, solver, nodal and strict E/B
comparison thresholds have not changed.

Prior strict-control history remains authoritative

The old fresh9 mm static GPU solve failed2733 strict component comparisons;
the old accepted additive CPU review separately accounted for the nodal
roundoff and retained those failures. This candidate does not turn that
failed recipe into a pass. Its fresh9 mm GPU mode retains the original strict
gate and may return3 with useful output. Further static cases require a new
explicit parent continuation referencing the preserved additive review.
The earlier 9→25 mm study was not converged at25 mm and remains unchanged.

Commands and boundaries

All commands use CPU7,nice15,one thread and testing-core7.lock. The CPU and
GPU targets compile using existing frozen GNU16.2.1/MFEM4.10 dependencies.
Frozen dynamic-library resolution succeeds. The GPU target has not run.

CPU binary: build/electric_exterior_cpu
GPU static binary: build/electric_exterior_gpu
Manufactured checks: build/exterior_manufactured

CLI: BINARY cpu-preflight|cpu-geometry|gpu FIXTURE NEW_OUTPUT WIDTH_MM MARGIN
WIDTH_MM=9 requires MARGIN=0; WIDTH_MM=0 is the fixed1.05 nested family.
The literal fixture and argv arrays are pinned in gpu-request.json. GPU mode
requires PIC_STATIC_DOMAIN_PARENT_GPU=1. It factors once, solves captured
phi and the existing unit-bias psi, gathers static fields, and exits. It does
not call GapPic::advance, Boris or any evolution loop. Outputs must be fresh.

GPU request status is prepared for review only. CPU controls permit review
of a bounded static GPU test; they do not authorize execution or production
integration. Root owns dispatch, resource monitoring and any later evolution.
