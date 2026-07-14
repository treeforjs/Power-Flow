# MHD Course Notes For Solver Direction

Source reviewed: `C:\Users\tsmit\Documents\Papers\porto_course.pdf`

The PDF text layer is badly encoded, so this note records the reliable course
structure visible from the extracted table of contents rather than detailed
equation transcription.

## Useful Numerical Themes

- Put MHD in conservative/hyperbolic form wherever possible. This matters for
  shocks, steep gradients, and energy accounting.
- Prefer finite-volume or finite-element conservative discretizations for the
  field and material equations. Finite differences are useful for prototypes,
  but they should not become the trusted foil model.
- Track consistency, stability, and convergence explicitly with verification
  tests. Visual plausibility is not enough.
- Treat weak solutions and conservative discretization carefully; shocks and
  discontinuities are expected in pulsed-power material response.
- Use high-resolution shock-capturing ideas for compressible material/plasma
  evolution: Riemann solvers, TVD/MUSCL-style reconstruction, slope limiters,
  entropy fixes, and positivity fixes.
- Keep the magnetic field solenoidal. The course highlights the standard
  choices: nonconservative treatments, constrained transport, and projection
  schemes. For this project, constrained-transport-like updates should be
  preferred when the EM solver is upgraded.
- Treat time integration as a core design choice, not an afterthought. Explicit
  stepping, operator splitting, fully implicit/semi-implicit schemes, IMEX
  schemes, Jacobian evaluation, linear solvers, and approximate implicit methods
  are all relevant.
- Code organization should separate equation modules, boundary conditions,
  spatial discretization, temporal discretization, and software design. This
  supports the AGATE-like framework shape while keeping the physics ALEGRA-like.

## Implications For This Repository

- The current finite-difference reduced MHD step should be treated as a scaffold.
  It is not the final numerical method for a trusted foil calculation.
- The next EM/MHD solver should be written as a module with an explicit equation
  set and conservation properties, not as another patch inside
  `ReducedMHDSolver`.
- Verification problems should be added before relying on Mykonos-like geometry:
  magnetic diffusion in a slab, current penetration/skin depth, Orszag-Tang or
  Brio-Wu-style MHD tests for later full-MHD modes, and thermal-conduction
  manufactured solutions.
- The result visualization should include conservation/error diagnostics:
  integrated current, magnetic energy, Joule energy, thermal energy, neutral
  inventory, and `div B` or the appropriate 2D equivalent when available.
- Mesh refinement should eventually be nonuniform or block-structured. Uniform
  20 um cells are acceptable for a first refined foil run but will not scale once
  velocity-space neutrals and species chemistry are fully active.
