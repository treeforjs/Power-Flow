# Framework Roadmap

This project should follow an AGATE-like framework shape while staying close to
ALEGRA-MHD for physics choices, coupling strategy, diagnostics, and model
hierarchy. In short: AGATE is the software ergonomics reference; ALEGRA is the
physics reference.

## Closest Framework Reference

Primary framework reference:

- C. Bard and J. Dorelli, *High-Performance Computational
  Magnetohydrodynamics with Python*, Computer Physics Communications, 2026.
  DOI/article page: `https://www.sciencedirect.com/science/article/pii/S0010465526000597`
  Preprint: `https://arxiv.org/abs/2503.20899`

Useful design lessons for this repository:

- Keep the user-facing layer Python-native.
- Separate high-level interfaces from numerical implementations.
- Make grid, state, equations, boundary conditions, solvers, diagnostics, and
  backends swappable components.
- Support NumPy first, then CuPy/CUDA and later compiled C++ kernels for heavy
  loops.
- Prefer standard benchmark problems and regression tests over visual-only
  validation.
- Move large field outputs toward HDF5/Zarr-style chunked diagnostics rather
  than keeping all sampled fields in memory.
- Do not import AGATE's heliophysics assumptions wholesale. Its ideal/Hall/CGL
  examples are useful as framework examples, but this project should remain
  focused on resistive pulsed-power foil/electrode physics.

## Sandia Physics Reference

Primary physics reference:

- SAND2003-4074, *ALEGRA-MHD: Version 4.0*, T. A. Haill,
  C. J. Garasi, and A. C. Robinson.
- Sandia ALEGRA Multiphysics capabilities page:
  `https://www.sandia.gov/alegra/capabilities/`
- OSTI 2585866, *Progress on Extended MHD Modeling in ALEGRA*,
  M. M. Crockatt and A. C. Robinson, 2024 ZFS Workshop.
  `https://www.osti.gov/servlets/purl/2585866`

Useful physics lessons for this repository:

- Do not prescribe the spatial current-density profile.
- Derive `J` from the magnetic/electromagnetic solve and impose only global
  drive constraints such as total current or a circuit model.
- Use implicit electromagnetic solves for stiffness from high conductivity.
- Treat `E`, `B`, and `J` consistently enough that Joule heating and Lorentz
  forces are compatible with the discrete field solve.
- Move toward coupled or IMEX time integration once hydro/material motion is
  more than a surrogate expansion model.
- Preserve ALEGRA-style naming in internal design notes where practical:
  transient magnetics, thermal conduction, conductivity models, Joule heating,
  Lorentz force, material/void limits, circuit coupling, and diagnostics.
- Keep the capability ladder ALEGRA-like: low magnetic Reynolds number,
  resistive MHD, generalized MHD, then HEDP extensions such as radiation,
  opacity, and two-temperature physics.

## Conductivity And Resistivity References

Useful closure reference:

- P. J. Dellar, *Lattice Boltzmann magnetohydrodynamics with
  current-dependent resistivity*, Journal of Computational Physics 237,
  115-131, 2013. DOI: `10.1016/j.jcp.2012.11.021`

This should inform optional anomalous/current-dependent resistivity closures,
not force the whole project into a lattice-Boltzmann formulation.

## Numerical Methods Reference

Local course overview:

- `C:\Users\tsmit\Documents\Papers\porto_course.pdf`

Useful numerical-method lessons:

- Use conservative finite-volume/finite-element formulations for trusted MHD
  work.
- Add high-resolution shock-capturing methods before modeling compressible
  plasma/material shocks.
- Keep magnetic-field solenoidality explicit in the method design, preferably
  with constrained-transport-like updates where appropriate.
- Treat time integration, splitting, implicit solves, Jacobians, and linear
  solvers as first-class framework components.

## Proposed Package Shape

Near-term target structure:

- `mhdlab.grid`: geometry, raster/nonuniform mesh, metric data
- `mhdlab.state`: material, EM, neutral, chemistry, and diagnostic state objects
- `mhdlab.equations`: reduced RMHD, magnetic diffusion, thermal conduction,
  kinetic neutrals, CR kinetics
- `mhdlab.solvers`: time integrators, linear/nonlinear solvers, drive/circuit
  coupling
- `mhdlab.backends`: NumPy, CuPy, and C++/CUDA implementations behind one API
- `mhdlab.io`: config loading, HDF5/Zarr diagnostics, run manifests
- `mhdlab.viz`: quick-look plots, animations, spectral comparisons

## Immediate Technical Direction

The current code remains a reduced prototype. The next significant refactor
should make the code more ALEGRA-like in physics structure rather than adding
more ad hoc physics:

1. Split `ReducedMHDSolver` into EM, thermal, and material-motion components.
2. Replace direct field-history accumulation with chunked diagnostic output.
3. Add standard MHD verification problems before trusting the foil case.
4. Add nonuniform or surface-refined grids so a 200 um foil can be resolved
   without exploding neutral velocity-space memory.
5. Move from scalar thermal-expansion displacement to a real momentum/material
   model once the EM and thermal pieces are stable.
