# ALEGRA-MHD Notes For This Prototype

Source reviewed: `C:\Users\tsmit\Documents\Papers\918373.pdf`

Report: SAND2003-4074, *ALEGRA-MHD: Version 4.0*, Thomas A. Haill,
Christopher J. Garasi, and Allen C. Robinson, Sandia National Laboratories,
November 2003.

Additional public capability reference:

- Sandia ALEGRA Multiphysics capabilities page:
  `https://www.sandia.gov/alegra/capabilities/`
- The page distinguishes ALEGRA-EM and ALEGRA-HEDP. For this project,
  ALEGRA-EM is the nearer baseline: electromechanics/MHD, continuum mechanics
  coupled to electromagnetic induction and Lorentz force, conductivity models,
  and coupled circuit behavior. ALEGRA-HEDP adds radiation transport, opacity,
  and two-temperature plasma options that are longer-term extensions.

## Relevant Model Structure

- ALEGRA-MHD couples hydrodynamics/solid mechanics, transient magnetics,
  thermal conduction, conductivity models, and optional radiation emission.
- Sandia's public capability list reinforces the same hierarchy: ALE mechanics,
  large deformation and shock handling, low magnetic Reynolds number MHD,
  resistive MHD, generalized MHD, Joule heating, thermal transport, material
  response models, and lumped circuit coupling.
- The 2D Cartesian mode supports in-plane magnetic field components with
  out-of-plane `Jz`, matching the foil/load geometry used here.
- The report treats current density as derived from the magnetic field or vector
  potential, not as a prescribed spatial profile. For this prototype, `Jz` should
  come from a magnetic diffusion/induction solve with the total circuit current as
  a constraint.
- Magnetic diffusion time control uses `D = 1 / (sigma * mu0)` and a cell-size
  criterion proportional to `h^2 / (4D)`. This matters strongly for refined foil
  meshes and high-conductivity material.
- Joule heating should be consistent with the discrete magnetics solve. ALEGRA
  also discusses limiting mixed material/void cell overheating; the current
  prototype does not yet have a multi-material cell model.
- Thermal conduction is modeled as `rho * Cv * dT/dt = div(k grad T)`. The
  prototype has only an explicit finite-difference version and should move toward
  an implicit solve as mesh refinement increases.
- Conductivity is material-state dependent. The report lists LMD, Spitzer,
  anomalous, Knoepfel, and tabular Sesame-style options. The current prototype
  still starts with a simple SS304 conductivity and needs state-dependent
  electrical/thermal conductivity models.

## Implications For The Power Flow Code

Immediate corrections:

- Do not use a user-forced edge current profile.
- Preserve the imposed total current, but derive local `Jz` from a reduced
  induction/magnetic-diffusion solve.
- Refine the foil mesh enough to resolve the 200 um foil thickness. The main
  Mykonos config now targets roughly 20 um cells, giving about ten cells through
  the foil.
- Keep a coarse preview config for fast debugging.
- Store diagnostic samples as `float32` by default to reduce memory pressure at
  1 ns cadence.

Next physics upgrades:

- Replace the current scalar pressure/expansion proxy with momentum and energy
  equations or a clearly labeled solid/thermal-expansion surrogate.
- Add state-dependent conductivity, beginning with a low-temperature
  metal/Knoepfel-like model and later moving to tabular or LMD/Sesame-style data.
- Add magnetic diffusion and Alfven time-step diagnostics to the run summary.
- Move thermal conduction to an implicit solve before using very fine uniform
  grids.
- Add nonuniform/surface-refined meshing so the foil skin/current penetration
  region can be refined without exploding the neutral velocity-space memory.
- Preserve an ALEGRA-like model ladder in the code: low magnetic Reynolds number
  mode, resistive MHD mode, then generalized/Hall-like MHD mode only after the
  reduced resistive model is verified.
- Add explicit circuit/mesh coupling so the load has an effective resistance and
  inductance rather than only a prescribed current waveform.
