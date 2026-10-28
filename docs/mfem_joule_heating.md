# MFEM Joule heating backend

MFEM's `miniapps/electromagnetics/joule.cpp` is the selected electrode-heating
backend. It replaces FLASH in the active workflow. The local MFEM 4.10 source
and MPI installation are already available in Ubuntu 24.04 WSL2. The existing
FLASH interfaces and prior verification records remain as historical/optional
inputs; acquiring FLASH is no longer a project prerequisite.

The target chain is MFEM electrode heating → conservative H2O desorption →
Raikou gap transport and resolved populations → observation → neutral
inventory/temperature inference. The gas geometry still determines the optical
path. Thermal penetration depth and electrode thickness are not optical lengths.

## Scope and inspected implementation

The stock Joule miniapp solves transient low-frequency electromagnetics and heat
conduction. It uses H1 potential, H(curl) electric field, H(div) magnetic/thermal
flux, L2 temperature and implicit time stepping. The coupling is one way:
electromagnetic fields heat the material, while temperature does not update the
electromagnetic coefficients. It has no fluid motion, moving surface, Hall
induction or self-consistent plasma/CR model. See the
[official description](https://mfem.org/electromagnetics/),
[v4.10 driver](https://github.com/mfem/mfem/blob/v4.10/miniapps/electromagnetics/joule.cpp)
and [solver](https://github.com/mfem/mfem/blob/v4.10/miniapps/electromagnetics/joule_solver.cpp).

This is a suitable starting point for the stationary solid electrode stage.
Qualification must establish its time/regime of validity. Melting, motion,
ablation, strong temperature-dependent properties and Hall physics require
separately implemented and verified extensions or external owner models.
Magnetized charged-particle transport remains with the gap-model owner.

Local source is `/home/tsmit/src/mfem-4.10`; the installed CMake package is
`/usr/local/lib/cmake/mfem`. The source is an extracted archive, so it has no
Git identity. Its three Joule source files were compared byte-for-byte with
official v4.10 files. Record file hashes separately from the installed library
identity; matching source files alone do not attest the library build.

## Stock example boundaries

The stock rod example is a demonstration, not a foil experiment:

- Its default mu=1, sigma=20*pi, capacity=1 and conductivity=.01 are not
  qualified SI SS304 properties. Temperature starts at zero. No ambient Kelvin
  offset is inferred from the name of the temperature field.
- Its exterior has electrical conductivity 1e-6 times the metal and thermal
  conductivity 1e6 times the metal. This is not physical vacuum.
- The terminal potential is +/-cos(omega*t), selected by the z coordinate.
  The stock example does not ingest the recovered current trace or close the
  experiment's circuit.
- Rod ends impose zero normal thermal flux; the outer thermal boundary has
  natural zero-temperature treatment. Define actual foil, contact and gap
  thermal boundaries explicitly.
- The driver assumes three-component fields and uses z in its drive. Use the
  supplied 3D mesh for smoke verification; do not infer unmodified 2D support.
- Use the long `--problem rod` option: v4.10 registers the short `-p` for both
  problem and visualization port. Use VisIt field output for native inspection;
  the stock individual `F` file path writes to a closed/wrong stream.

No demonstration result is an absolute material temperature or calibrated source.
The isolated consumer in `native/mfem_joule` builds the unchanged upstream driver
and solver against the installed library, without modifying its shared tree.

## Physical foil adaptation

An owned experiment driver based on Joule must supply the following before
surface temperatures are used for physical desorption:

1. A finite-width, three-dimensional or explicitly justified extruded geometry
   with both conductors, true foil thickness, edge radius, contacts and the
   exterior magnetic domain. Resolve skin depth and edge current crowding.
2. Sourced cold SS304 electrical conductivity, volumetric heat capacity and
   thermal conductivity in SI with a justified temperature range. Update
   coefficients, assembled operators and preconditioners consistently if
   temperature dependence is introduced. Stock one-way coupling must remain
   explicit until that update is qualified.
3. Shot-specific, bounded drive and registered time origin. Current drive,
   terminal voltage and circuit closure are different boundary prescriptions.
   Verify the resulting integral of current over each conductor cross-section;
   do not substitute the shot 17133 benchmark drive for spectroscopy shots.
4. Independently closed magnetic, thermal and boundary energy balances. Check
   direct integrated sigma*|E|^2 against the projected heating field. Record
   actual linear-solver residuals and completion; process exit alone is
   insufficient scientific verification.
5. One-sided metal-temperature trace on each emitting face, evaluated with
   explicit face quadrature and physical area. L2 temperature is discontinuous;
   a cell-centre value, neighbouring vacuum value or global bulk average cannot
   silently replace that trace. Refine before nonlinear release is evaluated.

The [surface contract](mfem_surface_contract.md) admits explicitly evaluated
samples with units, temperature basis, face/patch support, sampling semantics
and source provenance. It does not guess SI units or parse arbitrary opaque
native fields. Positive-depth samples require a justified boundary-proxy
decision before desorption.

## Verification ladder

1. Build/run the unchanged 3D stock rod on MPI and read native saved fields;
   check finite data and nonzero heating, retaining its demonstration status.
2. Zero-drive invariant and homogeneous Ohmic conductor: check imposed/current
   response, direct heating and thermal energy increase.
3. Analytic magnetic diffusion and independent thermal conduction, using the
   existing [diffusion](flash_heating_verification.md) and
   [conduction](conductive_heating_verification.md) reference APIs and the
   [MFEM profile comparator](mfem_heating_verification.md). Match point,
   cell-average and face-trace quantities explicitly.
4. Energy balance and independent mesh/time/order refinement; check solver
   tolerances and exterior-domain sensitivity.
5. Finite-width physical foil and the two separate source/heating hypotheses,
   followed by projected-feature convergence below 1% and numerical uncertainty
   below 10% of any claimed spectral contrast.

Passing the stock smoke establishes native execution and data readback only.
No physical foil temperature, water-release amount or experimental neutral fit
is established by it. New evidence is recorded separately from the preserved
FLASH-era setup and comparison records.

## Native build and saved-field readback

The isolated consumer compiled successfully against the installed library and
ran the unchanged v4.10 rod example on two MPI ranks. Its native reader loaded
`T`, `w`, `E`, `B`, `Phi` and `F` from the initial and final VisIt collections.
All saved field coefficients and sampled scalar fields were finite. Final
native time was 0.2, with positive volume integrals of temperature and Joule
heating. A small negative temperature undershoot was retained, not clipped.
This is execution/readback evidence, not positivity, convergence or energy
qualification of a physical foil.

```powershell
python tools/probe_mfem_joule.py --wsl-distro Ubuntu-24.04 --source-root /home/tsmit/src/mfem-4.10 --build-dir /home/tsmit/.local/share/power-flow/mfem-joule-build --output-dir runs/new_mfem_joule_demo
python tools/export_mfem_surface.py --demo --output-dir runs/new_mfem_surface_demo
```

The first command saves exact commands, exits, native logs, snapshot/source/
executable hashes and inspection results. Its output directory must be new.
The second command exercises the independent SI table handoff with authored
synthetic temperatures; it does not convert the native rod into a physical foil.
The native reader currently reports volume diagnostics. A planar material-face
evaluator producing qualified one-sided DG trace samples remains required for
the physical case.
