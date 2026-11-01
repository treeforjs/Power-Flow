# Conditional 3-D foil mechanics and conservative thermal exchange

This extension replaces independent thickness rods with a transient 3-D Q1
solid per foil. It adds conservative reference-energy projection and a signed
thermoelastic exchange packet, plus a local-cell coupled temperature corrector.
The parent still owns native conductive heat integration and field solves.

## Geometry and support scenarios

The reference foil is 5 mm wide, 0.2 mm thick and 25 mm axial height; inner
surfaces are y=±2.5 mm. Feed is z=+12.5 mm, the nominal 5×5 mm short footprint
is at z=-12.5 mm, and the observation plane is z=0. The electrical short record
does not specify its structural thickness, joint stiffness or mechanical anchor.
No bridge volume or fixture stiffness is fabricated.

Three explicit scenarios use this geometry: all faces traction-free; the entire
short-end foil plane fixed in xyz; or both axial end planes fixed in xyz. They
are a support bracket for comparison, not a proof that displacement lies between
their results or a model of measured supports. Each foil is a disconnected body.
For the paired-source benchmark, lower motion is the explicit y-reflection of
the upper solution under symmetric conditions; lower native bulk is still absent.

The solver uses 2×2×2 Gauss integration, reference lumped mass and average-
acceleration Newmark. Thermal eigenstrain and E(T) retain the v1 data bounds.
The new strip-producer datum ν=0.29 is held constant as an explicit nominal
scenario, not a measured elevated-temperature function. Magnetic loading,
plasticity, buckling/fracture and phase motion remain separate. A nominal
von-Mises/proof-stress screen rejects conditional work transfer above the
declared elastic screen; passing the screen is not foil-batch certification.

Full-integration Q1 solids can be too stiff in thin-element bending. Patch,
wave and time-integrator tests establish numerical identities, not application
bending convergence. Source trials compare mesh and time steps explicitly.

## Geometry and source conservation

`deform_hex_mesh` calculates 3-D cell Jacobians, Q1 volumes, rho_ref*Vref/Vcurrent,
boundary triangle normals/areas and Q1 surface-area vectors. Gauss, corner and
center Jacobian checks reject inversion; these samples are not a global
self-contact proof. Boundary vector-area closure and divergence-theorem volume
are checked independently. Material face IDs identify axial endcaps separately
from perimeter source patches.

`midplane_gap` intersects the moved inner surfaces with physical z=0, without
extrapolating beyond either face. The reported fixed-voltage `E=V/g` change is
a parallel-plate shape sensitivity only. It contains no edge, sheath, inductive
or current-feedback solution and does not authorize a moving-field handoff.

Source projection evaluates the separately sealed nominal caloric model on
the matched native temperatures, remaps extensive reference cell energy, then
inverts energy to T. Full 25 mm depth enters once. This conserves the selected
model's mapped energy; it is not the original constant-property heating run's
energy ledger. Native T/Az/rate remain hash-bound to the same 256×32 run.

## Work and thermal-stage interface

Engineering Voigt strains use [xx,yy,zz,2xy,2yz,2xz] and conjugate physical
stresses. The mechanical solver returns the path integral over each accepted
subcycle, `Wstress=∫sigma_average:Delta_epsilon dVref`. For free or immobile
homogeneous supports, `DeltaK+Wstress=0` to the discrete solve tolerance.

Recoverable elastic energy Phi is not the elastic internal energy. Under the
declared Helmholtz model the latter is `Uel=Phi-T*dPhi/dT`, which may be negative
as an excess energy. The signed debit from base caloric energy is
`DeltaUel-Wstress`, while the distinct algorithmic strain-energy reservoir debit
is `DeltaPhi-Wstress`. The code exports both and never substitutes one for the
other. The detailed derivation is in `solid-thermoelastic-work-derivation.md`.

`thermoelastic_exchange` checks mesh/quadrature/Voigt/material/reference identity,
recomputes constitutive stresses and local path work, checks kinetic balance,
and exports per-cell signed joules with state/mesh/material/source/transaction
identity and interval-end thermal stage. A positive debit produces a negative
thermal load. No species inventory is changed. This is an **internal reversible
exchange**, never authored external heat, solid Joule heat or PIC electric work.
Native integration needs a dedicated signed exchange ledger entry and matching
cell projection; it must not hide the packet in an external heat-source budget.

`coupled_cell_heat_step` repeatedly solves end-stage mechanics and temperature
until every cell and the combined internal-plus-kinetic energy close against
supplied cell heat. Its diagonal V*Ceps update is a quasi-Newton preconditioner,
not the full coupled Jacobian. Only the actual residual accepts a trial.
The function mutates no accepted state; the parent must preflight and commit
thermal, geometry, emission and field states atomically. This executable local
caloric test omits conduction and is not the existing native MFEM heat driver.

The caloric C ABI is loaded read-only after its source seal, header, implementation,
scenario and compiled library hashes match. It uses the owner's explicitly
conditional Cp-as-Cv reference-volume law with mechanics' reference density.
No alternate density default, new property fit or high-temperature extension is
introduced. The interface is a conditional constitutive model, not a physical EOS.

## Differentiability and admissibility

The v1 cumulative expansion and E tables are continuous but piecewise linear.
Their derivatives jump at coefficient knots, so `Phi-T*Phi_T` can jump for a
strained state. A physical-form thermal transaction must not treat that fitting
artifact as latent heat. The new port rejects any path crossing or landing on
a coefficient knot. The source demonstration therefore uses the first smooth
segment, with all observed source faces below a 368 K screening cap. The broader
773.15 K mechanics-data bound is unchanged and is not a melting temperature.
The separate nominal melting sources disagree and remain explicitly recorded;
neither is treated as a calibrated foil phase diagram.

## Running the bounded verification

```sh
bash tools/test_ss304_thermal_expansion.sh -m pytest -q \
  tests/test_solid_foil_mechanics.py tests/test_solid_foil_geometry.py \
  tests/test_thermoelastic_energy.py tests/test_solid_foil_coupling.py \
  tests/test_solid_foil_coupling_review.py
bash tools/test_ss304_thermal_expansion.sh tools/run_solid_foil_supports.py \
  --output /absolute/path/to/new-report-directory
```

Use the shared CPU7 lock, frozen userspace, nice15 and one thread. The report
includes source/property bindings, actual reference-mass/energy mapping errors,
support/gap/conditional-field comparisons, mesh/time comparisons, signed work
packets, corrected temperatures, white plots and saved displacement geometry.
No production state, old-run process, dependency or shared checkout is modified.
