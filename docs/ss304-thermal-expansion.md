# SS304 solid expansion and geometry interface

This implementation supplies bounded solid eigenstrain, static and transient
elastic benchmarks, conservative geometry/mapping controls, and a local
transaction interface. It does **not** qualify deformation of the actual
connected feed/foil/short assembly or a coupled production field solve.

## Material and strain

`data/materials/ss304-v1.json` freezes producer tables for Type 304, UNS S30400,
EN 1.4301, including source URLs, page locators and downloaded-byte SHA256s.
304L is not substituted. Elastic/expansion evaluation rejects outside
293.15–773.15 K. Density is assigned to the explicitly disclosed 293.15 K
reference. Conductivity and heat capacity in this record are room-temperature
data only; they do not define a temperature-dependent heating law or caloric EOS.

The published expansion coefficients are **means from 20°C**. At each knot,
`s(T)=mean_alpha(T)*(T-293.15)`. Linear interpolation of cumulative `s` yields
`lambda_th=(1+s(T))/(1+s(T_ref))`; the engineering eigenstrain is `lambda_th-1`.
Equivalently, integrate current-length tangent `alpha=s'/(1+s)` and exponentiate.
The caller must supply the stress-free `T_ref`; the matched-run benchmark uses
300 K. A different reference also changes reference density consistently.

`SolidExpansionMaterial` validates finite arrays and increasing tables.
`load_ss304_material` rejects a different grade/schema and records the material
hash. No endpoint holding, clipping or solid elasticity through melting occurs.
The separate producer TP304 melting interval 1672.15–1727.15 K does not extend
the elastic-data domain. Vaporization and latent heat are not implemented.

## Mechanics implemented

`static_bar` solves the exact equilibrium of piecewise-constant uniaxial elements:
`N=-sum(eps_th*h)/sum(h/(E*A))` for fixed ends; `N=0` for a free end.
The free case removes a rigid translation by setting the mass-center displacement
to zero. `fixed_left` and `fixed_both` are separate explicit benchmark conditions.
For uniform heating, free strain is thermal strain and stress is zero;
uniaxial full restraint gives `sigma=-E*eps_th`.

`transient_bar` uses linear elements, lumped **reference** mass and Newmark
average acceleration (beta=1/4, gamma=1/2) with the endpoint temperature-dependent
modulus and thermal nodal force. The free rigid mode is retained. Initial
displacement and velocity default to zero, not a hot equilibrium. No damping,
Lorentz force or magnetic pressure is added. The 2% total-strain guard is a
numerical small-strain restriction, not a yield criterion.

`isotropic_elastic_stress` provides a separate 3-D Hooke-law trial tensor, with
caller-declared Poisson ratio. Fully restrained isotropic stress is
`-E*eps_th/(1-2*nu)` on each diagonal. No grade-qualified `nu(T)` is silently
invented. Rods have a constant reference cross-section and omit Poisson response,
bending, buckling, fracture, contact, plasticity and fixture loads.

At 300 K the extensional crossing times are roughly 40 ns through thickness,
1 µs across width and 5 µs axially. Static expansion is an equilibrium benchmark;
it cannot be interpreted as actual 60 ns foil motion. Even a temperature-admitted
elastic trial can exceed the annealed-bar proof-stress proxy. That proxy is a
warning, not a certified foil yield curve or plastic return map.

## Source admission

`load_matched_history` pins the new 256×32 run's `source_run.json` SHA256,
verifies all consumed files, and checks exact matched temperature/Az/rate clocks.
Upper bulk temperatures already include 300 K. Eight directly exported face
tables are increments; the loader adds 300 K exactly once. Lower bulk was not
exported. It remains absent by default; explicit mirror reconstruction checks
all lower face histories and remains labeled unverified bulk.

Admission checks every available cell and surface patch **before** averaging.
The first invalid observation ends the entire admitted prefix, even if a later
frame cools. Recorded surfaces exceed 773.15 K at 17.622896529 ns; the last common
safe bulk/native-field frame is 17.413910165 ns. Cell averages and exported faces
still cannot bound unobserved nodal/subcell maxima.

The matched native solve contains no mechanical supports. Its electrical feed,
short and thermal insulation do not specify displacement or traction.
The intended target remains width 5 mm, thickness 0.2 mm, gap 5 mm and height
25 mm, with feed z=+12.5 mm, short z=-12.5 mm and observation z=0.

## Geometry, mapping and energy ownership

`affine_foil_geometry` returns reference coordinates, transformed vertices,
six face areas/normals (Nanson formula), positive determinant, physical volume
and `rho=rho_ref/det(F)`. Mass is invariant; reference density must not receive
a second thermal correction. The two-foil projected-y gap rejects overlap.
Own-center isotropic free expansion narrows the gap by `thickness*(lambda-1)`.
An alternative translation changes that result and must be declared.

`conservative_remap_1d` integrates piecewise cell overlaps on matching complete
material-coordinate domains. It checks coverage, positivity and conservation
at fixed relative tolerance 5e-13. Use extensive mass/energy inventories or
their proper densities. A remap of prescribed temperature is only a sampling
choice; conserving its integral is **not** conserving thermal energy.

`transform_field` distinguishes material vectors (`F v`), covectors (`F^-T E`)
and flux/Piola vectors (`F B/J`). Their geometric integral tests do not replace
a compatible Maxwell solve on the new mesh. `qualify_field_handoff` reconstructs
geometry and transfer evidence, admits only affine benchmarks and leaves
`production_field_solve_authorized=false`. Nonuniform transient rod geometry
cannot pass this field gate.

Six-face area for the two foils is 5.24e-4 m². The familiar 5.20e-4 m²
perimeter-times-height source area excludes axial endcaps. Patch masks must
retain this distinction; depth and volume must each enter only once.

`GeometryPort` binds material/source hashes, state/mesh revisions, accepted
`[t,t+dt]`, stage time and transaction identity. Prepare does not mutate accepted
geometry; commit checks revisions; rollback discards prepared proposals.
Duplicate, stale and forged proposals fail. Empty species event/parent lists
state that this port mutates no species. The caller owns global atomic acceptance.
The heat-feedback v1 stage contract is enforced: the trial geometry consumes
interval-end thermal state (`stage_time_s == t+dt`); emission for that interval
uses the prior accepted interval-start temperature. Middle-stage and start-stage
geometry proposals cannot be committed as an end-stage state.
Current replay protection is in-process; checkpoint persistence and rollback of
accepted commits are not implemented.

This port explicitly rejects physical coupling. Kinetic and elastic trial
energies are diagnostic only. Thermal/mechanical work and heat debit are `null`,
not fabricated zeros: a consistent thermoelastic caloric law, work exchange and
parent atomic coordinator must be qualified before activating deformation in a
physical heat/EM pilot. A prescribed thermal history is not feedback.

The separate `ss304-thermal-properties-20261029` work owns extended-temperature
enthalpy/transport/phase laws. Its temperature range does not automatically
extend this component's selected **mechanics-data** bound of 773.15 K. This
bound is not a melting temperature. See `interfaces/ss304-geometry-port-v1.json`
for the concrete common-port binding and unresolved work port.

## Reproduction and integration

From this checkout, use the frozen dependency stack and shared CPU7 lock:

```sh
bash tools/test_ss304_thermal_expansion.sh -m pytest -q \
  tests/test_thermal_expansion.py tests/test_thermal_expansion_source.py \
  tests/test_thermal_expansion_mapping.py tests/test_thermal_expansion_port.py
bash tools/test_ss304_thermal_expansion.sh tools/run_ss304_thermal_expansion.py \
  --output /absolute/path/to/new-report-directory
```

The launcher uses CPU7, nice15, one thread, no GPU and the frozen loader.
A held lock fails quickly; queue and retry without moving to another CPU.
No native rebuild or dependency update is required. The CLI refuses to overwrite
an output directory. It saves input/output manifests, white plots, admission
decisions, and source-driven free thickness-bar **trial** results at the edge
and midwidth. Mesh/time refinements hold the native source fixed and do not
establish source-grid convergence or physical convergence.

The tests cover zero change, integrated mean-alpha reference composition,
uniform free and constrained analytical stresses, nonuniform compliance,
second-order spatial/time wave limits, rigid motion, mass and unforced energy,
conservative mesh overlap, geometric line/flux invariants, source-hash tampering,
temperature/phase rejection, and transaction rollback/replay. Parent integration
should consume the readiness flags and keep physical deformation disabled until
the enumerated support, constitutive, work and mesh qualifications are supplied.
