Bounded native electron/ion screening verification — 2026-10-29

Accepted evidence: evidence-003/results.json, 88 checks, 13 records.
Runtime: /home/treefor/projects/runtimes/pic-plasma-screening-20261029
This is an isolated verification component. It does not validate H2O breakdown,
kinetic Debye equilibration, a foil sheath, or the earlier electron-only halo.

Native code and ownership
Seven native/interface files were copied unchanged from the existing Hall
prototype, including GapPic, Q1ScalarGather and hall::CoupledPic. This reuses
native signed-particle insertion, Q1 charge deposition/gather, MFEM Poisson,
Boris transport, charge collection, native residual limits and the checked
species/ID owner. It adds no numerical shielding multiplier and changes no
native source. The copied source hashes, compiler/header/static-library inputs
and executable are bound in the manifests. Only the isolated benchmark is new.
No production checkpoint is read or advanced. No GPU or shared source changes.

Controls actually executed
1. Zero-density capacitor: native vacuum field agrees with the exact uniform
   field to 7.30e-11 V/m maximum absolute error.
2. Co-located equal-weight electrons and positive ions: signed particle charge
   and deposited charge cancel exactly; the biased vacuum field is bitwise
   unchanged. Each species has 144 macros, 1,777,680 physical particles and
   charges -/+2.848157358729119e-13 C. A nonzero gross inventory cannot be inferred
   from the zero net charge; the species totals are reported separately.
3. Closed common-velocity neutral pairs: 32 real native advances at 10 fs,
   zero applied E/B. Each species count, charge and kinetic energy is unchanged;
   native field energy, global Gauss defect and nodal continuity defect remain
   zero. This is a transport/cancellation control, not a warm screening history.
4. Linear Debye slab: three dimensional parameter choices,
   (n_e [cm^-3], Te parameter [eV]) = (1e10,2), (1e12,2), (1e10,20).
   Each has nx=2 and ny=16/32/64/128, with lambda_D/dy=4/8/16/32.
   The solution is invariant in x and z. dx=lambda_D/2 is held fixed because
   there is no x variation to resolve; only the normal y direction is refined.

Static linear-response formulation and independent reference
Let Te denote k_B T_e/e in eV, n0 the physical density in m^-3, immobile singly
charged positive ions n_i=n0, and a fixed total electron count. The assumed
small-amplitude isothermal Maxwell–Boltzmann susceptibility is
  n_e(y) = n0 [1 + (phi(y)-mean(phi))/Te],
  rho(y) = e n_i - e n_e,
  lambda_De = sqrt(epsilon0 Te / (n0 e)),
  -epsilon0 phi'' = rho.
The chemical-potential/mean subtraction makes the linear response invariant to
a scalar potential shift and preserves the total electron count. Both species
counts remain fixed through the diagnostic iteration to <=1.23e-14 relative.
The test uses L=4 lambda_De, phi(0)=+V0, phi(L)=-V0, V0=0.01 Te volts. By
antisymmetry mean(phi)=0. The independent continuum solution is
  phi(y) = V0 sinh[(L/2-y)/lambda_De] / sinh[L/(2 lambda_De)],
  Ey(y) = V0/lambda_De cosh[(L/2-y)/lambda_De] / sinh[L/(2 lambda_De)].
Hence Ey(mid)/Ey(vacuum)=2/sinh(2)=0.5514411295435665.

At every iteration, electron and ion physical number weights are constructed at
2x2 Gauss points in each native Q1 cell. The current native potential supplies
the next density response. Actual positive and negative particles are inserted
through the existing atomic authored verification API, and actual GapPic
charge deposition/Poisson solves produce the potential. Relaxation 0.35 applies
only to the nonlinear-solver guess. It is not physical damping or an electric
field attenuation factor. Each iteration is a fresh field-only diagnostic state,
not a source event in a continuous accepted physical history. All 12 cases
converged in 53 iterations to a relative fixed-point defect below 2e-12.

The response particles have zero velocities; Maxwellian thermal susceptibility
is assumed, not sampled or evolved. Do not advance these field-only states as a
physical thermal distribution. Do not add this susceptibility on top of kinetic
PIC electrons: that would count the electron response twice. In a coupled PIC
application, all species move and their signed charge/current provide screening.

Observed numerical refinement
Across all three cases, potential L2 error decreases 0.00442585 -> 0.00110809 ->
0.000277127 -> 0.0000692881, with observed finest order 1.99987. Electric-field
L2 error decreases 0.0619577 -> 0.0310933 -> 0.0155610 -> 0.00778230, finest
order 0.999667. These are the expected Q1 potential/gradient orders for this
smooth linear slab. Finest midplane field ratio is 0.551506775; the remaining
relative discrepancy from the continuum ratio is about 1.19e-4. This establishes
this bounded operator/closure check, not convergence of the foil application.

An initial diagnostic attempt stopped at a dimensionful |Ex|<1e-8 V/m test.
At the higher-density finest case, |Ex|=1.719e-8 V/m while the physical field
scale V0/lambda_D is about 1902 V/m. The accepted criterion records both values
and requires |Ex|/(V0/lambda_D)<1e-9; measured maximum is 9.04e-12. The native
Poisson residual and solver limits were not changed. evidence-001 is explicitly
partial; evidence-002 passed before runtime-resource recording was added.

Applicability diagnostics
applicability-parameter-grid.csv varies density 1e3/1e6/1e9/1e12 cm^-3 and Te
parameter 1/2/10/1000 eV. It reports lambda_D/25um, omega_pe, 1/omega_pe,
2*pi/omega_pe, dt*omega_pe for 10 fs and the captured 4.7762916 fs step, and
physical particles per Debye sphere. These are hypothetical diagnostics, not
Te/density fits to the non-Maxwellian halo. A small plasma timestep does not
resolve an under-resolved Debye length or establish equilibration. Static
iteration has no physical timestep; dynamic screening time refinement remains
untested. Strong bias, anisotropic/multibeam EEDF, magnetization, finite ion
response, collisions, emission/sheaths and moving surfaces can invalidate this
scalar static closure. The actual 110 kV foil bias is far outside its 1% linear
perturbation amplitude. No claim is made that ionization will fix the halo.

Primary references, checked 2026-10-29:
A. Beresnyak, 2023 NRL Plasma Formulary, Fundamental Plasma Parameters, printed
p.28: https://www.nrl.navy.mil/Portals/38/PDF%20Files/NRL_Plasma_Formulary_2023.pdf
R. Fitzpatrick, University of Texas at Austin, Debye Shielding:
https://farside.ph.utexas.edu/teaching/plasma/lectures/node7.html
The latter uses a responding two-species example; this test instead holds ions
immobile and uses the electron-only length derived above.

Joint-driver ports and remaining gaps
See driver-interface-inventory.json for exact paths/methods. GapPic already
solves one signed charge RHS from all stored particles. electric_field returns
only the Poisson/Dirichlet field. advance adds external E and replaces B with
its external callback at the force midpoint. A full Maxwell adapter that already
owns total E must not be added to that self-field blindly. One admitted total
field/particle owner is required. Native 2D Poisson has no spatial z response.

Native horizontal/vertical transfer arrays are integrated edge charge in C;
dividing by dt gives edge current in A, not yet a compatible 3-vector current
density. Jz and a geometry/volume-aware compatible field transfer are not exposed
by this API. Species-resolved continuity/work and proper units/stage times must
be retained. Midpoint push observers audit particle electric work but do not
supply a complete bias/circuit/Poynting or electrode heat ledger.

The newer compatible::ports::FieldAdapter has prepare/prepared/commit/abort,
SI field/current/charge/energy samples, binding and common-clock identities.
Its current scope explicitly accepts only periodic 1D3V fixed geometry with its
own electron/proton inventory, and rejects electrode/open boundaries, reactions,
solid heating and circuit/driver/electrode-energy ports. It cannot be treated as
an admitted foil field driver or co-evolve a second copy of PIC electrons.

Native wall charge/energy collection is aggregate by electrode. Finite-foil
sink face/normal/patch metadata need an admitted heat/source coupling port.
add_authored_batch is an external authored reservoir; internal H2O reactions
must instead atomically debit parents and replace/add all products through the
reaction transaction owner. These tests do not add that mechanism. Multi-species
water masses/charge states, qualified high-energy cross sections and product
kinematics, actual evolving material/source law, finite-electrode EM boundaries,
common-clock rollback/restart and coupling/time refinement remain integration
work owned by the joint driver.

Execution and reproduction
Native accepted case execution: CPU7 only, nice15, OMP_NUM_THREADS=1, frozen
userspace, 0.252 s native wall time. Resource details are in evidence-003.
Build/run with a NEW output directory (no production state is touched):
  /home/treefor/projects/runtimes/pic-plasma-screening-20261029/build-and-run.sh \
    /home/treefor/sim-workspace/powerflow/outputs/pic-plasma-screening-20261029/evidence-review
The wrapper holds testing-core7.lock, compiles -j1 and uses the frozen userspace.
Maximum diagnostic mesh is 3x129 nodes and 2048 macros; each response has at
most 160 iterations; there are exactly 12 slab cases and 32 neutral transport
steps. No new dependencies, full run, GPU work or field/source corrections.
