# Energy pathways, velocity distributions and Hall transport

Planning requirements for the finite-width desorption study. These define the
intended coupling; no water kinetic solver, Hall gap calculation or all-molecular
population solver is implemented by this document.

## Predict distributions before assigning thermal speeds

For each species and retained internal state, the transport quantity is a velocity
distribution, or particles sampling that distribution. Resolve drift, anisotropy,
thermal spread and reaction recoil separately. Only a Maxwellian population has
the usual single-temperature speed distribution. If a temperature moment is useful,
compute it after subtracting that species' local mean velocity:

`kB T_trans,s = (m_s/3) <|v-u_s|^2>`.

Report parallel/perpendicular temperature moments where the distribution is
anisotropic. A moment does not establish Maxwellian equilibrium, and rotational,
vibrational and translational temperatures are not interchangeable. The molecular
model evolves explicit internal populations rather than inferring all velocities
from a temperature fitted to band intensities.

## Energy and momentum at every transfer

1. **Electrode:** account for resistive heating, electron/lattice exchange where
   separately represented, conduction, material phase change and deposited particle
   energy. Desorption uses the physical surface/lattice temperature.
2. **Surface:** maintain finite adsorbate inventory and account for net surface
   energy exchange and emitted translational/internal energy. A rate activation
   barrier is not automatically the kinetic energy available to an emitted particle;
   do not subtract it again from the outgoing speed. Adsorbed and gas-phase energies
   require a consistent thermodynamic reference.
3. **Thermal launch:** for complete translational accommodation, sample the outward
   flux distribution `P_emit(v) proportional to vn exp[-m v^2/(2kB Ts)]`, `vn>0`.
   Its mean translational energy is `2kB Ts`; a volume Maxwellian has `3kB Ts/2`.
   Internal-state accommodation needs a separate model. Molecular rotation and
   vibration are not automatically equilibrated by this velocity sampling.
4. **Reactions:** conserve nuclei, charge, momentum and energy. Use species formation
   energies and level energies on a common reference. Initial/final internal states,
   dissociation/ionization energy, outgoing electron energies, fragment recoil and
   emitted/absorbed photons must all enter the event budget. Thresholds, binding
   energies and level energies must not be counted twice.
5. **Subsequent evolution:** collisions redistribute internal and kinetic energy;
   radiation, electron detachment and wall interactions have explicit sinks/sources.
   Retain nonthermal fragments until calculated collisions thermalize them.

For a specified two-fragment kinetic-energy release K in the product center-of-mass
frame, momentum conservation gives `K1=K*m2/(m1+m2)` and `K2=K*m1/(m1+m2)`.
Electron-impact fragmentation usually includes an outgoing electron and possibly
more products, so it needs the corresponding multi-body distribution. Ionization
also needs primary/secondary-electron energy sharing. A total cross section does
not determine any of these partitions. Record unavailable differential or branching
data as missing, with identified uncertainties or justified approximations.

Prompt H emission can arise from fast fragments. [Makarov et al.](https://doi.org/10.1029/2002JA009353)
measured distinct H(2p) kinetic-energy components from electron-impact H2O; those
condition-specific measurements are useful benchmarks, not universal product speeds.
Use an energy and momentum ledger per event and across the electrode/surface/gap
interfaces. Verify the accumulated ledgers independently of emitted spectra.

## Electrons and ions in strong fields

Both electrons and molecular/atomic ions obey

`dp/dt = q (E + v cross B)`, `dK/dt = q v dot E`.

Use relativistic momentum where the electron energies require it. Magnetic force
changes direction without directly doing work; inductive electric fields from a
time-varying magnetic field can transfer energy. Reactions create charged products
with their own mass, charge, internal state and momentum. Thereafter their orbits
respond to the fields; ionization does not instantly make a neutral follow a
prescribed drift trajectory.

Electrons and ions have different gyrofrequencies and gyroradii. Evaluate
`omega_cs=|q_s|B/(gamma_s m_s)`, `rho_s=p_perp/(|q_s|B)`, and the collision-to-orbit
and field/boundary scale ratios locally and by energy. Slower ion response does not
establish that ions and electrons are equally magnetized. Ions can be accelerated
by fields and inherit reaction recoil; their speed is an outcome, not a fixed
fraction of the electron speed. Neutrals have no direct Lorentz force but still
exchange momentum and energy with the charged populations.

Keep three velocity components even in a 2D spatial cross-section. For x=foil
width, y=gap normal, z=current, normal E and width-directed B produce an out-of-plane
E cross B drift. A 2D3V model retains that velocity but has no physical z end losses
under translational invariance. Quantify finite foil-length residence times with
an explicit validated closure or selected 3D calculations.

## Hall current and field ownership

In its applicable guiding-center limit, `v_E = E cross B / B^2` is shared by
magnetized species. Electron-ion relative motion produces Hall current. For a
quasineutral single-positive-ion plasma, an illustrative generalized Ohm law is

`E + ui cross B = eta J + (J cross B)/(e ne) - grad(pe)/(e ne) + electron inertia`.

The Hall contribution is nondissipative: `J dot (J cross B)=0`. It still changes
field evolution, trajectories, residence times and where electric work/collisional
energy deposition occur. Use a pressure tensor where the electron distribution
requires one. In the multispecies water plasma, especially with negative ions or
sheaths, form `rho_q=sum(q_s n_s)` and `J=sum(q_s n_s u_s)` explicitly. Do not impose
the single-ion quasineutral formula throughout the gap.

Separate-species kinetic motion plus self-consistent current/field evolution can
produce Hall behavior without an added particle force. Do not add an extra Hall
force to those orbits. A fluid region can use a consistent Hall closure instead;
the [FLASH extended-MHD guide](https://flash.rochester.edu/site/flashcode/user_support/flash_ug_devel/node107.html)
documents Hall support and wave/timestep constraints. Its plasma closure needs
independent justification before use in cold metal. The conductor and gap can
have different valid closures.

Assign one owner to each field update and to each contribution to charge, current
and energy deposition. The interface must distinguish applied, self-generated and
inductive fields, their time staggering, units and coordinate conventions. Do not
duplicate the same field, Hall response, Joule work or particle heating in two
solvers. Start with prescribed-field replay as a conditional benchmark; predictive
Hall feedback requires current/field coupling rather than particle tracing alone.

## Mesh and verification choices

The geometry/mesh frontend, FLASH grid and kinetic mesh need not be identical.
Define face normals/areas, material boundaries, interpolation, charge/current
deposition and conservative transfer explicitly. Cell size must resolve the scales
required by the chosen field/sheath model; particle timestep or a validated reduced
integrator must resolve the required orbit dynamics. An explicit PIC method and
an implicit/guiding-center method have different resolution requirements. A mesh
alone supplies none of the collision or orbit physics.

Before selecting reductions, benchmark pure-B orbits and energy conservation,
crossed-field drift, ion/electron mass response, collision energy loss, reaction
recoil, wall events, and a Hall/current/field comparison. Test transfer consistency,
finite-length losses, convergence and total energy/momentum/charge accounting.
The [MITL drift-kinetic comparison](https://doi.org/10.1103/PhysRevAccelBeams.26.090403)
illustrates useful reductions and the importance of particle-wall interactions.

Keep this energy/transport treatment identical between the finite-width source and
its release-rate-matched uniform counterpart. Otherwise changed Hall confinement
or fragment energy could be mistaken for evidence of spatially varying desorption.

Additional implementation references: [kinetic particle pushers](https://warpx.readthedocs.io/en/latest/theory/kinetic_particles.html),
[collision kinematics](https://warpx.readthedocs.io/en/latest/theory/multiphysics/collisions.html),
and [flux-weighted launch distributions](https://warpx.readthedocs.io/en/latest/usage/python.html).
