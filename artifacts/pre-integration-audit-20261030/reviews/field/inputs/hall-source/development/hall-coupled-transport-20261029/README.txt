Hall-coupled transport: bounded CPU prototype, 2026-10-29

This directory implements a manufactured bridge from signed kinetic currents to
self-consistent electron/ion Poisson response and a linear electromagnetic Hall
mode. It is a reviewable research milestone. It does not reproduce a CHICAGO
application, couple the heating history to plasma magnetic feedback, or supply
water reaction data. All tests start from their declared manufactured state.

Source base: treeforjs/raikou a7e44f7e9c3f7bd7d653d7e5a31cf3d850ee8e7d.
The only changes outside this directory are add_authored_batch and a diagnostic
force-midpoint observer in include/raikou/mfem_gap_pic.h and src/mfem_gap_pic.cpp.
The existing relativistic Boris force law, Q1 deposition/gather, Poisson solver,
residual tolerances and boundary treatment are retained. The observer is empty
by default and must not throw or reenter the solver. Raw native checkpoints
retain their existing ABI; the adapter adds a separate, versioned envelope.

Minimum formulation and ownership

For each charged species s, advance the existing kinetic equations:
  dx/dt = u/gamma,  du/dt = (q_s/m_s) [E + (u/gamma) cross B].
  rho = sum_s q_s n_s,  J_s = q_s integral(v f_s dv),  J = sum_s J_s.
The electrostatic stage solves one common height-averaged 2D Poisson RHS for
ALL signed species before the next force gather. Its trajectories have three
velocity and position components. The manufactured wave is invariant in x,z,
so its out-of-plane current is meaningful in this restricted geometry. Spatial
self-field variation along z is absent from the existing application.

No extra Hall velocity belongs in these orbits. For singly charged, quasineutral
species J = en(ui-ue). The electron momentum equation gives
  E + ui cross B = J cross B/(en) - div(Pe)/(en)
                  - (me/e) D_t ue + Re/(en).
Here Re is the momentum transferred INTO electrons by collisions. The Hall term
expresses electron/ion relative motion; it is not a second force or drift to add
to Boris. When both species share a balanced E cross B drift, their equal-density
charge currents cancel. During a cold start, magnetized electrons and much slower
ions can carry a signed transverse relative current. This transient is the first
control here; it is not a collisional steady conductivity measurement.

A fully electromagnetic stage also needs
  epsilon0 dE/dt = curl(B)/mu0 - J,  dB/dt = -curl(E),
  div(epsilon0 E)=rho, div(B)=0,
with compatible charge/current deposition, field interpolation and boundary or
circuit work. Adding mobile ions to Poisson alone does not implement these
induction equations. For a split prescribed background plus plasma perturbation,
the background's Maxwell residual and supplied energy must remain explicit ports.

Energy is owned once: the PIC particles own electron and ion translation; the
zero-bias Poisson control owns Usc; the Maxwell control owns electric and magnetic
perturbation energy. Constant prescribed B contributes no direct work. All
populations here are externally authored, with no collision reference-energy,
neutral, material or radiation reservoir. Initial field energy is part of the
initial state. No source or loss occurs in the closed wave after initialization.
Wall tests count incident kinetic energy once and do not interpret it as a model
of deposited material heat.

At the actual Boris midpoint, the observer uses the total E already gathered.
Let a = dt*q*E/(2m), u_minus=u0+a, u_plus=u1-a. It records
  D = dt/2 [(u_minus+u0)/(gamma_minus+gamma0)
            +(u1+u_plus)/(gamma1+gamma_plus)],  W = w*q*E dot D.
It audits both signed and absolute accumulated (Delta K-W) defects. In the
zero-bias wave this is self-field work. It does NOT provide the missing separate
bias/induction/axial-correction work ledgers for the production field adapter.
No numerical residual is deposited as physical heat.

Executable controls

1. Signed E cross B current, physical proton/electron mass ratio.
   E=(1e5,0,0) V/m, B=(0,0,+/-10) T; n_e=n_i=1e8 cm^-3.
   Both species start cold. Average over 20 electron gyroperiods, approximately
   71.45 ps, while the ion gyro angle is only 0.06844 rad. Compare with the cold
   nonrelativistic analytic average -E/B [1-sin(Omega_s T)/(Omega_s T)]. The actual
   native push is relativistic; v/c is small enough for that reference here.
   Three time resolutions retain dt<=10 fs and gyroangle<=0.05. A separately
   balanced common-drift launch checks zero relative current. Reversing B reverses
   the signed currents. Neither result supplies an instability growth rate.

2. Native self-consistent magnetized Poisson wave.
   A grounded 1 mm x 1 mm x 1 mm extruded box; nx=2, ny=16/32/64. No wall loss.
   omega_pe=1e12 rad/s, n_e=n_i=3.142077833705378e14 cm^-3, Bx=5.6856301113 T.
   Quiet-start electron displacement delta*sin(2*pi*y/L), delta=1 nm, ions
   initially undisplaced. Ion mass ratios 25 (manufactured) and 1836.15267
   (physical proton/electron ratio) are separately labeled. No particle thinning.
   Uniform physical density is held fixed as samples/cell are increased 2/4/8.
   E_y produces y motion and z Lorentz transport. Duration 4 ps; dt=10/5/2.5 fs.

   The independent small-amplitude cold solution uses tau=omega_pe*t, M=mi/me,
   b=|Omega_e|/omega_pe, and displacement amplitudes xi=(xi_e,xi_i):
     xi'' + A xi = (b^2*delta,0),
     A = [[1+b^2,-1],[-1/M,1/M+b^2/M^2]].
   The two eigenvalues of this 2x2 matrix give cos(sqrt(A)*tau) exactly.
   vz_s=-Omega_s*(xi_s-xi_s(0)); E_y=en*(xi_e-xi_i)*sin(k*y)/epsilon0.
   phi=en*(xi_e-xi_i)*(cos(k*y)-1)/(epsilon0*k) satisfies both grounded ends.
   Nonlinear and relativistic corrections are small for the chosen amplitude;
   mesh, time and sampling convergence are tested independently. Maximum sampled
   waveform errors supplement the final-time errors to avoid phase cancellation.

   Signed species Q1 endpoint charge and two-half-drift edge transfers are audited
   independently against each species' continuity equation and the native total
   path transfers. Errors use gross absolute charge, not a nearly cancelling net
   charge. Currents are velocity moments with signed q; species density uses
   physical number weights. Source/boundary number is not inferred from net charge.

3. Cold two-fluid Maxwell Fourier mode, physical mass ratio.
   This is a linear, single-mode, periodic 1D model, not a full EM particle solver.
   It uses MFEM DenseMatrixInverse for an implicit-midpoint 8-real-variable solve.
   The spatial operator is the Fourier symbol of a staggered periodic curl,
   k_h=2*sin(k*dx/2)/dx, with 16/32/64/128 samples/wavelength; the full spatial grid
   is not stored. These restrictions and the Fourier invariant subspace are
   essential to interpreting the refinement result.

   Normalize time by omega_pe, distance by de=c/omega_pe, E by me*c*omega_pe/e,
   B perturbation by me*omega_pe/e, velocity by c. For helicity h=+1 the complex
   amplitudes represent Ey-iEz, By-iBz and vy-i vz; h=-1 reverses polarization.
   X=(E,B,ve,vi) obeys
     E'=h*k_h*B + ve-vi; B'=-h*k_h*E;
     ve'=-E-i*h*b*ve; vi'=E/M+i*h*b*vi/M.
   This includes displacement current, mobile ions, electron inertia and both
   Maxwell curl updates. The equilibrium Bx has b=0.5 (2.8428150557 T) at
   omega_pe=1e12. Initial |E|=1e-6, delta B/B0 approximately 2.07e-5.
   rho and div J are identically zero within this transverse linear subspace;
   the zero continuity result here does not qualify a general EM deposition.

   The independently solved continuum dispersion relation is
     omega^2 - k^2 - omega*[1/(omega-b)+1/(M*(omega+b/M))] = 0,
   using the positive branch below b. Midpoint has exact discrete frequency
     omega_num=2/dt*atan(omega(k_h)*dt/2).
   The Hall-MHD limit with cold mobile ions, negligible electron inertia and
   negligible displacement current is
     omega_H=.5*[b*k^2+sqrt(b^2*k^4+4*b^2*k^2/M)].
   The difference decreases as k*de decreases. At k*de=0.2 it is about 5%, so
   replacing the full model by that limit is not silently treated as exact.
   Field+species energy is (|E|^2+|B|^2+|ve|^2+M|vi|^2)/2. Midpoint's electric
   work uses the SAME midpoint J as Ampere; magnetic rotation contributes zero.
   Disabling J in Ampere or disabling Faraday deliberately breaks the coupled
   physical model; their changed mode and energy defects are negative controls.

Species, authored products and restart

coupled_pic.hpp owns a checked species registry and ID-keyed event/reference
metadata. An authored batch validates all products, stages allocations/ledgers,
and commits before any solve. An invalid second product leaves particles, IDs,
ledgers and cached fields unchanged. A valid collocated positive-ion/electron
pair cancels nodally at insertion; an electron-only intermediate solve has a
spurious nonzero field and is retained as a negative control. This is the charge
transaction needed for an eventual ion/secondary event, not an ionization model.

The separate hall-authored-pic-envelope-v1 records species/event IDs, initial
reference positions, accumulated work and activity, the immutable registry and
native checkpoint checksum. It checks the same executable on restart. FNV-1a
checksums detect mismatched files; SHA256 source, dependency and executable hashes
are additionally saved in the run manifest. They are provenance checks, not an
adversarial file-authentication system. Restart must use a fresh adapter and
retain its envelope. Native checkpoint callbacks are reattached explicitly.
Wrong species registries and corrupted native files are rejected. The split-run
continuation is compared with an uninterrupted run, including work ledgers.

The prototype's mutable native() access is for this controlled harness; application
integration must prevent insertions that bypass the sidecar owner. The authored
batch API intentionally cannot debit neutral/primary particles or label its
kinetic energy as a chemical transfer. A future apply_reaction_batch must stage
all reactant/product weights, source provenance, IDs, reference/internal/recoil
energy, RNG state and field invalidation together. It must debit a parent once
and create secondary charge and ion charge in the same accepted transaction.

Reproduce using the pinned machine paths (fresh output name each time)

flock -n /home/treefor/sim-workspace/powerflow/outputs/pic-profile/testing-core.lock taskset -c 6 nice -n 15 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CUDA_VISIBLE_DEVICES= PYTHONDONTWRITEBYTECODE=1 /home/treefor/projects/runtimes/arch-base-20261028/run-frozen /usr/bin/python /home/treefor/projects/runtimes/hall-coupled-transport-20261029/source/development/hall-coupled-transport-20261029/run.py run-new

run.py compiles against the frozen MFEM/HYPRE/METIS/GKlib/MPI CPU libraries, runs
the benchmark and original PIC foundation suite, applies controls.json limits,
and renders white-background plots with the frozen Power-Flow Python environment.
No CUDA library is linked by these commands. It uses one CPU6 job, nice15, one
thread, an exclusive shared lease, 2 GiB address-space and 180 CPU-second limits
per process. The populated step cap, gyro/plasma/displacement guards, and native
Poisson residual guards are retained. controls.json was written before the first
run; its numerical limits were not relaxed. Failed/early run folders are retained.

Results and intended integration are in the accompanying review.txt and
next-integration.txt. The complete raw waveform JSON, numerical checks, commands,
input hashes, executable hash and checkpoints are in the reported run directory.
