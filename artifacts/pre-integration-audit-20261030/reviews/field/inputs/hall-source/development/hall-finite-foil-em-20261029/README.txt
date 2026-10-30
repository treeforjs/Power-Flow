Current-controlled finite-foil EM / native transport / thermal milestone

The axial backend now evolves total 3D E/B, native charged particles, and actual
Joule-driven conductivity/enthalpy feedback with a C19 TOTAL-current feed port.
Terminal voltage is predicted. Source or particle charge is never added to a
second native Poisson field. The original 2D prescribed-voltage cases remain
operator fixtures; they are not the selected physical drive.

Read native-binding.txt for entry points and ownership. Read
phase-and-source-integration.txt before connecting actual water SSA or source-v2.
The source-v2 midpoint/path-work gate is unresolved by that source package's
accepted-start wrapper and is not hidden by this milestone.

Geometry, discretization and drive
  Two 5 mm by 0.2 mm foils span z=-12.5..+12.5 mm, separated by a 5 mm gap. The pilot
  uses a distributed impressed current branch at the positive-z feed and an ideal
  PEC short sheet at negative z. Its steady short-return profile is explicitly
  the reflected uniform feed-contact quadrature. Experimental lead/short shape
  and its resistance/heat are not supplied or inferred.

  A 3D Yee incidence complex and positive SI Hodge volumes retain all derivatives.
  Ex/Ey/Ez currents feed Ampere; Faraday evolves Bx/By/Bz. The scalar potential
  has phi=0 at the short and includes all three gradient components. Initial
  magnetostatics uses tangential A=0 on the finite exterior with interior Coulomb
  tests. No sampled-Az/backward-rate replay or extra Hall velocity is used.
  The physical core B orientation is preserved; genuine end fields are allowed.

  The first supported C19 knot is 7.49788055492e-10s,6012.07136047 A. The chosen
  initialization is a uniform 300 K DC scenario, not recovered experimental
  prehistory. The sealed measured trace is interpolated inside original support;
  prehistory and post-support extrapolation reject. Clock registration remains
  an experimental uncertainty. The conditional voltage pulse is not imposed.

  Current is prescribed only at the feed branch. Solid, plasma, displacement,
  transverse transfer, charge accumulation and outer terms respond through the
  same field solve. Terminal voltage is minus the E line integral along the unit
  feed path, positive lower-minus-upper. Source power is I*V once. An impedance
  boundary has explicit outgoing loss and a separate static matched-start
  boundary supply. It does not prescribe an evolving background B independently.

Native full-time transport and energy
  The one canonical arena supplies a charged GapParticle view. Only native
  pic_boris supplies momentum. The qualified work displacement and line-integrated
  gather/deposition are a distinct orbit map; they must not be confused with the
  old GapPic D-K-D phase. Its full-time x/u endpoint is explicitly typed and permits
  variable dt. Layer-owned IDs are excluded from ordinary motion.
  All three current components and endpoint charge use exact piecewise polynomial
  path quadrature. The nonlinear midpoint field/particle iteration replays the
  same immutable accepted arena; no duplicate accepted inventory exists.

  Chemistry(dt/2), EM(dt), chemistry(dt/2) is qualified with manufactured angular
  operators at physical stages t_n and t_n+dt. SSA waiting time remains local
  exposure; positions cannot move during chemistry. Source tokens, neutral
  anchors, energy, momentum and nodewise charge have explicit admission checks.
  This demonstrates the composition and phase contract, not actual water rates
  or actual-water SSA accuracy. Root owns those seals and their combined tests.

  Field updates retain canonical high+low words so sub-ULP changes of the large
  DC field survive. Both words are checkpointed and included in physical energy;
  no corrective energy is added. An independent 80-digit quadratic-energy check
  reconstructs the SI Hodge weights and agrees with the work ledger to 9.57e-26 J.

Thermal coupling
  The copied sealed native v1 thermal owner uses evaluated SS304 solid properties
  on 300..1500 K, rho_ref 8000kg/m3, and an explicit fixed-volume Cp-as-Cv approximation.
  One owner per foil and axial slice receives the SAME actual midpoint Joule
  power that debits the field. Height=dz is applied once. Accepted temperature
  changes the conductivity assembled for the next EM step. This lagged thermal
  coupling is not claimed generally second order.
  Axial heat diffusion, moving geometry, short-sheet heating and vapor are absent.
  The root's source-v2 ABI/material/source owners are separate; do not substitute
  this v1 copy or compile a second source kernel. Existing mechanics/data limits
  remain in force.

Accepted evidence
  finite-final-002:70 checks for2D fields, native transport, surface paths and
    ownership/rejection contracts.
  axial-phase-002:72 checks for3D fields, feed/short, actual native particles,
    actual thermal feedback, variable-step composition, source path coordinates,
    rejection atomicity, refinements, fresh-process restart and independent energy.
  New milestone total 142 checks; prior periodic 321 and initial 162 checks remain
  preserved in their earlier source/evidence milestones.

  The 32-step10 fs native electrothermal pilot spans 0.32 ps with an authored
  electron/proton seed density 1e12 cm^-3. It transfers 4.15812104e-7 J of Joule heat,
  raises maximum T by 2.1590913 microK, and changes sigma by 0.0034247362 S/m. The
  changed sigma alters the actual next Maxwell solution by 2.1234e-5 V/m. Maximum
  B change is 1.47336e-6 T. Maximum particle/field work defect is 1.48e-22 J.
  The thermal numerical residual is reported separately; the combined ledger's
  maximum step defect 2.97e-15 J is dominated by native temperature roundoff.

  dt 10/5/2.5 fs refinements of the same 0.16 ps interval reduce successive E/B
  differences by about4.0. Variable chemistry/EM composition tests range 10 fs to
  0.625 fs and reduce the proper-velocity difference by 4.00032 under subdivision.
  Ten checkpoint files, including field low words and native thermal owners,
  reproduce exactly in a fresh process. Root still owns full source/RNG/event-
  journal and crash-consistent whole-driver checkpointing.

Physical gates
  Spatial sweeps use nx 45/90/180 at fixed 9 mm width, axial 4/8/16 cells with fixed
  6.25 mm air padding per end, and 9/17/25 mm boxes at 0.2 mm transverse spacing.
  Sampled central Bx is 0.492345,0.690871 and 0.729353 T across 9/17/25 mm boxes. The
  17-to25 mm difference is 5.57%; these small exteriors are NOT converged. Axial
  refinement also remains required. Reuse/extend the root's wide graded exterior
  rather than reverting production to the small operator box. No axial-instability,
  experimental-breakdown or long-run readiness is claimed.

  Actual source-v2 layer integration still needs a trial-field stage preserving
  accepted thermal-time provenance, conjugate path work, canonical exclusion,
  handoff remainder and collection/source ledgers. The reference split coordinator
  rejects active layers so they cannot accidentally receive ordinary motion.
  Root owns the full actual-water/source/arena coordinator and GPU/long-run release.

Reproduce under the existing resource agreement
  From the private source checkout, run each command sequentially:
  env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CUDA_VISIBLE_DEVICES= taskset -c6 nice -n15 flock /home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-profile/testing-core.lock python development/hall-finite-foil-em-20261029/run.py NEW-2D-RUN
  env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CUDA_VISIBLE_DEVICES= taskset -c6 nice -n15 flock /home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-profile/testing-core.lock python development/hall-finite-foil-em-20261029/run_axial.py NEW-AXIAL-RUN
  Each build/run uses frozen userspace/dependencies, CPU 6/nice 15, one thread,
 2 GiB,180 CPU seconds and240 wall seconds per native command. New run names prevent
  overwriting evidence. No GPU, dependency upgrades, shared-repository edits,
  old-run signals, source thinning, or numerical-limit relaxation occurred.
