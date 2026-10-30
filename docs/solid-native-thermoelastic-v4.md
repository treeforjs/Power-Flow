# Signed mechanics on the sealed heat/source v2 API

The separate v4 package rebases the dedicated signed mechanical cell-energy port
onto the sealed `electrode_heat_v2` API. V3 and the heat owner's live source remain
unchanged. `vendor/pristine` preserves the complete sealed v2 snapshot; the exact
source kernel is compiled once. The working material callbacks, surface binding,
water law and their seals remain byte-identical. The v4 executable links the
existing frozen MFEM archive without rebuilding dependencies.

## Combined energy and state ownership

Native heat receives the signed mechanical transfer
`Qmechanical = Wstress - delta(Uelastic)`, with
`Uelastic = Phi - T*dPhi/dT`. The thermal source input is
`Qexternal = Qwall,e - Kwater - Uinternal,water - delta(Uads)`.

The certified electron equation stays
`Qwall,e + delta(Ke) + delta(Ref_e) + delta(Image) = Wfield`.
Only `Qwall,e` enters the solid heat equation. Field work remains a separate
weighted current-path quantity. It is not Joule heat or a material reference
change. The source's particle ledger is birth metadata; continuation charges
only its returned incremental exchange.

Gas H2O reference energy remains zero. The adsorbed state owns the negative
reservoir `Uads(N) = -Eb*q*N*(1-alpha*N/(2*N_ML))`. Removing actual weighted water
raises that reservoir by the sealed finite-step difference. Native heat loses
actual gas kinetic energy, explicit internal energy and this reservoir increase
once. No activation-barrier, mean `2kT`, extra transported reference or invented
recoil heat is added. Canonical hi/lo adsorbed and gas counts are preserved.

For stationary homogeneous supports, mechanics checks
`delta(Ksolid) = -Wstress`. The joint first law is therefore

`delta(H) + delta(Uelastic) + delta(Ksolid) + delta(Ke) + delta(Ref_e)
 + delta(Image) + Kwater + Uinternal,water + delta(Uads) - Wfield = 0`.

Each native node independently checks
`deltaH - mechanicalLoad - electronWallLoad - waterLoad + conductionOut = residual`.
The three input channels are assembled separately, including tiny cancelling
source contributions. Conservative reference-cell mapping remains the exact
native lumped corner-enthalpy partition. The corrector validates each cell's
signed mechanical load and the global sum, using explicit local thermal
representability bounds. Independent tests reconstruct energy from mechanical
tensors, native temperatures, proper velocities, weighted paths and canonical
water inventories rather than relying only on reported residuals.

## Common preparation, acceptance and restart

Each corrector iteration prepares all source actions and water draws from the
same accepted-start temperature and explicit total E/B stage. The source
candidate must remain bit-identical while the mechanical heat guess changes.
One native participant owns heat, each bound electron layer/origin/status, gas
packets, adsorbed inventory, weighted current paths and cumulative charge and
energy terms. It retains all nineteen v2 thermal cumulative fields and keeps
mechanical exchange separate. Active layers must reach the common accepted end
clock; ordinary PIC is not advanced by this component. Same-face water patches
cannot overlap and duplicate the physical reservoir.

A joint checkpoint stores that complete native payload with displacement,
velocity, initial mechanical identity, proof-screen identity, reference-mesh
and material/source hashes, and the last consumed source stage/actions/draws.
Restart validates a fresh native participant and mechanical state before either
accepted component is replaced. Candidate and accepted native digests cover
all source/gas/path state. Lost commit acknowledgements are reconciled with the
complete candidate digest; an unprovable outcome stops the owner and requires
checkpoint recovery. This is a logical acceptance barrier and atomic common
checkpoint publication, not a distributed crash-atomic commit protocol.

## Admission and component evidence

Mechanics remains within 298.2–773.15 K, with the unchanged nominal proof-stress,
2% principal-strain, positive-Jacobian and positive-capacity guards. Thermodynamic
coefficient-knot crossings remain rejected. Reviewed water with remaining
inventory requires both accepted and candidate patch temperatures in its
300–1300 K source range. The 298.2–300 K interval requires the explicitly selected
conditional analytic water continuation; no source bound is silently relaxed.

The combined demonstration uses actual sealed source begin/advance calls and
actual oxide-water transfers on disjoint reference-surface patches, along with
native conduction and Q1 mechanics. Its explicit electron weight is a component
event input, not a predicted count or current. Prescribed one-sided fields do
not constitute a field solve. The clean-electron and oxide-water patches are
separate component hypotheses, not a measured SS304 microscopic calibration.
Joule and authored heat remain zero.

The source surface is frozen and stationary in this approximation. Its actual
emission momentum is retained as a fixed-support impulse diagnostic; it is not
applied to the evolving foil as an invented recoil force or kinetic heat. There
is no moving-wall source coupling, deformed field solve, plasticity, ablation or
physical compressible EOS. Same-particle boundary collection without intervening
PIC transport is an ownership test, not a claimed return trajectory. Actual
layer recollection is exercised separately through the sealed layer dynamics.

The actual evolving-material pilot is audited separately. Its accepted upper
and lower native checkpoints, physical coordinates and node permutation are
bound before selecting an early supported prefix. A new nominal-density scenario
matches its 8000 kg/m3 reference assumption without changing the baseline SS304
record. Prescribed-temperature geometry comparisons retain explicit support
brackets and do not infer missing short/joint stiffness or claim mechanical
feedback was present in that pilot.

## Qualified evidence

The frozen CPU7 run passed 166 native v4 checks, all 8,731 sealed-v2 regression
checks against the working derivative, 11 new joint Python checks, 418 existing
Python checks, and 8 actual-pilot projection/admission checks. The signed
four-step component histories each passed complete rollback and bitwise restart
continuation. Their maximum joint residual was 4.4072e-21 J and maximum native
nodal residual was 1.3538e-23 J; both had nonzero actual electron and water
transfers with zero Joule/authored heat.

The observed material-feedback pilot supplied 244 accepted two-foil states.
Its first open constitutive-segment prefix contains 60 states, ending at
12.428157577 ns. All 59 intervals passed the conditional mechanics comparison
on the 32x4x2 mesh per foil. At the last stage, minimum gap reductions were
22.541192 nm for free foils and 23.284622 nm for both-end-clamped foils. Maximum
sampled von Mises stresses were 32.811453 MPa and 27.107708 MPa, respectively;
the nominal proof ratios remained below 0.211. Each foil retains 0.0002 kg of
reference mass. These numbers are a coarse prescribed-temperature comparison,
not spatial convergence or a deformed electromagnetic solve. The observed
checkpoint and executable bytes are hashed; a formal owner source-to-build/run
seal is not established. Raw checkpoints, coordinate certificates, restart
states and detailed limits are preserved with the v4 handoff.
