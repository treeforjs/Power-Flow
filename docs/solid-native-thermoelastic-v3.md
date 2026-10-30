# Native conductive heat and thermoelastic work, v3

This separate integration package binds the sealed electrode heat API v1 to
three-dimensional Q1 mechanics. Its dedicated signed mechanical cell-energy
channel represents reversible internal transfer. The native ABI accepts zero
Joule, authored, wall and species heat; it does not implement an EM solve.

The frozen MFEM archive is non-PIC. A persistent native executable links it
unchanged; the Python coordinator calls its C ABI through a framed, 17-digit
roundtrip protocol. No dependency is rebuilt.

The working native header is derived from a byte-preserved API v1 snapshot.
`vendor/pristine-manifest.json` records that baseline. The heat owner's live
source is never edited. Thermal material v2 is separately sealed and selects
q1 `nominal_1800`; the new adapter admits 298.2 K. Negative reference enthalpy
below 300 K is valid. Neither the old adapter nor the 773.15 K mechanics guard
is relaxed.

## Energy and mapping

The conditional constitutive energy is
`U = U0(T) + Phi(epsilon,T) - T * dPhi/dT` on the reference volume. For each
mechanical cell, the signed load to native heat is
`Qmechanical = Wstress - delta(Uelastic)`. Existing subcycle stress work and
kinetic-energy closure are checked before this load is admitted. Recoverable
strain energy alone is not used as thermal debit.

Native heat stores lumped nodal enthalpy and solves backward-Euler conduction
with harmonic conductivity frozen at acceptedtemperature. Each signed cell-J
is distributed equally to its four nodes. Full axial height is included once
in cell-J and divided once in the native two-dimensional load assembly.
The native per-node ledger checks
`delta(H) - Qmechanical + Qconduction_out = residual`.

Each reference XY cell owns one quarter of its four corner enthalpies. That
exact discrete native energy is distributed along the solid z partition by
reference length fractions, then inverted to a uniform mechanical cell
temperature. The first implementation requires identical complete XY grids.
It checks every mapped cell and the global sum, including signed cancelling
loads. Caloric inversion has an explicit representable-temperature roundoff
bound. This projection homogenizes subcelltemperature. Reaggregating nodal
loads does not reproduce the original cell debit; no false local cell
conservation identity is asserted. Local validation uses the native nodal
balance and the original mechanical cell balance, plus conservative transfer.

## Corrector and transactions

Each iteration starts from the same accepted heat state, displacement and
velocity. The native solver uses the current signed cell-J guess; mechanics
uses its trial caloric cell temperatures; a new work load follows. Acceptance
requires every cell-load mismatch and its signed sum to converge, native nodal
balance to pass, and the total change in caloric plus elastic internal plus
kinetic energy to close. No trial source is committed before that barrier.

Base state, revision, source, material, reference mesh, stage time and sequence
are bound at both APIs. Native transaction identities are deterministic from
the accepted state digest and next sequence. Discarded iterations can retry;
committed replay and stale trials fail. Root caller transaction labels are
included in receipts. There is no stochastic or species state in this ABI.

Joint checkpoints publish one directory containing the native checkpoint and
checksummed mechanical manifest. Reload validates an independent native handle
before replacing either accepted component. It binds the native executable,
material runtime, mesh, initial profile, constitutive reference and boundary
scenario. It is intended for the same qualified frozen runtime, not a portable
binary format. The worker acknowledgement can be lost after native commit:
the coordinator queries the accepted state and reconciles the deterministic
candidate. If it cannot prove the outcome, it stops that owner and requires a
fresh joint instance loaded from a checkpoint. This is a logical acceptance
barrier with atomic checkpoint publication, not a distributed crash-atomic
commit protocol.

## Demonstration and limits

`tools/run_native_thermoelastic_v3.py` uses a manufactured 2 by 1 by 2 micrometre
solid, a smooth temperature gradient near 299.5 K, an affine initial displacement matching the mean thermal strain, and explicit
initial expansion or compression velocity. Native conduction and signed internal work
are active with zero external heat. It checks discarded trials, accepted
steps, and bitwise-identical continuation after a joint restart. This is a
numerical binding benchmark, not a rerun of the 5 mm foil heating history.

All existing proof-stress, two-percent principal strain, positive-Jacobian,
material range and coefficient-knot guards remain. Piecewise coefficient
knot crossings are rejected because the current derivatives would create
artificial internal-energy/entropy jumps. The conditional Cp-as-Cv caloric
model is not a physical compressible EOS. The reference heat mesh is fixed;
there is no moving native mesh, deformed Maxwell solve, plasticity, ablation,
measured fixture, or production coupled breakdown claim.
