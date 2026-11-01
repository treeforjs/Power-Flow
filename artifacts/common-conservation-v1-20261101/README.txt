Independent field/thermal stock operands v1

Qualified component; full system and production physics remain unqualified.

field_thermal_assessment(actual Hall Owner, actual JointProposal) evaluates
compensated endpoint EM energy using exact products, actual Q1 nodal thermal
energy under the pinned nominal SS304 material, and terminal drive/exterior
work from the actual midpoint operands. No reported component residual is
summed. The reconstructed nodal ordering must match the accepted Coupler's
public stored_energy_J() bit for bit before candidate evaluation is permitted.
Reference-enthalpy change and the stable material increment are both retained;
their difference is reported explicitly. A lower-level FieldHeatProposal
reader requires the actual declared dt, not end minus begin (which can round).
Root must retain/preflight the actual Hall publication capability; matching
configuration/revision strings alone do not confer authority.

SimQueue671: 27 normal and 27 UBSan controls with deterministic replay. Actual
Hall fixture: eight thermal slices, 416 nodes, frozen userspace, one explicitly
borrowed core4, preemption freeze/resume respected. This is functional
qualification, not a timing benchmark. Independent vacuum balance was
2.962157776198373e-15 J within the unchanged 6.568373383520938e-12 J allowance.
The source/material/gas combined bound is implemented separately in
common-material-stocks-v1-20261101; this component alone does not establish a
populated water-plasma conservation result.

Controls cover signed low-word EM cross terms, exact tiny terms, incomplete
state/nonfinite/negative/inexact-subnormal metric refusals, nonuniform actual
material temperatures, real nodal topology and evaluation-order parity,
material support limits, actual Hall proposal and read-only endpoint stocks,
actual terminal/exterior fluxes, and wrong declared interval refusal.

Initial harness failures659/663 were link/include omissions; their logs are
preserved as historical setup evidence. Final qualification uses source and
input pins in INPUTS.json. No external runtime or running simulation was
modified. Numerical binaries, objects and frozen dependencies stay outside git.
