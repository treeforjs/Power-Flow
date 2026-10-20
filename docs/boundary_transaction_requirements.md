# Boundary transaction requirements for a coupled host

The [thermal boundary interface](temkin_boundary_interface.md) supplies immutable
release histories and weighted launch candidates. A host must still connect
those candidates to its live wall inventory, wall thermodynamics and gas state
through one accepted transaction. This document specifies that missing contract.
It does not implement wall cooling, moving surfaces or host transactions.

## Inputs and ownership

Before constructing a candidate, the host must identify the committed step and
wall state, selected source interval, reservoir/current/width patch, physical
area, exact emitted state and mass, and the proposed time interval. The source
and its temperature qualification, parameter provenance, state-energy reference,
geometry and supplied variates must be retained by hash or equivalent immutable
reference. A stored release history is evidence for a candidate; it is not
authority to deplete an independently evolved live reservoir twice.

The host must define the following quantities explicitly:

- Which wall inventory owns each adsorbed population, its nuclei and charge,
  and how emitted molecular states share that physical inventory. Splitting one
  adsorbate population among rovibrational states must not duplicate coverage.
- Whether an adsorbate's mass is already included in the wall mass and energy
  state, or is held in a separately conserved surface inventory. Removal must
  occur from that owner exactly once.
- The local material surface, outward normal, actual patch area and wall
  velocity at each emission event. A temperature sampled inside the material
  still requires the separately documented boundary-temperature qualification.
- The wall and gas thermodynamic references, available thermal energy and
  admissible post-release equation-of-state domain.
- Whether an external constraint supplies impulse/work to hold the wall fixed
  or at prescribed velocity, or whether finite-wall recoil is solved.

For a resolved planar source, physical counts are released number per area
times the actual current/width patch area. An additional current-length factor
would duplicate mass. A width-redistributed control moves those same launches;
the resolved and control arrays are alternative hypotheses, not two sources to
inject together.

## Conservative candidate ledger

For packet i, let wᵢ be its physical particle count, mᵢ its exact-state mass and
cᵢ its velocity in the instantaneous wall rest frame. Counts may be unequal and
noninteger. The host must use weighted sums, not the number of packet records:

\[
N_{r,p}=\sum_{i\in(r,p)}w_i,\quad
M=\sum_i w_i m_i,\quad
\mathbf P_{\rm rest}=\sum_iw_i m_i\mathbf c_i,\quad
K_{\rm rest}=\sum_i\tfrac12w_i m_i|\mathbf c_i|^2.
\]

The inventory decrease must equal these counts separately for every owning
reservoir and patch. Its associated mass, nuclei and charge transfers must
match the gas increments and any explicitly modeled surface reaction products.
Rounding tolerances must be stated in physical units and scaled to the actual
sum, with compensated or scaled accumulation where needed. Positive rare
components cannot be dropped to make a global residual small.

The source moment predicts a mean outgoing translational energy of 2kBT per
particle. Finite packet samples carry their actual K_rest instead. Preserve
both quantities and their sampling residual, but use the actual applied gas
energy in the accepted conservation balance. The residual is not another heat
term. A host that instead uses a moment-conserving representation must qualify
that representation separately; these packet interfaces do not rescale
velocities to enforce the analytic expectation.

## Adsorption and internal-energy references

The activation barrier in the release rate is a kinetic parameter. It does not
by itself specify the thermodynamic energy of the bound reservoir, even when
that barrier depends on coverage. The host needs an independently justified
adsorption energy reference or reservoir state function, together with emitted
state formation/excitation energies in a consistent reference.

For example, take gas-ground energy as zero and let removing a molecule from
its bound reference cost a positive D. For a constant binding cost, removal
raises the remaining reservoir reference energy by D, because a negative bound
energy has disappeared. Let I be the emitted internal energy above that same
gas-ground reference. With a stationary externally supported wall, the thermal
debit is

\[
Q_{\rm thermal}=K_{\rm rest}+D+I,\qquad
\Delta U_{\rm ads}=D.
\]

The total wall change is then −K_rest−I and the gas gains K_rest+I. D is an
internal redistribution between the wall's thermal and bound-reference
accounts, not an extra loss from the combined wall/gas system. A different
consistent reference changes the component entries, not the total balance.

For a coverage-dependent thermodynamic reservoir, use the actual state
difference U_ads(after)−U_ads(before), including any specified site or mixture
interactions. Do not infer that state function by substituting the Temkin rate
barrier for D. The supplied desorption and internal-energy parameters must
state whether their endpoints already include molecular formation, electronic
or rovibrational energy; the same contribution must not be added again by the
gas state catalog or wall equation of state. Missing adsorption/internal
references prevent a closed energy transaction even when the release rate and
thermal velocities are numerically well defined.

## Stationary and moving-wall momentum accounting

The current launch velocities have a wall-rest thermal interpretation. A
stationary constrained wall receives reaction impulse −P_rest; the external
constraint supplies the opposite impulse to hold it fixed. Its support work is
zero at zero wall velocity. A freely recoiling finite-mass wall requires a
coupled momentum and kinetic-energy update instead.

For a common prescribed wall velocity U, nonrelativistic frame conversion gives

\[
\mathbf v_i=\mathbf c_i+\mathbf U,\qquad
\mathbf P_{\rm lab}=\mathbf P_{\rm rest}+M\mathbf U,
\]
\[
K_{\rm lab}=K_{\rm rest}+\mathbf U\!\cdot\!\mathbf P_{\rm rest}
               +\tfrac12M|\mathbf U|^2.
\]

With different local wall velocities, apply these identities packet by packet;
one mean U cannot replace the weighted cross terms. Geometry and wall velocity
must be evaluated at the corresponding emission times. Moving-wall event
positions, swept areas and collision with a moving material interface also
require host definitions; a constant launch-plane coordinate does not provide
them.

A moving-wall host must define the departing adsorbate's pre-emission mass,
momentum and energy reference. Under one explicit prescribed-U convention,
that mass initially translates with the wall. Removing it subtracts M U and
M|U|²/2 from the host's advected momentum and kinetic energy. The relative
emission additionally transfers impulse P_rest to the gas. The support applies
the opposite reaction to maintain U, doing work U·P_rest on the wall/gas system.
With the rest-frame thermal debit above,

\[
\Delta E_{\rm wall}=-K_{\rm rest}-I-\tfrac12M|\mathbf U|^2,
\quad \Delta E_{\rm gas}=K_{\rm lab}+I,
\quad \Delta E_{\rm wall}+\Delta E_{\rm gas}
       =\mathbf U\!\cdot\!\mathbf P_{\rm rest}.
\]

This convention separates advected mass removal, recoil impulse and external
work. Simply debiting K_lab plus the rest-frame adsorption/internal terms from
wall thermal energy would mix those accounts. A freely moving wall has no
prescribed support work: the host must solve the remaining wall mass, velocity
and recoil energy consistently. No moving-wall implementation is supplied by
the present boundary modules.

## Atomic accept, rollback and replay

The coupled host must perform the following operations as one transaction:

1. Read a committed state/version and form a release trial from its live
   inventory and proposed wall-temperature history. Reject stale source or
   geometry associations before launching packets.
2. Produce immutable candidate launches from the supplied variates. Retain
   component counts, inventory changes, actual/expected energies, frame
   convention and all proposed mass, momentum, internal and thermal transfers.
3. Trial the wall energy/momentum update and gas insertion together. Check
   inventory nonnegativity, exact state membership, finite supported values,
   wall thermodynamic admissibility, event geometry and conservation residuals.
   Include modeled external work and reaction products in those checks.
4. Commit the wall inventory, wall mass/energy/momentum, gas particles and
   external-transfer ledger together. A durable transaction identifier must
   prevent inserting or depleting the same accepted candidate twice.
5. On rejection, discard every trial effect, including gas insertion, source
   depletion, wall updates and external-ledger entries. Restore the supplied
   variate position and committed step identity. Replaying the same candidate
   against the same committed state must reproduce identical arrays and ledger.

Reducing a rejected time step changes the release trial and conditional-time
law for that interval. Reevaluate both using a documented variate/component
association; do not keep emission times or inventory changes mapped for the
rejected interval. Restart records must distinguish accepted transactions from
pending candidates. Parallel ownership and reduction order must likewise avoid
double transfer at patch or rank boundaries.

Minimum acceptance exercises are a successful stationary transfer, insufficient
thermal energy, stale or depleted inventory, repeated transaction identifier,
an invalid state/geometry, rejection after provisional gas insertion, identical
replay, and restart immediately before and after acceptance. A prescribed-U
exercise must also close the separate advective, support-work and thermal
accounts. These are requirements for the integrating owner, not completed host
tests.

## Manufactured two-packet arithmetic acceptance case

Use one patch of area 10⁻¹² m², outward normal +y, and two distinct pools:

| Quantity | Packet/pool 1 | Packet/pool 2 |
| --- | ---: | ---: |
| Physical packet weight | 2 | 3 |
| Initial pool particle count | 10 | 8 |
| Remaining pool particle count | 8 | 5 |
| Exact-state mass (kg) | 2×10⁻²⁶ | 4×10⁻²⁶ |
| Wall-rest velocity (m/s) | (1000,2000,0) | (−500,1000,500) |
| Source temperature (K) | 300 | 600 |
| Independent constant bound-removal cost (eV/particle) | 0.2 | 0.4 |
| Gas internal energy above its ground reference (eV/particle) | 0.1 | 0.05 |

The initial coverages are (10¹³,8×10¹²) m⁻² and accepted coverages are
(8×10¹²,5×10¹²) m⁻². The weighted removed mass is 1.6×10⁻²⁵ kg. Direct sums give
P_rest=(−2×10⁻²³,2×10⁻²²,6×10⁻²³) kg m/s and K_rest=1.9×10⁻¹⁹ J.
Using the exact SI conversion 1 eV=1.602176634×10⁻¹⁹ J:

| Stationary-wall ledger | Energy (J) |
| --- | ---: |
| Bound-reference increase D=1.6 eV | +2.5634826144×10⁻¹⁹ |
| Outgoing internal I=0.35 eV | +5.607618219×10⁻²⁰ |
| Expected translational energy 2kBΣwT | 6.6271152×10⁻²⁰ |
| Actual minus expected translational energy | +1.23728848×10⁻¹⁹ |
| Actual wall thermal debit K_rest+D+I | 5.0242444363×10⁻¹⁹ |
| Total wall change, including +D | −2.4607618219×10⁻¹⁹ |
| Actual gas gain K_rest+I | +2.4607618219×10⁻¹⁹ |

A pre-step thermal account of 2×10⁻¹⁸ J remains positive at
1.49757555637×10⁻¹⁸ J. Subject to the other host checks, the candidate can be
accepted with the stated coverage and mass changes. Debiting the expected
instead of actual launch energy would leave an unbalanced 1.23728848×10⁻¹⁹ J.
With available thermal energy only 2.51212221815×10⁻¹⁹ J, this same candidate
must be rejected and every state must remain at its pre-step value.

For the separate prescribed-U accounting check, set U=(400,−200,100) m/s.
The lab velocities are (1400,1800,100) and (−100,800,600) m/s, giving
P_lab=(4.4×10⁻²³,1.68×10⁻²²,7.6×10⁻²³) kg m/s and
K_lab=1.648×10⁻¹⁹ J. The removed advected kinetic energy is 1.68×10⁻²⁰ J;
support work U·P_rest is −4.2×10⁻²⁰ J. The total wall change is
−2.6287618219×10⁻¹⁹ J and the gas gains +2.2087618219×10⁻¹⁹ J. Their sum
equals the support work input: the support extracts 4.2×10⁻²⁰ J.
The thermal debit remains the rest-frame value under this specified convention.

These numbers are independently calculated with exact rational arithmetic, not
sampled from a physical species or a calibrated heating case. The calculation
checks the ledger identities; it does not implement acceptance or rollback.
The recorded verification command was:

```sh
/workspace/Power-Flow/.venv/bin/python \
  /workspace/research_work/mitl-upstream-20261018T1858/boundary_transaction_arithmetic.py
```

It exited 0 at 2026-10-18T19:16:59.364069+00:00. The calculation source SHA256 is
`e8ead3f5edd50afd1c187073d7aa6b6bffd8f84359aa872e9058861fba079a68`;
its complete output is `boundary-transaction-arithmetic.log` in that verification
directory, SHA256
`e344e55239441b8f39561e06667c044a2ceb014f4aaa3c747210ca67c52bb27d`.
The command, working directory, timestamps, exit and log hash are recorded in
the same directory's `commands.jsonl`.
