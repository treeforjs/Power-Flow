# V4 source/heat/mechanics thermodynamic review

This independent source review treats v3 as immutable. It reviews the sealed native heat/source v2 package and its proposed v4 integration. Tests are prepared for the central CPU7 owner; this reviewer has not run numerical tests or builds.

## Source identity and admitted scope

The reviewed seal is `pic-electrode-heat-feedback-20261029/source-seal-v2.json`, SHA256 `6a802752324a697ff2e6ede68e088ddccd13e7983158a754da5df513b7ddac32`. The five requested snapshot files were independently hashed and matched their seal entries:

| File | SHA256 |
|---|---|
| API-v2.txt | b1d24b5f51feec3022a2a46cc6d132549251482757387b722dea47483eda716d |
| electrode_heat_feedback.hpp | 8e53182362af46ee991e87df668f1fc7959dcc47adfca56f4ab4e379ae165fc5 |
| surface_source_binding.hpp | 2439d148c2d2fa55902a29880a14750bd37231c468ac3b775578eea089d82482 |
| water_scenario.hpp | 99e82198fbc1766c29fbaf4a196118f322d45e70b1b2007ffa6f414656833489 |
| test_binding.cpp | f1062fd3b9c6f4b7a64ecc4e38266cebbc35bed4476c72050803215a39538bdd |

The thermal owner provides the reviewed 298.2 K cold branch with reference energy zero at 300 K. Mechanics independently remains restricted to 298.2–773.15 K in this integration and must stay within one open smooth constitutive segment for the conditional internal-energy port. This is not a physical compressible EOS qualification.

The water source is a conditional Fe2O3(0001)-like surface surrogate. Its selected parameters and Eyring prefactor are not a measured SS304 surface calibration. With remaining water, `Reviewed300To1300` intersects the mechanical support at 300–773.15 K. Using 298.2–300 K requires the explicit `ConditionalAnalytic298To3000` alternative and its unchanged ideal oxide hypothesis. The sealed source evaluates accepted-start patch mean temperature; trial-stage admission should also prevent acceptance of an unsupported remaining-water state. Mechanical limits apply to actual native nodes, not just source patch means.

The electron law additionally needs explicit material phase, locally planar barrier, and equilibrium electron-supply admission. Its E/B and source temperature belong to the same accepted-start stage, with exact field/common/thermal/mesh/material identities. Clean 4.3 eV electron emission and the oxide-water surrogate are a jointly declared scenario, not one measured microscopic surface model.

## Independent combined energy derivation

All quantities below are already weighted, full-height joules. Positive heat is into the solid. Let H be base caloric reference energy and let Uel = Phi - T Phi_T be the conditional elastic internal energy. Let Ks be mechanical kinetic energy. Actual quadrature/substep stress work is Wstress.

For stationary homogeneous supports or free mechanical boundaries:

    Delta Ks = -Wstress
    Qmech = Wstress - Delta Uel

The certified electron source supplies five distinct operands:

    Qwall,e + Delta Ke + Delta Re + Delta Image = Wfield

Only Qwall,e is deposited in heat. Electron kinetic change, reference change, image change, and field work are distinct owner ledgers. They must not also be entered through legacy incoming/outgoing/Nottingham/authored heat fields.

For water, with gas reference zero:

    Qwater = -(Kwater + Uinternal,water + Delta Uads)
    Uads(N) = -Eb0 qe [N - alpha_ref N^2/(2 N_ML)]
    Delta Uads = Eb0 qe dN [1 - alpha_ref (theta_before + theta_after)/2]

Removing water makes the negative adsorbed reservoir less negative, so Delta Uads is positive. It is a reservoir-state change, not gas reference carried by the packet. `packet.reference_J` stays exactly zero. The activation barrier governs rate only and is not another energy debit. The diagnostic mean 2 kT is not actual packet kinetic energy and is not charged.

Native heat therefore obeys

    Delta H = Qmech + Qwall,e - Kwater - Uinternal,water - Delta Uads + QJoule

where QJoule is present only when a separately admitted electrical owner supplies it. In the zero-Joule v4 pilot it is zero. Substitution gives the independently checkable joint residual

    Rjoint = Delta H + Delta Uel + Delta Ks
           + Delta Ke + Delta Re + Delta Image
           + Kwater + Uinternal,water + Delta Uads
           - Wfield - QJoule.

The intended result is Rjoint = 0 within separately reported native, mechanical, source, mapping, and arithmetic bounds. Do not add Qwall,e a second time to this final expression. If a field solver owns electromagnetic energy, its energy change and external circuit work must close with Wfield once. With prescribed accepted fields, Wfield is external work in this subsystem balance; the pilot cannot claim whole-driver closed-system conservation.

## Local identities and spatial scope

The native seven-ledger interface has the correct local sign convention:

    deltaH_node - mech_node - external_node + conduction_out_node = residual_node
    external_node = electron_node + water_node.

Summed electron nodal load must equal sum Qwall,e. Summed water nodal load must equal minus the sum of actual packet K, internal energy, and adsorbed reference change. Summed mechanical nodal load must equal the signed full-height cell-J input once. Conduction cancels globally for insulated boundaries, within arithmetic error.

The mechanical cell identity is Delta Uel_cell - Wstress_cell + Qmech_cell = 0. This and the native nodal heat identity do not establish a common cellwise total-energy PDE balance: nodal-to-cell caloric averaging homogenizes the mechanical temperature, and mechanical energy can cross cell boundaries through stress power. Local residuals and the separate global balance should be reported with their actual discretizations.

## Absolute state and checkpoint ownership

The accepted total subsystem energy can be reconstructed as H + Uel + Ks + Ke + Re + Image + Kgas + Uinternal,gas + Uads, with initial reference offsets retained. Subtract accepted accumulated external field work for a conservation diagnostic. Do not reconstruct tiny changes by subtracting two large scalar totals alone: use stable per-transfer Delta Uads, canonical hi/lo counts, compensated accumulation, and an explicit absolute-state reconstruction bound.

There is one parent electron population. An Active source-layer electron is excluded from ordinary PIC advance. A PICHandedOff electron is transferred once with current proper velocity and physical origin plus layer displacement; its image tail is restored to zero already. Recollected removes the same electron. The birth particle ledger is metadata; only returned incremental exchange is charged on continuation.

For an independently reconstructed current layer state:

- Active: Ke comes from current weighted relativistic proper velocity; Re = weight qe phi; Image = -weight qe^2/(16 pi eps0 x), x the outward normal distance.
- PICHandedOff: Ke is current; Re remains weight qe phi; Image = 0 after tail restoration.
- Recollected: the current particle, electron reference, and image contributions are zero.

Checkpoint native thermal state, source-layer payload/origin/status/identity, physical charge/reference/image ownership, canonical adsorbed and gas counts, actual gas packets, weighted current paths, accepted event watermarks, and external work in the same transaction as mechanical u/v. Current paths and LayerExchange values are already weighted. Zero-duration quantum bridges retain displacement/current/work without an invented flight time. Retry or corrector iteration must reproduce the same source candidate from accepted-start state and caller draws; it must not commit another event or advance RNG.

Native v4 getters were requested to expose work function and outward normal alongside layer state, and separate accepted/trial counts for all copied arrays. The native owner confirmed these additions. Reservoir output is five entries per patch; cumulative output is 23 entries (19 sealed source terms, mechanical exchange, particle charge, electrode charge, and electron count).

## Explicit physical limitations

Water support impulse is minus actual weighted m*u. The sealed scenario assumes a fixed infinite-mass support, whose recoil work is zero. A moving finite foil would require relative-frame emission, momentum transfer, and actual mechanical support work. The current frozen source surface must be described as such; no recoil kinetic heat may be invented to fill the missing model.

Collection of a retained PICHandedOff electron without a parent transport trajectory is a manufactured boundary-event check, even if its incoming kinetic energy matches the retained state. It is distinct from the source layer's directly integrated magnetic recollection. Whole PIC/current transport accuracy is not established by a collection-port test.

The sealed tests use an authored electron packet weight of 1e9 to make thermal changes resolvable. That is a component interface fixture, not a predicted source event count or measured current. Any inherited .3 magnetic criterion is also a manufactured override; nominal selection remains .2 unless separately admitted. No authored heat may stand in for an actual source exchange in the v4 integration evidence.

The new v4 factory instead declares a component event weight of 1000, also explicitly not a predicted event count. Electron and water hypotheses occupy disjoint halves of the LowY face. Its prescribed field is a component input, and no field solve is claimed.

## Planned independent checks

1. Reconstruct the joint balance from physical native states and actual packet/path outputs, separately from the implementation's reported residual.
2. Verify each nodal ledger identity and the electron/water/mechanical sums; check weighted q E dot dx against Wfield for every path.
3. Check signed electron energy/reference/image changes against current layer state, and water packet reference zero, kinetic energy from proper velocity, and stable adsorbed reference change.
4. Verify accepted-start source values and identical source candidates across mechanics corrector retries; no duplicate birth debit or source count.
5. Reject unsupported cold reviewed-water state, missing source admissions, stale field/material identities, and source/body mismatch without changing any accepted owner.
6. Restart heat, layers, gas, current paths, inventory, and mechanics together, then compare the next actual trial; rejected staged restart must preserve the old state.

## Review-driven changes and test handoff

The first v4 draft reconstructed electron and water nodal channels by subtracting a large mechanical/combined load. This can erase tiny source terms through cancellation. The native owner changed the implementation to deposit each channel directly into a separate vector; the combined source may still round within the solver's declared bound, while the channel diagnostics retain the actual source operands.

The native owner also confirmed initial and candidate water-range admission, and strengthened restart checks for gas kinetic energy/proper velocity, velocity/proper-velocity consistency, fixed-support impulse, relativistic conversion, source-layer normal, live electron charge/count, and path work versus the cumulative field-work ledger. A rehashed but inconsistent gas-kinetic-energy checkpoint fixture was added by that owner. These confirmations are implementation status, not evidence that the tests passed.

`tests/test_solid_foil_native_v4_review.py` now contains six actual-worker cases selected by `SS304_NATIVE_V4_EXECUTABLE`: two signed mechanical histories using actual birth, continuation, water and conduction; three missing source-admission cases; and reviewed cold-water rejection versus explicit conditional cold-water admission. Its `audit_trial` reconstructs tensor Uel, reference-mass kinetic energy, analytic stable caloric increments, electron reference/image/kinetic state, qE·dx path work, actual gas proper-velocity energy and support impulse, Decimal hi/lo water conservation and reservoir changes, and all seven native nodal channels. The tests require the central CPU7/nice15 environment and do not start builds. Numerical execution remains with the root test owner.
