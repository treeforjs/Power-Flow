# SS304 thermal/vapor v2 handoff

The new isolated native source is `/home/treefor/projects/runtimes/ss304-thermal-properties-v2-20261029/source`; header `include/ss304_thermal_v2.hpp`, library `build/libss304v2.so`. The sealed1800K release is unchanged. Source/evidence bindings are in `source-manifest.json`; callbacks, ownership and restart fields are in `thermal-material-port-v2.json`.

## Ready low-bound correction

Select **q1 nominal_1800** for298.2-1800K. It preserves h=u_model=0 at300K and is bitwise identical to the old release on300-1800K (1507 compared states). Mills' explicit298.2K table row and Ho-Chu's293/300K transport rows support the new reviewed cold interval. At298.2K,h=-866.079237344J/kg. Cold patches can cool below300K without clamping or raising initialT. Negative reference enthalpy is valid. This needs the new adapter/revision; do not just relax the old callback guard.

## Conditional hot alternative

Select **q2 source_conditioned_liquid_3000** for298.2-3000K. Above1800K, Kim's estimated304L liquidCp794.96J/kg/K extends energy continuously. Liquidk uses Kim's source shape normalized to admittedHoChu1800K conductivity; sigma follows the anchored constant-Lorenz model. These are explicit scenario assumptions. Normalization is1.59844: at3000K selectedk35.2201W/m/K versus Kim absolute22.034. Report/test this discrepancy rather than presenting it as measured304.

The1800-3000K vapor alternative matches Kim's7.40568MJ/kg latent estimate and pressure at3000K, uses an ideal congruent monatomic vapor, and enforces Kirchhoff and Clausius-Clapeyron identities. The source pressure slope alone implies6.51985MJ/kg; both source and model pressure curves are reported. The alternative pressure is0.803Pa at1800K versus the rawsource4.595Pa. This is a thermodynamically closed conditional alternative, not the missing fullEOS4273.

At3000K, selected condensedu is2.153727855MJ/kg, sigma480516S/m. Constant3000K/vacuum/accommodation1 for60ns gives2.5121mg/m2 metal loss and18.04J/m2 cooling in this model. The corresponding metal column is2.735e19atoms/m2; each extra1eV/atom costs4.382J/m2. `quantified-material-budgets.json` contains temperature points, cold/hot margins, pressurework, active-layer mass fractions, gas transport assumptions and source discrepancies. These are not a prediction of the evolving temperature trace.

## Concrete ownership interface

Heat owner holds one canonical metalmass/energy inventory, separate from adsorbedH2O. Gas owns actual weighted emitted/incident packets and caller-owned sampling. `prepare_transfer` prepares a reversible transaction using actual signeddm and actual packetK:

`E_packet=dm*u0_atom+K_actual`, `m_next=m-dm`, `E_next=E+Q_other-E_packet`, `T_next=inverse(E_next/m_next)`.

No extra mean2RT or latent-vaporization debit is added. Parent preflights/commits heat,gas,geometry and producer watermarks together. Replay/stale revision, unsupportedenergy or depletion leaves the accepted base/output untouched. Pure signed condensation is also checked. After mass loss use `(m_current/V_ref)*u` and `(m_current/V_ref)*C_app`, exposed by `evaluate_inventory`; retaining initialrho_ref*u would create energy. Intrinsick/sigma is not an effective-medium law. Geometry owner receives removed referencevolume dm/rho_ref and owns recession/remap/work. The independent773.15K mechanics guard is unchanged.

A supplied positive collisiondiameter parameter enables a strictly neutral hard-sphere gas thermal coefficient; no measured diameter or ionized-gas conductivity is invented. Liquid transport is never applied to dilute/plasma gas. Congruent metal vapor is an approximation: Kim's pure-element pressures predict a Cr-enriched vapor. Real selectiveFe/Cr/Ni inventories, rates and charge/excitation energies belong to the gas/plasma owner. An explicitly open escapedmetal reservoir can balance mass/energy as a simpler bounded alternative; do not invent a SS304 molecule in the H2O network.

## Evidence and next actual run

**17,601 assertions passed**: cold references,1507 bitwise old/new states,2507 broad-range inverses, derivative/phase identities, source-conditioned transport, signed flux, energy/mass closure, rollback/replay, inventory-weighted heat capacity and unsupportedtrial energy budgets. Numeric limits test arithmetic, not property uncertainty. CPU7/nice15/one-thread lease; noGPU/MFEM/PIC solve or dependency changes. Production wasT/frozen1 with the original identity before/after.

Parent should first reevaluate the actual thermal/EM trajectory with sigma(T) and the new caloric law, retaining interval/refinement diagnostics. `assess_energy(m,E_trial)` returns signed cold/hot margins; an above3000K trial reports its exact unsupportedexcessJ and aborts, retaining the heat. The old constant-property corner trajectory is not a forecast for this model. SourceEOS4273 remains unavailable as a numericgrid: the manuscript directs requests to SESAME; no request was sent. Compressible expansion, critical/high-density gas, plasma/WDM transport and actual-foil calibration are not supplied by this alternative.
