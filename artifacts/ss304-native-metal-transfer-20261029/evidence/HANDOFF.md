# Native thermal-metal-vapor transfer qualification

The isolated package is `/home/treefor/projects/runtimes/ss304-native-metal-transfer-20261029/source`. Entry header: `metal_transfer.hpp`. The implementation forks only the native arena surface API/journal support, and adds the immutable adapter plus exact binary arithmetic. The original shared native code, native heat header, and sealed v2 material release remain unchanged. `arena-surface-extension.patch` is the reviewable arena diff. `thermal-metal-vapor-port-v1.json` gives exact API and ownership rules.

## Qualified behavior

**75 native checks,190 existing arena checks and171 independent highprecision assertions passed.** Actual CPU MFEM, actual native GapPic carriers, native drift/first-hit collection and the existing PicArena are used. There are no GPU, production, dependency or shared-checkout changes.

One canonical compensated solid-atom Count supplies the actual material mass. One exact dyadic total-energy record and three exact support-impulse components survive sub-ULP transfers and restart. Gas packets/carriers belong exclusively to PicArena/GapPic after publication. The component keeps no duplicate gas population or second solid inventory.

The native heat Coupler is an ephemeral immutable trial evaluator reconstructed from that canonical stock. The qualified thermal fixture is one uniform1x1 block. Current mass determines capacity/reference-volume energy. The load is **E_next-m_next*u(T_old)**, so removing material carries its old sensible/phase energy once. The cache's authored_heat field is this numerical load representation, not an additional physical energy reservoir.

For every actual outgoing or incoming carrier, the transfer uses native proper-velocity kinetic energy plus registered internal and reference energy. Five caller-owned open-unit draws can generate flux-weighted half-Maxwellian proposals through the existing heat API's native kinematic conversion. Actual sampledK is retained; the mean2kT is not added or substituted. The particle registry binds the source's atomic reference energy. Fusion/vaporization and field work are not independently debited again.

Fixed support receives the exact opposite3D proper momentum. Its work is0 because it is explicitly immobile. Moving recoil/support work is rejected for another admitted mechanical closure; the separate773.15K mechanical-data guard is unchanged.

## Real native collection and atomicity

`prepare_native_collection` consumes a typed staged `PicArena::prepare_transport` token and its **actual first-hit species sink receipts**. It credits incomingK/internal/reference and material count, then publishes the original native transport token with the material snapshot. It never recreates a removed particle or debits gas again. Consumed sink IDs persist and are replay checked. In this bounded one-block primitive all new sinks must belong to this upper-back-face metal proxy; H2O, unrelated sinks and escapes are refused for parent multi-owner routing.

Native absorbed-energy diagnostics remain throughput receipts. When the root credits a consumed sink to heat it must not also count itsK as an external boundary-energy loss. Exact represented per-atom kinetic currency can differ from aggregate floating diagnostic sums by primitive rounding; the independent160-digit kinetic oracle bounds this explicitly. Such a difference is numerical accounting, not another physical heat source.

Prepare leaves both accepted owners unchanged. Publication uses the arena's nonallocating native swap followed by the canonical shared-snapshot publication. Root owns whole-coordinator preflight and exposure. Later failures require reverse-order rollback. Faults at all native staging points, after both stages, and after gas publication restore IDs, particles, fields, clocks and material. Generation guards prevent stale-token resurrection.

Restart binds compiled native-source identity, sealed material identity, geometry/axial frame, source configuration and actual arena prefix. The compact material journal replays mass, actual nativeK, internal/reference energy and support momentum against accepted arena stages. Rehashed changes to energy, impulse or consumed-sink identity fail before acceptance.

## Bounds and evidence

The sealed material remains298.2-3000K; emission requires1800-3000K. Negative reference energy below300K is retained. Actual MFEM checks cover heating/cooling through fusion and cooling300->299K. A single-atom debit below inventory and total-energy ULPs is preserved. Heat that would cross an exact qualified endpoint is rejected even when adding it leaves rounded total energy unchanged. Unrepresentable weighted particle mass/energy and nonfinite native heat activity are refused without a source floor.

The conditional surrogate is `SS304_CONGRUENT_NEUTRAL`, with no H/O stoichiometry or inherited water collision data. Any registered collision channel touching it is refused. Physical internal energy is the ground-state0 of this scenario; nonzero internal-energy integration tests are explicitly manufactured. No metal chemistry, excitation, ionization or full alloy EOS is added.

`transfer/native-tests.json` contains native audit receipts. `existing/arena-tests.json` preserves the190 water/arena regressions. `highprecision-audit.json` independently uses exact Python Fractions for zero represented-binary mass/energy/momentum discrepancy and160-digit Decimal kinetics with an8ULP arithmetic limit. The observed worst fraction is0.06401 of that kinetic limit. Numeric checks do not change solver/scientific guards or validate actual foil properties.

## Root integration limits

Use the supplied actual-native xy geometry and explicit axial center; both centered and0..25mm z frames are tested, and the actual native upper foil/depth/stage must match. Root owns HK rate selection, sample timing, all-source event ownership, EM heat admission, global clocks, checkpoint manifest and multi-owner sink distribution.

This package qualifies a uniform-block transfer primitive. The full spatial heat owner still needs per-control-volume canonical mass/energy and a conservative source/deposition/remap mapping. It must adopt a single canonical stock rather than add a second accepted energy inventory alongside this one. This package does not qualify a nonuniform full-foil solve or60ns coupled history. Geometry owner must handle the volume associated with removed/condensed mass and mechanical work. Production remains manually stopped and cgroup frozen.
