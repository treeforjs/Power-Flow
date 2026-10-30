The sealed material, spatial-transfer and expansion-v4 components satisfy their
published bounded contracts. They do **not** establish readiness of a single
deforming, evaporating, source/EM/native C19 driver. This audit identifies five
P1 integration gates and three P2 validation/applicability gaps. These are
cross-component requirements, not newly discovered failures of the sealed tests.
Later source/EM/native work can close a gate only with an explicitly bound common
implementation and evidence; it cannot inherit closure from these component seals.

P1 means a conservation, ownership or admitted-domain requirement for the combined
driver. P2 means missing convergence, property applicability or evidence needed
before physical interpretation. No simulations, builds, component suites or new
coupled tests were run. No source/model/guard was changed. The only numerical work
reconstructed existing receipts using exact Fraction arithmetic and their existing
published bounds:524 checks, CPU7/nice15/thread1,0.168s exclusive lease.

| Component | Verified evidence | What it does not establish |
|---|---|---|
| Material v2, seal `1a1845e0…` |12 source hashes and library matched;17,601 existing checks;300K reference zero, valid negative cold h, phase/latent/inverse identities, conditional vapor closure and current-mass capacity | Compressible EOS, measured hot foil/vapor transport, mechanics above its own range, or actual evolving vapor trajectory |
| Spatial material, source `0ed9e7d0…`, commit `f2d3f93` |10 source hashes matched;639 existing assertions; one exact mass/E stock per cell,2x2 local samples and actual first hits,3x2 conduction after mass loss, local guards, rollback/restart | Production typed metal birth, moving geometry, variable-mass v4 mechanics, or a second compatible Q1 heat store |
| Expansion v4, source commit `85ba2a1`, final evidence `7ff9faa` |106 source hashes and four binaries matched;437 Python,166 native and8,731 unchanged heat-v2 checks; signed work, native conduction, source partition, rollback/restart | Deformed EM solve, moving-wall recoil/work, hot/liquid mechanics, whole native PIC/gas transport or a physical compressible EOS |
| Exact gas-count review, source `b8a9a27e…` |47 targeted cases,1,287 oracle assertions in each of normal/UBSan builds; exact imports/debits/credits/leaves and lossy-export refusal | Repair of a separate source Count2 debit or qualification of the new full-time source/motion driver |

The v4 joint histories use a2µm×1µm×2µm component around700K, prescribed fields and
an explicitly authored event weight1000. That weight is not a predicted current
or event count, and the fixture is not the C19 foil geometry.

**MG-01 — P1: bind one physical mass/reference-volume convention and one thermal store.**

The v4 joint demonstration uses reference density7897.403049272394kg/m³. The
spatial transfer and solid-only pilot use the separately declared nominal8000
scenario. The actual-pilot geometry comparison also selects a newly hashed8000
scenario, but is temperature driven and does not close the native caloric ledger.
These are distinct qualified scenarios; the same material header seal does not
make their masses or density assumptions interchangeable.

For a reference cell V0, the material convention requires
`E=m_current*u_model(T)`, `C=m_current*C_app(T)` and `rho_current=m_current/V_current`.
Deformation with fixed mass gives `V_current=J*V0`, hence `rho_current=rho0/J`.
Energy density per reference volume is `(m_current/V0)*u`; multiplying initial
rho by u after metal loss would restore removed energy. Multiplying by a thermal
density and deformed volume inconsistently would also change mass. Apply the
invariant height exactly once to full-height energies, mass and physical areas.

The material owner exposes this current-mass law, and the spatial owner uses it.
V4 mechanics instead assembles nodal masses from `rho_ref*V0/8`, its geometry
adapter reconstructs cell mass as `rho_ref*V0`, and its native heat callback uses
fixed rho. V4's constructor correctly rejects density mismatch. It has no
per-cell metal-loss binding. Running its accepted Q1 caloric state alongside the
spatial FV energy state would create two solid energy stores.

Required before combined admission: one bound density/stress-free reference and
one accepted cell mass/E owner; a conservative extensive mass/E map into the
chosen thermal/mechanical discretization; common handling of removed material
volume and mechanical inertia. The source's adsorbed water inventory remains
separate from condensed SS304, and must not be silently treated as bulk metal.

Evidence: [v2 inventory law](/home/treefor/projects/runtimes/ss304-thermal-properties-v2-20261029/source/include/ss304_thermal_v2.hpp:74),
[v4 density guard](/home/treefor/projects/runtimes/ss304-thermal-expansion-20261029/source/mhdlab/solid_foil_native_coupling_v4.py:71),
[fixed mechanical masses](/home/treefor/projects/runtimes/ss304-thermal-expansion-20261029/source/mhdlab/solid_foil_mechanics.py:195),
[geometry mass reconstruction](/home/treefor/projects/runtimes/ss304-thermal-expansion-20261029/source/mhdlab/solid_foil_geometry.py:133),
[joint caloric receipt](/home/treefor/Documents/Codex/2026-10-28/we/outputs/ss304-thermal-expansion-20261029/foil-support-v4/qualified-final/expansion-checkpoint/joint.json),
[spatial ownership contract](/home/treefor/Documents/Codex/2026-10-28/we/outputs/ss304-spatial-metal-transfer-20261029/spatial-material-port-v1.txt).

**MG-02 — P1: the admitted mechanical trajectory ends before the hot-material/vapor regime.**

Material q1 supports298.2–1800K; q2 extends conditionally to3000K. V4 mechanics
has the narrower298.2–773.15K envelope, excludes coefficient-knot landing and
crossing, and requires positive effective capacity, nominal proof stress,2%
principal strain and positive sampled Jacobians. The open-segment restriction
also excludes the coefficient-table endpoints. It is not enough to check only
T<773.15K or source patch means; actual thermal nodes/mechanical cells must pass.

In the actual-pilot geometry comparison only59 intervals through12.428157577ns
are admitted. Its final native range is300.569674–370.569113K. The next actual
stage crosses373.15K, so the comparison stops. The conditional evaporation sampler
starts at1800K, with no overlap with v4's active mechanical range. There is no
liquid, plastic, free-surface, recession or vapor mechanics closure in these seals.
The separate64x8/125ps solid-only pilot stops at49.105419ns when a local3000K
trial is unsupported, despite large remaining global thermal headroom.

Required: preserve local and coefficient-segment refusal and define the combined
driver's admitted scope. Neither timestep refinement, a global temperature mean,
the q2 caloric alternative nor numerical residual closure extends mechanics or
the3000K material ceiling. This audit supplies no replacement constitutive model.

Evidence: [knot refusal](/home/treefor/projects/runtimes/ss304-thermal-expansion-20261029/source/mhdlab/thermoelastic_energy.py:159),
[crossing refusal](/home/treefor/projects/runtimes/ss304-thermal-expansion-20261029/source/mhdlab/thermoelastic_energy.py:208),
[actual supported geometry prefix](/home/treefor/Documents/Codex/2026-10-28/we/outputs/ss304-thermal-expansion-20261029/foil-support-v4/actual-pilot-geometry-summary.json),
[solid-only material stop](/home/treefor/Documents/Codex/2026-10-28/we/outputs/ss304-material-feedback-pilot-20261029/HANDOFF.txt).

**MG-03 — P1: deformed geometry, boundary current and mechanical/field work are not one coupled operator.**

The Q1 geometry utility verifies reference/current volumes, mass, oriented surface
area closure and sampled Jacobians. It does not prove global injectivity or check
self-contact. The actual geometry comparison gives minimum gap changes of
−22.541192nm free and−23.284622nm with both ends clamped. Its fixed-voltage
parallel-plate factors are shape diagnostics, not C19 field solutions. C19 is
current driven; a simultaneous independent voltage constraint must not be inferred.

V4 explicitly leaves native heat geometry fixed, uses prescribed fields, performs
no deformed EM solve and keeps source surfaces stationary. Source support impulse
is diagnostic with zero fixed-support work; it is not applied to the moving
finite foil. The geometry comparison includes no magnetic force. Those explicit
component approximations do not close moving-wall momentum/work or Maxwell-stress
coupling. Intrinsic k and sigma are not automatically porosity or gas transport
laws after material loss; their physical geometry, volume and measure must remain
consistent in heat, Joule current, fields and remap.

Required: a common geometry trial with physical patches/areas/normals/wall motion,
species first-hit boundaries, reference/current maps and geometry revision; field
and current/charge conservation on that same geometry; mechanical traction/work
and moving-wall source exchange under one owner. Do not add a parallel-plate
correction, arbitrary recoil heat, duplicate Joule/field work or a new effective
transport law to close the missing terms.

Evidence: [v4 explicit readiness flags](/home/treefor/projects/runtimes/ss304-thermal-expansion-20261029/source/docs/interfaces/solid-native-thermoelastic-port-v4.json),
[stationary-source declaration](/home/treefor/projects/runtimes/ss304-thermal-expansion-20261029/source/mhdlab/solid_foil_native_coupling_v4.py:93),
[sampled geometry and self-contact limit](/home/treefor/projects/runtimes/ss304-thermal-expansion-20261029/source/mhdlab/solid_foil_geometry.py:94),
[actual geometry comparison](/home/treefor/Documents/Codex/2026-10-28/we/outputs/ss304-thermal-expansion-20261029/foil-support-v4/actual-pilot-geometry-summary.json).

**MG-04 — P1: reconcile source inventories, native particles and the new count codec before common publication.**

V4 owns local accepted layer states, paths, water packets and Count2 gas inventory
inside its source/heat participant. The new native arena is intended to own the
actual particle IDs and exact gas counts. Both cannot become independent accepted
populations. Active source-layer particles must be excluded from ordinary push;
handoff/return/collection must refer to the same native ID and remaining interval.

The spatial seal qualifies actual v6.2 first hits but labels birth as authored
genesis fixtures and rejects that scope for production. Its production GasPort is
an external owner-attestation contract, not proof supplied by a scope string.
The legacy fixture Count2 pair parser cannot be relabelled as the exact backend's
`{unit_exponent,positive_carriers}` codec. Use actual positive native leaves for
count/current/mass/E; do not use rounded `view().weight` or truncate extended gas
counts back to two limbs.

The fresh exact reconstruction of all eight v4 history intervals found **zero**
source-debit/packet and gas-Count2/packet discrepancies. This verifies those short
histories only. Count2 cannot represent every heterogeneous future debit; common
preflight must prove the separate source-pool delta equals actual newborn weights
or refuse the trial. Exact gas arithmetic does not repair a lost source debit.
Native absorbed-K receipt counters are throughput, not a second energy sink after
the same collection has credited the solid. Metal packet reference+actual K enters
once; there is no extra mean2kT, latent or field-work debit.
For metal, `m_atom=Mmean/NA`, `R_atom=m_atom*u0[J/kg]` and each actual native leaf
contributes `weight*(K_atom+I_atom+R_atom)` once. In contrast, water packet K and
the electron LayerExchange/current-path operands are already weighted; multiplying
those by packet weight again would duplicate energy, charge or work.

Required: root's source/native audit must bind the production typed birth and
collection receipts, exact source pool delta, one ID/motion owner and once-only
mass/E/momentum/current transfers to the same native trial. This audit does not
assert whether a later full-time implementation has already supplied that proof.

Evidence: [v4 local population storage](/home/treefor/projects/runtimes/ss304-thermal-expansion-20261029/source/native/mfem_thermoelastic_v4/thermoelastic_v4.cpp:14),
[v4 Count2 gas update](/home/treefor/projects/runtimes/ss304-thermal-expansion-20261029/source/native/mfem_thermoelastic_v4/thermoelastic_v4.cpp:112),
[spatial GasPort scope](/home/treefor/projects/runtimes/ss304-spatial-metal-transfer-20261029/source/spatial_material.hpp:56),
[exact-count review](/home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-profile/exact-count-independent-review-20261029/REVIEW.md),
[fresh existing-history count reconstruction](/home/treefor/Documents/Codex/2026-10-28/we/outputs/ss304-material-geometry-readiness-audit-20261029/evidence-audit.json).

**MG-05 — P1: component rollback/restart is not the full common accepted boundary.**

V4 successfully restarts heat/source/gas/path plus mechanical u/v, proof and caller
provenance. Spatial material successfully replays mass/E/support/receipts and
conduction against a restored v6.2 journal. The solid-only pilot has byte-identical
upper/lower/field/row restart evidence for its bounded eight-step prefix. These
are different ownership sets, clocks and schemas, not one interchangeable snapshot.

Required: root's combined boundary must include canonical material mass/E, chosen
thermal store, mechanics/geometry and density/reference identities, native exact
particles/IDs/counts, source-layer transitions, producer watermarks and both RNG
streams, fields/current/work/circuit state and common clock. Prepare/retry must
leave every accepted owner unchanged; common preflight precedes publication and
restart must reject partial stages. Lost acknowledgements require candidate-digest
reconciliation or whole-checkpoint recovery. No participant-only restart or old
schema relabelling is evidence of full-driver rollback or crash consistency.

Evidence: [v4 common checkpoint protocol](/home/treefor/projects/runtimes/ss304-thermal-expansion-20261029/source/mhdlab/solid_foil_native_coupling_v4.py:53),
[spatial checkpoint boundary](/home/treefor/projects/runtimes/ss304-spatial-metal-transfer-20261029/source/spatial_material.hpp:213),
[solid-only restart proof](/home/treefor/Documents/Codex/2026-10-28/we/outputs/ss304-material-feedback-pilot-20261029/restart-proof.json),
[v4 restart scope and acknowledgement limit](/home/treefor/Documents/Codex/2026-10-28/we/outputs/ss304-thermal-expansion-20261029/foil-support-v4/native-independent-review.md).

**MG-06 — P2: mapping and convergence evidence remain component specific.**

V4 conservatively maps native nodal caloric energy to homogenized mechanical cell
temperatures and signed cell work back to Q1 nodes. Its native nodal heat identity
and mechanical cell identity are verified separately. This is not a common
cellwise total-energy PDE. The actual-pilot geometry comparison instead uses Q1
corner-average T and explicitly claims no caloric feedback. Temperature averaging
is not a substitute for an extensive mass/E map, particularly with nonlinear cp,
latent heat or mass loss. The FV spatial scheme and Q1 heat scheme require an
explicit conservative interface if both discretizations participate.

The spatial3x2 time comparison reduces peak difference from0.601359K to0.132388K
against a64-piece reference. The solid-only64x8 interval comparison changes peak
T by up to3.281053K; the32-mesh comparison differs by up to71.418209K at242 common
times. No corner/spatial convergence is established. V4's two four-step component
histories and59-step prescribed-T geometry prefix do not establish coupled
source/thermal/field/mechanical convergence or the evolving physical event rate.
Required evidence is convergence of the actual combined discretizations and
quantities, with local phase/energy margins and unchanged guards retained.

Evidence: [native-to-solid map contract](/home/treefor/projects/runtimes/ss304-thermal-expansion-20261029/source/docs/interfaces/solid-native-thermoelastic-port-v4.json),
[spatial refinement](/home/treefor/Documents/Codex/2026-10-28/we/outputs/ss304-spatial-metal-transfer-20261029/disposition.json),
[solid-only time/mesh comparison](/home/treefor/Documents/Codex/2026-10-28/we/outputs/ss304-material-feedback-pilot-20261029/results.json).

**MG-07 — P2: hot condensed/vapor laws are conditional and are not a plasma/EOS model.**

The material model uses fixed-reference Cp-as-Cv energy with h=u_model=0 at300K;
h(298.2K)=−866.079237344J/kg is valid. Melting latent heat is included once in the
declared1672–1727K phase band, with276320J/kg fusion energy. V4 adds conditional
`Phi−T*Phi_T` elastic internal energy and signed
`Qmechanical=Wstress−DeltaUelastic`, not an extra latent or recoverable-energy
heat source. Positive fixed-strain capacity is a guard, not physical EOS validation.

Above1800K q2 uses conditional liquid cp794.96J/kg/K. Selected liquid k is anchored
to Ho-Chu1800K and follows Kim's shape, with normalization1.59844; at3000K selected
k35.2201W/m/K differs from the raw Kim absolute22.034. Sigma follows an anchored
Lorenz model. Vapor is a congruent neutral monatomic proxy: the selected latent
estimate7.40568MJ/kg differs from the raw pressure-slope implied6.51985MJ/kg.
The closed alternative changes pressure (0.803Pa versus raw4.595Pa at1800K).
Its constant3000K/vacuum/accommodation1 loss/cooling illustration is not an
evolving C19 vapor calculation or proof that vapor cooling prevents a hot stop.

No selective Fe/Cr/Ni inventories, metal ionization/excitation/collision data,
measured accommodation, compressible liquid/gas EOS4273, plasticity or calibrated
actual-foil surface model is supplied. Never reuse H2O collision data for the
proxy metal or apply condensed sigma/k to dilute plasma gas. Required physical
interpretation must retain the exact selected scenario and discrepancies; no
guard extension or replacement property assumption is supplied by this audit.

Evidence: [material handoff and uncertainty](/home/treefor/Documents/Codex/2026-10-28/we/outputs/ss304-thermal-properties-v2-20261029/HANDOFF.md),
[material verification](/home/treefor/Documents/Codex/2026-10-28/we/outputs/ss304-thermal-properties-v2-20261029/verification.json),
[conditional thermoelastic capacity](/home/treefor/projects/runtimes/ss304-thermal-expansion-20261029/source/mhdlab/thermoelastic_energy.py:93).

**MG-08 — P2: actual-pilot geometry is a conditional comparison with incomplete build/run provenance.**

The geometry audit preserves244 actual thermal states and independently verifies
upper/lower checkpoints and MFEM node-order certificates. Completion binary and
material-header hashes match the observed snapshot, but its formal source-to-build/
run seal remains explicitly unverified. The admitted59-interval geometry prefix
uses a nominal8000kg/m³ scenario, Poisson0.29, free/end-clamped support brackets,
and300K stress-free zero u/v at the first C19 time. Actual fixtures, foil density,
surface microphysics and missing prehistory are not reconstructed or calibrated.

It prescribes T, does not project/debit native caloric energy, includes no magnetic
force and gives no production field handoff. It is evidence of consistent geometry
under stated inputs, not a calibrated coupled prediction. The independent
solid-only pilot is reduced z-invariant MQS with all water, emission, plasma,
metal vapor and mechanics off. Its terminal-work identity correctly separates
Joule heat from backward-Euler magnetic damping; damping is not extra heat.

Required before upgrading that evidence: close source/build/run provenance and
bind the actual selected initial/support/density/drive identities. Retain the
conditional status where measurements or physical support are absent; do not
turn shape comparisons or missing prehistory into implicit production inputs.

Evidence: [actual-pilot geometry provenance](/home/treefor/Documents/Codex/2026-10-28/we/outputs/ss304-thermal-expansion-20261029/foil-support-v4/actual-pilot-geometry-summary.json),
[solid-only formulation and inactive physics](/home/treefor/Documents/Codex/2026-10-28/we/outputs/ss304-material-feedback-pilot-20261029/component-contract.json).

The existing local and global bookkeeping is verified within its declared scope.
For v4, only certified electron Qwall enters solid heat; electron K/reference/
image and field work remain distinct. Water gas reference is zero, its negative
adsorbed reference rises on removal, and actual packet K/internal/adsorbed changes
are charged once. In the existing v4 runs Joule and authored heat are zero. The
fresh reconstruction obtains maximum joint residual4.407121087e−21J and maximum
native nodal residual1.353713267e−23J, within unchanged owner bounds. Two nonzero-
displacement, zero-duration quantum bridges are retained, with weighted qE·dx
work checked directly. These facts do not supply field-solver/circuit conservation:
when EM becomes an accepted owner, its work/Joule/boundary ports must cancel once
in the combined ledger. Prescribed-field Wfield cannot also remain an independent
external input for that same internal transfer.

All source hashes, evidence paths, reconstruction results and process identity
are in [evidence-audit.json](/home/treefor/Documents/Codex/2026-10-28/we/outputs/ss304-material-geometry-readiness-audit-20261029/evidence-audit.json).
Production PID278902/start5288339 remains T/frozen1, original executable SHA
`cddd9f1867b50f471b6e8ce49e698c45414b568acb6e42c5a4364d5c5fb37797`.
No integration test, GPU work, production signal, deployment, model assumption or
guard change was made.
