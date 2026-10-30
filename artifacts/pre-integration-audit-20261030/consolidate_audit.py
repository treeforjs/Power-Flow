#!/usr/bin/env python3
"""Preserve the pre-integration audit; performs no numerical solve or launch."""
from pathlib import Path
import datetime
import hashlib
import json
import shutil

O = Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs')
D = O / 'pre-integration-audit-20261029'

# Each disposition separates an implementation defect from an unimplemented
# capability or an explicitly conditional scientific model.
rows = [
 ('A01', 'defect', 'P1', 'Legacy entry point can invoke native Poisson in full-time mode',
  'Empty legacy compound preparation reaches stage/solve_fields and can publish without updating the full-time boundary. The intended new coordinator avoids the route, but the public mode boundary does not.',
  'Refuse every legacy entry point before callbacks/mutation; retain zero-call traps for the forbidden solver and legacy advance.', 'native;contract'),
 ('A02', 'defect', 'P1', 'Proposed field, source and geometry can change after their seals are captured',
  'Final common preflight does not compare mutable proposed field/heat to immutable native final seals. Narrow source probes also accept changed particle weight, spectrum, reused owner address and mutated mesh with old metadata.',
  'Use opaque immutable candidates, live owner/generation and geometry identities, and recompute final native/field/heat/source/RNG/time bindings before allocation-free publication. Mutation/fault tests must leave every accepted owner unchanged.', 'contract;source'),
 ('A03', 'defect', 'P1', 'Private scope identity, lifetime and reentry are incomplete',
  'Intermediate001 shares an outer identity between sibling scopes and retains a raw owner pointer. Factories do not reserve exposure slots; per-scope builder flags allow sibling reentry. Exposure receipt allocation can throw after partial private publication. Candidate002 fixes only part of this and is not substituted into the review.',
  'One outer reservation/build authority; exact scope keys and weak owner lifetime; prepare all allocating receipts before private swaps. Test destroyed-owner, sibling, reentry, stale and throwing paths.', 'native;contract'),
 ('A04', 'defect', 'P1', 'Native stochastic and closed-source seals are weaker than the advertised boundary',
  'Native admission does not independently bind first RNG to accepted RNG, both stages to one owner, and final RNG to actual second completion. Empty source receipts permit a changed source hash. Motion JSON does not fully validate schema and unique complete path/endpoints. Root supplies stronger normal-path checks, which do not repair the primitive.',
  'Bind the complete RNG chain and unchanged closed-source seal natively; admit only a typed genuine motion receipt with complete census, geometry, path, time and policy checks.', 'native;contract'),
 ('A05', 'defect', 'P1', 'A nonempty full-time checkpoint cannot save or restore',
  'Both routes call the legacy journal validator, which rejects full_time_outer_begin and other full-time entries before semantic replay. Merely extending its allowed names leaves incorrect step/previous-dt accounting. Fresh-genesis restart evidence did not cover this.',
  'Implement the exact bounded full-time grammar and accepted-step accounting. After component fixes, require genuine one-step fresh-process restart, split/full parity and corruption/refusal evidence.', 'native;contract'),
 ('A06', 'defect', 'P1', 'Restart and genesis admission need checked types and bounds',
  'Raw option bytes are copied into a bool-bearing C++ object; malformed representations can cause undefined behavior. JSON integer conversions can truncate/wrap before replay. Collection caps are incomplete before copying. Actual option operands are not independently bound in the native prefix; authored aggregate charge/energy can overflow.',
  'Use an explicit typed option codec, exact nonnegative integer/range checks, preallocation shape/count/depth limits, direct native-option commitments and finite aggregate ledger checks. Preserve semantic replay and existing caps.', 'native;contract'),
 ('A07', 'defect', 'P1', 'Source reservoir debit and packet allocation are not exact',
  'Actual component preparation emits about 1e-30 particles from {1e15,-0.01} with unchanged inventory limbs. A separate unit release into ten packets exceeds its exact debit by 5/2^55. Exact gas arithmetic preserves the upstream excess rather than correcting it.',
  'Use the existing exact-count representation across the complete source debit/newborn credit transaction, or explicitly refuse unrepresentable debits. Exact-sum all actual returned leaves; do not introduce another count backend.', 'source;material'),
 ('A08', 'defect_or_admission_gap', 'P1', 'Tiny positive source contributions can silently lose charge or heat',
  'A roughly 1e-310 positive-weight birth has rounded-zero charge and wall heat. Separately, an exposure-product scalar underflows despite a finite log expected count and current_underflow=false. V3 documents the latter; the external count owner must honor it.',
  'Preserve the required weighted quantities in an admitted representation or refuse with explicit unresolved scope. Use finite-log exposure admission; never floor, resample, or silently treat positive flux as zero.', 'source'),
 ('A09', 'defect_or_missing_guard', 'P1', 'Initial absolute field constraints are not enforced',
  'Incremental Gauss checks can preserve an inconsistent initial charge/field state. A narrow actual Fields::prepare probe accepts and preserves divB about 50 T/m while energy and Gauss-increment checks pass. This is manufactured evidence, not a measurement of the paused run.',
  'Admit initial total E against actual canonical charge and conductor/boundary ownership; admit initial divB and check preservation. Reject inconsistent states without inserting a second self-field solver.', 'field;contract'),
 ('A10', 'defect', 'P2', 'Particle/source midpoint differs from the compensated Maxwell midpoint',
  'The joint port averages high field words while the Maxwell receipt includes compensation. The operator probe found 9,354 differing E edges, maximum 1.4552e-11 V/m. Small magnitude does not establish the exact shared-field contract.',
  'Publish one authoritative total midpoint from the field owner for particles, sources and work. Test nonzero compensation; do not loosen the work guard.', 'field;contract'),
 ('A11', 'unimplemented', 'P1', 'There is no whole-owner atomic restart or long-run history scheme',
  'The native, parent, field, thermal, source/count/RNG and future forcing configurations do not yet form one atomic disk snapshot. Five journal rows per empty step and a 4096-row ceiling allow at most 819 steps, 8.19 ps at 10 fs, before earlier event/byte/headroom limits. About 5.93 million 10 fs steps are needed from the supported start to 60 ns.',
  'Bind exact drive/configuration and all accepted owners in one manifest, restore privately, and use validated bounded history anchors plus external receipts. Do not just raise caps or drop history. Bound transient forks and memory before allocation.', 'native;contract;material;environment'),
 ('A12', 'unimplemented', 'P1', 'The staged common driver is resident-only',
  'It explicitly excludes source births, active microscopic layers, collection and escape. Source arrival/count/RNG/ID ownership and newborn collision age are external. No full joint step or physical source-coupled trajectory has been accepted.',
  'Implement one canonical population/source-layer mapping, finite inventories, actual physical birth times and remaining flight, first-hit collections/escapes, current and exactly-once energy exchange. Qualify a stated arrival/weight model and active-layer collision scope.', 'source;contract;field'),
 ('A13', 'scientific_coverage_gap', 'P1', 'The water mechanism is a conditional subset, not full breakdown',
  'Generic v6 omits low-energy elastic/angular and rovibrational transport, several attachment/ionization branches and fragment/ion continuation. New v7 metadata maps only part of the positive 2a1 hole and is not the integrated bridge; its generic handoff also mixes old and current metadata. Product conservation does not validate imposed secondary spectra.',
  'Version and qualify supported exclusive channels and product distributions; implement native channel-set admission and physical collision-site/current treatment for capture. Assess omissions on actual densities/EEDF without inventing rates or renormalizing positive missing branches.', 'chemistry;water-data'),
 ('A14', 'model_and_resource_gap', 'P1', 'Broad-population chemistry has a model seam and evidence ceiling',
  'The selected cross-section sum jumps 3.0329-fold across adjacent binary64 energies at 1 keV without overlapping double counts; this is not a measured ensemble-rate or EEDF jump. The evidence collector refuses 257 energies by default and 1025 at its maximum; a captured 46,726-carrier pre-kick fixture evaluated against stationary H2O has 46,726 distinct relative energies.',
  'Resolve supported coverage and quantify the model seam rather than smoothing it arbitrarily. Add bounded streaming evidence with equivalence checks, and converge macroparticle noise/weights independently of exact count arithmetic.', 'chemistry'),
 ('A15', 'conditional_physics_and_missing_integration', 'P1', 'Thermal and field-assisted emission remain essential but not fully admitted',
  'The unified supply/barrier law includes both mechanisms. Historical extracting probes have unsupported outgoing phases; tangent and oblique diagnostics do not qualify a general replacement. Unit return accommodation and clean 4.3 eV stationary-surface assumptions are explicit closures, not demonstrated actual heated layered-electrode properties.',
  'Retain both mechanisms and one source integral. Evaluate actual one-sided total midpoint E/B and accepted surface material; bound or qualify unsupported phases and return behavior, and test held-T time error. Do not add independent thermionic/FN populations or impose shielding.', 'emission;source'),
 ('A16', 'unimplemented_and_conditional_physics', 'P1', 'Material, expansion and field geometry do not yet share one mass/energy owner',
  'Qualified scenarios use distinct reference densities and thermal stores. The reviewed Hall thermal port uses the separate evaluated-solid 298.2–1500 K scope, not the material q2 extension to 3000 K. Expansion has fixed-mass/reference-geometry restrictions and no deformed EM solve. Its actual prescribed-T comparison ends at 12.428 ns; the separate solid thermal pilot refuses near 49.105 ns at 3000 K. These are component limits, not a predicted failure time of the future coupled model.',
  'Choose consistent declared density/reference conventions for each material and phase, with one extensive mass/energy owner; conserve geometry/current/wall work and removed mass. Preserve local material/mechanical guards. The newly requested layered, melt-triggered energy-limited release is a separate closure to implement and qualify, not a silent extension of old seals. Released metal/oxide needs explicit species or a removed-material reservoir with reference energy; the current water catalog cannot supply those kinetics.', 'material;surface-request'),
 ('A17', 'data_and_validation_gap', 'P2', 'Drive, RGA and observation assumptions remain conditional',
  'The new selected driver uses C19 current recovered from published average-display vertices, with unknown trigger registration and terminal/AK voltage as outputs; old recalled 220 kV feed/110 kV midpoint are conditional comparisons. Earlier B comparisons differ about 12–14% and do not validate the new 3D field. The published/vector-figure RGA recovery has 50 mass channels, not raw acquisition or calibrated neutral fractions. Numerical density slabs are not calibrated slit response.',
  'Bind current bytes/clock/topology, compare field and terminal measurements without independent conflicting current/voltage constraints, retain RGA calibration/composition uncertainty, and keep actual 3D collision volume separate from cm^-3 display/slit sampling.', 'root-data;field;chemistry;water-data'),
 ('A18', 'environment_admission_gap', 'P1', 'Pinned userspace does not freeze every host dependency or guard',
  'Selected frozen userspace, native and CUDA hashes match. The host kernel release is guarded, but GPU module/build identity is observed rather than admitted. Python assert guards can be disabled by inherited optimization. CPU leases/thread settings require an external wrapper; compiler/MPI overrides and complete launch seals need explicit admission.',
  'Use non-assert checks, a sanitized recorded environment, explicit CPU/GPU leases, complete source/build/runtime pins and expected kernel/module/device identity. A changed host must refuse pending requalification, not silently refresh pins.', 'environment'),
 ('A19', 'validation_gap', 'P2', 'No converged continuous coupled plasma result exists',
  'Legacy parent 442 checks still pass on intermediate001 with no physical trajectory advancement. Large component assertion totals, the paused 26.4787 ns electron-only run, and solid-only trials do not establish plasma shielding, breakdown, sheath/grid/particle convergence or a full 60 ns history.',
  'After defect closure, progress through component adversarial checks, an exact-version common-owner short test, physical pilot/refinement, then bounded production qualification. Record actual accepted range, limitations, ledgers and white cm^-3 figures. No long-run launch yet.', 'all'),
]

findings = [{'id':a,'kind':b,'priority':c,'title':d,'evidence':e,'closure':f,'reviews':g.split(';'),'status':'open'} for a,b,c,d,e,f,g in rows]
(D/'findings.json').write_text(json.dumps({'schema':'consolidated-pre-integration-audit-v1','disposition':'HOLD_FULL_INTEGRATION','findings':findings},indent=2)+'\n')

intro = '''PRE-INTEGRATION AUDIT — 2026-10-29 work

Disposition: HOLD full integration tests and production launch. The requested
audit found concrete correctness defects, unimplemented common-driver functions,
and scientific assumptions that must remain explicit. No full joint test ran.
The original electron-only job remains stopped/frozen and untouched.

This review covers the exact intermediate001 native source, staged root
coordinator, sealed Hall/source/material/chemistry components, launch contracts,
and existing scientific evidence. It does not claim an exhaustive proof or
qualify any physical result. Narrow component counterexamples and read-only
source inspection are identified separately in the individual reports.

The reviews below preserve their original identifiers and priorities. This
consolidation groups overlapping findings; 19 groups are not 19 independent
runtime failures. A declared approximation is not relabeled as a code defect.
No current source seal was silently replaced with a partially corrected draft.

The latest request adds a reduced energy-limited metal-release model and a
layered initial surface. It does not remove the numerical defects or authorize
unlimited energy, invented surface composition, or immediate full ionization.
See surface-model-request.txt and the separate material appendix when present.

'''
lines=[intro]
for r in findings:
 lines += [f"{r['id']} [{r['priority']}; {r['kind']}] {r['title']}\n",r['evidence']+'\n','Closure: '+r['closure']+'\n','Evidence reviews: '+', '.join(r['reviews'])+'\n\n']
lines += ['''CORRECTIVE ORDER AND RELEASE GATES

1. Preserve audited inputs and reports. Fix deterministic mode, token, final
   seal, exact-count, field-constraint, checkpoint codec and launcher defects in
   separate versioned sources. Keep old seals and the paused run intact.
2. Require focused reproducer closure, legacy parity, failure/rollback cases,
   sanitizer coverage where relevant, and read-only independent review. Full
   integration stays held while these prerequisites remain open.
3. Complete common source/collection/channel-set ownership and full restart;
   qualify only the explicitly supported physics and source/material scope.
   A closed parent fixture must never be advertised as water breakdown.
4. Admit a short physical coupled pilot only after the reviewed exact versions
   have passed the preceding gates. Refine time, mesh, particles, source model
   and exterior domain; audit mass, H/O, charge, momentum, energy and circuit
   transfer. Diagnose screening from these quantities, not cloud appearance.
5. Before a long run, demonstrate bounded history/storage and restart cost and
   a measured performance envelope. Confirm the supported material/chemistry
   range can cover the planned interval; do not bypass a guard for 60 ns.

Historical evidence is still useful: exact gas-count arithmetic, the explicit
two-half-exposure clock, component work/heat conservation and the 442 legacy
parent controls remain valid within their stated scope. None repairs the
confirmed common-driver issues by implication.

Q1 gather cache: the explicit invalidate() correction is already committed in
raikou-mfem a7e44f7e9c3f7bd7d653d7e5a31cf3d850ee8e7d. The current production code
does not call the unversioned renumbering operation. No duplicate patch or
replacement of the frozen paused executable was made during this audit.

root-inputs.json pins selected root inputs. Each evidence review carries its
own exact source line/hash manifest and scope. preserved-production.json records
the verified stopped process identity; it is not permission to signal it.
''']
(D/'AUDIT.txt').write_text(''.join(lines))
print('Wrote 19 grouped findings and audit draft; final evidence index is a separate step.')
