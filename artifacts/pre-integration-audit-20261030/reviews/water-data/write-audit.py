import hashlib
import json
import shutil
import subprocess
from pathlib import Path

R = Path('/home/treefor/projects/runtimes/h2o-collision-data-20261029')
W = Path('/home/treefor/sim-workspace/powerflow/outputs/h2o-collision-data-20261029')
A = Path(__file__).resolve().parent
O = W.parent

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def load(p):
    return json.loads(Path(p).read_text())

def write(name, value):
    (A / name).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')

def ref(name, start, end=None):
    return {'path': str(R / name), 'start': start, 'end': end or start}

water = load(R / 'data/collisions/h2o-gas-20261029-deep-valence.json')
v7 = load(W / 'coordination-handoff-v7.json')
h = load(W / 'fragment-attempt-001/validation.json')
snapshot = load(A / 'v8-working-snapshot.json')
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=R, text=True).strip() == v7['remote_commit']
for name, want in snapshot['files_sha256'].items():
    assert sha(R / name) == want, name
assert sha(v7['source_archive']) == v7['source_archive_sha256']
assert water['full_fragmentation_ready'] is False

write('observable-catalog.json', [dict(
    **{k: c.get(k) for k in ['id', 'observable', 'model', 'support_eV',
        'zero_below_eV', 'unit_SI', 'ordinate_scale_to_SI', 'source', 'source_locator',
        'use', 'exclusive_products', 'collision_process_admission', 'interpolation', 'lineage']},
    table_points=len(c.get('energy_eV', []))) for c in water['channels']])

findings = [
    {'id': 'D01', 'severity': 'medium', 'class': 'confirmed_metadata_inconsistency',
     'finding': 'Generic api-handoff.json combines a v7 dataset/current_commit/current_native_checks with unprefixed v1 commit/native_checks/archive fields; outer3_coverage_warning still says Auger missing. Historical outer3 model exclusions remain correct for that subset, but the generic handoff is not a single-version import manifest.',
     'evidence': [{'path': str(W / 'api-handoff.json'), 'lines': [4, 19, 28, 94, 142, 147]},
                  {'path': str(W / 'coordination-handoff-v7.json')}],
     'impact': 'An importer choosing unprefixed fields can bind the wrong commit, tests or archive to the current dataset.',
     'remedy': 'Use the immutable coordination-handoff-v7.json and its exact contract/archive hashes for consolidation. A normalized generic handoff should be a separately reviewed metadata revision; none was silently edited in this audit.'},
    {'id': 'B01', 'severity': 'blocking_for_full_breakdown', 'class': 'explicit_missing_physics',
     'finding': 'Only 27% of the modeled 2a1 primary-vacancy rate has exclusive products; prompt Hplus58% and Oplus15% stay positive and unmapped. Hplus does not determine OH versus H+O neutral companions; Oplus does not determine H2 versus two H. Oxygen1a1 core and direct projectile double/triple channels remain unqualified.',
     'evidence': [ref('include/raikou/h2o_deep_valence.h', 7, 18), ref('include/raikou/h2o_deep_valence.h', 90, 95), {'path': str(R / 'data/collisions/h2o-missing-rate-ledger-v7.json')}],
     'remedy': 'Keep the two positive remainders and all unknown rates explicit. Source-qualified absolute exclusive event rates, neutral completion and states are required. Never renormalize the mapped27% or count the old full2a1 diagnostic a second time.'},
    {'id': 'B02', 'severity': 'blocking_for_validated_full_physics', 'class': 'declared_conditional_model_limits',
     'finding': 'Exact event conservation does not validate equal-energy primary/secondary sharing, ground states, isotropic breakup, fixed hole lines or a collapsed isolated-site cascade. The high-velocity matrix inferred using proton data is applied conditionally to electron BEB holes, not measured over every1..10keV electron energy. No proton collision rate was imported.',
     'evidence': [ref('include/raikou/h2o_outer_valence.h', 44, 74), ref('include/raikou/h2o_deep_valence.h', 35, 37), ref('include/raikou/h2o_deep_valence.h', 98, 133), {'path': str(R / 'data/collisions/h2o-collision-port-v7.json')}],
     'important_budget_eV': {'OHplus_Hplus_and_Auger': 5.16456, 'Oplus_Hplus_H_and_Auger': 0.163505},
     'remedy': 'Require explicit model opt-in and declare no field work/quenching/lifetime resolution at the collision site. Calibrated spectral/internal/angular data need a new model revision; especially do not describe the0.163505eV endpoint budget as measured KER/Auger spectroscopy.'},
    {'id': 'B03', 'severity': 'blocking_for_full_breakdown', 'class': 'unimplemented_fragment_continuation',
     'finding': 'Published v7 contains water-target kernels, not continuation for all produced fragments. The separate uncommitted v8 candidate adds only ground-protium H1s ionization,14..10000eV. OH, O, H2, molecular/atomic ion and anion excitation, ionization, detachment, recombination, charge exchange and heavy-particle families are not supplied by that addition.',
     'evidence': [ref('include/raikou/hydrogen_fragment_ionization.h', 10, 27), ref('include/raikou/hydrogen_fragment_ionization.h', 52, 67), {'path': str(A / 'v8-working-snapshot.json')}, {'path': str(W / 'fragment-continuation-research.json')}],
     'remedy': 'Preserve H as a separate target/density with oneH debit, zeroH2O debit and net electron gain1. Do not deploy this candidate during the audit. H2Oplus dissociative recombination is only a primary-source research lead; no rate/kernel was admitted. A thermal300K recombination coefficient cannot replace energy-resolved nonthermalPIC data.'},
    {'id': 'B04', 'severity': 'blocking_for_complete_transport', 'class': 'explicit_angular_and_internal_coverage_gaps',
     'finding': 'Full-angle event transport is conditional free-atom IAM only at1..10keV. Molecular low-energy integral/momentum data do not complete angles or rotational energy transfer. SongTable2 covers2..100eV with10..130degree angles and unresolved rotation; no full numerical angular extension was admitted. Other attachment resonances, Ominus products and most excitation/dissociation state outcomes remain unavailable.',
     'evidence': [ref('include/raikou/h2o_elastic_angular.h', 40, 73), ref('include/raikou/collision_data.h', 150, 178), {'path': str(R / 'data/collisions/h2o-missing-rate-ledger-v6.json')}],
     'remedy': 'Retain angular/energy gaps. IAM is incoherentO+2H, not measured molecularDCS; source atomic accuracy is not a molecular uncertainty bound. The1.9057882329984954% mesh-integral discrepancy against rounded atomic totals is only a numerical binding diagnostic.'},
    {'id': 'B05', 'severity': 'blocking_for_measured_background_initialization', 'class': 'unidentified_experimental_input',
     'finding': 'The corrected figure recovery has50 nominal mass channels, not450 independent measurements. RGA channels include ionizer cracking and species overlap; peak16/17/18 does not establish neutralO/OH/H2O fractions. Instrument settings, sensitivity/cracking matrix, raw timing, gas temperature and independent total-gauge calibration were not recovered.',
     'evidence': [{'path': str(O / 'rga-channel-recovery/figure_2_3_rga_recovery.json')}, {'path': str((O / '../work/smith-dissertation-2025.pdf').resolve()), 'PDF_page_one_based': 44, 'printed_page': 19}, {'url': 'https://www.thinksrs.com/downloads/pdfs/applicationnotes/Residual%20Gas%20Analysis%20Basics.pdf', 'role': 'generic manufacturer measurement formalism, not identification of the actual Mykonos instrument'}],
     'conditional_density_example_cm3': 5e-5 * 133.32236842105263 / (1.380649e-23 * 300) / 1e6,
     'remedy': 'The approximate5e-5Torr context with ALL-H2O/300K is a named scenario, not a calibrated neutral mixture, exact boundary or surface coverage. Preserve uncertainty/identifiability; never initialize separate radicals merely from RGA fragment peaks.'},
    {'id': 'G01', 'severity': 'mandatory_integration_gate', 'class': 'host_owned_unimplemented_here',
     'finding': 'These pure kernels stage candidates only. They do not own collisions, accepted state, macro availability, exact source-pool debit, RNG, moving-target transforms, collision volume, current/energy deposition, rollback or restart. Water component checks do not audit another owner implementation.',
     'evidence': [ref('include/raikou/collision_data.h', 2, 3), ref('include/raikou/collision_data.h', 127, 148), {'path': str(R / 'data/collisions/h2o-collision-port-v1.json')}, {'path': str(R / 'data/collisions/h2o-collision-port-v7.json')}],
     'remedy': 'Use onePIC/MCC electron owner, one coherent source-rounded state/mass catalog, immutable dataset/site/state/mesh/time/event identity and actual3D volume. Commit matchedW source/electron/parent/product/sink/current/energy updates atomically. Count final reference/internal energies once; intermediate virtual-hole losses are not additional heat.'},
    {'id': 'G02', 'severity': 'mandatory_integration_gate', 'class': 'overlap_and_dispatch_policy_required',
     'finding': '39 stored cross-section observables are not39 disjoint event channels. Inclusive fragment/OH yields, effective momentum, total scattering, aggregate inelastic and stopping moments overlap selected events. The measured parent model and outer3 parent model both admit exactly1000eV; global model dispatch is host-owned.',
     'evidence': [ref('include/raikou/collision_data.h', 150, 178), ref('include/raikou/collision_data.h', 180, 199), {'path': str(A / 'observable-catalog.json')}, {'path': str(A / 'conditional-model-catalog.json')}],
     'remedy': 'Select one parent model at the shared endpoint and declare a single hazard partition. Treat OHv0..4 as oneDEA family branching, not five copies of its rate. Do not add inclusiveOH plusA1B1, IAM plus molecular elastic, original fullBEB plus orbital branches, or old2a1 plus the new partition. ReducedEMT/bend export is not a completePIC or waterF0 mechanism.'},
    {'id': 'L01', 'severity': 'explicit_limitation', 'class': 'honest_refusal_not_defect',
     'finding': 'Strict unsupported-energy and unselected-kernel exceptions, positive unmapped rates, missing-rate nulls and false full-mechanism flags are deliberate. O1S remains quarantined for a factor100 unit/text discrepancy. No rate is invented by subtraction or unsupported tail fill.',
     'evidence': [ref('include/raikou/collision_data.h', 113, 119), ref('include/raikou/collision_data.h', 145, 178), {'path': str(R / 'data/collisions/h2o-missing-rate-ledger-v7.json')}],
     'remedy': 'Preserve these refusals. A valid restricted component test or conditional comparison must name its scope; it cannot be presented as full breakdown.'},
]

rates = v7['at8084eV']
outer = load(W / 'event-attempt-008/validation.json')['tests']['outer_valence']['at8084eV']
four = sum(outer[k] for k in ['parent_m2', 'OHp_H_m2', 'Hp_OH_m2', 'unmapped_2a1_m2'])
missing = rates['positive_unmapped_prompt_Hp_m2'] + rates['positive_unmapped_prompt_Op_m2']
audit = {
    'schema': 'water-physics-applicability-audit-v1',
    'status': 'sealed_component_audit_full_integration_on_hold',
    'date': '2026-10-29',
    'published_water_commit': v7['remote_commit'],
    'published_water_dataset': v7['dataset'], 'published_water_dataset_sha256': v7['dataset_sha256'],
    'published_water_contract': v7['contract'], 'published_water_contract_sha256': v7['contract_sha256'],
    'published_water_source_archive': v7['source_archive'], 'published_water_source_archive_sha256': v7['source_archive_sha256'],
    'published_water_conditional_models': 10, 'stored_observables_including_diagnostics': 39,
    'native_assertions_in_existing_v7_receipt': 54190, 'atomic_source_binding_checks': 1810,
    'separate_H_working_candidate': {'dataset': h['dataset'], 'sha256': h['dataset_sha256'], 'checks': 515,
        'status': 'uncommitted_separate_snapshot_preserved_no_integration', 'snapshot': str(A / 'v8-working-snapshot.json')},
    'at_8084_470604111418_eV': {'rates_m2': rates, 'unmapped_fraction_of_2a1': 0.73,
        'unmapped_fraction_of_original_four_valence_BEB_only': missing / four,
        'mapped_fraction_of_original_four_valence_BEB_only': 1 - missing / four,
        'physical_all_shell_coverage_fraction': None,
        'warning': 'This fraction is a conditional partition of an8electron/fourvalence model, never all-shell/full-mechanism accuracy.'},
    'findings': findings,
    'numerical_defect_conclusion': 'No new numerical kernel defect confirmed by this static/hash audit; this is not an independent proof of every native kernel. The existing tests establish their declared computational invariants, not physical completeness or full-driver correctness.',
    'full_breakdown_ready': False,
    'new_native_tests_or_full_integration_runs': 0,
    'GPU_used': False, 'physical_time_advanced': 0, 'active_run_modified': False,
    'scope_boundary': 'Read-only source/provenance/receipt audit; writes confined to this separate audit package. No source/model substitution, commit, push, launch or production change.',
    'consolidation_entrypoint': str(A / 'handoff.json')
}
write('audit.json', audit)

# Bind exact external inputs as well as all immutable release contracts/datasets.
paths = set()
for pattern in ['h2o-gas-20261029*.json', 'h2o-collision-port-v*.json', 'h2o-missing-rate-ledger-v*.json', 'H-1s-ionization-20261029*.json']:
    paths.update((R / 'data/collisions').glob(pattern))
paths.update((R / 'include/raikou').glob('h2o*.h'))
paths.update([R / 'include/raikou/collision_data.h', R / 'include/raikou/hydrogen_fragment_ionization.h', W / 'api-handoff.json', W / 'coordination-handoff-v7.json', Path(v7['source_archive']), W / 'fragment-continuation-research.json'])
for attempt in ['event-attempt-008', 'fragment-attempt-001']:
    paths.update([W / attempt / 'validation.json', W / attempt / 'seal.json'])
for pattern in ['*.pdf', '*binding*.json']:
    paths.update((W / 'source-documents').glob(pattern))
paths.update((W / 'source-documents/nist-srd64-20261029').glob('*'))
rga = O / 'rga-channel-recovery'
paths.update([rga / 'figure_2_3_rga_recovery.json', rga / 'figure_2_3_rga_mass_channels.csv', rga / 'figure_2_3_rga_peak_channels.csv'])
thesis = (O.parent / 'work/smith-dissertation-2025.pdf')
assert sha(thesis) == '1d7b1f9c191a37e94f1cc3634587750f1af35b07eb13e0b94c60f3da47495ed7'
paths.add(thesis)
for name, expected in water['sources']['nist_srd64_v5']['csv_sha256'].items():
    assert sha(W / 'source-documents/nist-srd64-20261029' / name) == expected
for source, name in [('song2021', 'song-2021.pdf'), ('munoz2007', 'munoz-2007.pdf'), ('ram2010', 'ram-2010-water-dea.pdf'), ('murakami2012', 'murakami-2012-water.pdf'), ('montenegro2015', 'montenegro-2015-2a1.pdf'), ('melton1972', 'melton-1972.pdf')]:
    assert sha(W / 'source-documents' / name) == water['sources'][source]['pdf_sha256']
write('input-bindings.json', {str(p.resolve()): {'sha256': sha(p), 'bytes': p.stat().st_size} for p in sorted(paths) if p.is_file()})
shutil.copyfile(W / 'api-handoff.json', A / 'generic-handoff-as-audited.json')
write('web-verification.json', {
    'date': '2026-10-29',
    'primary_pages_read': [
        {'url': 'https://physics.nist.gov/PhysRefData/Ionization/intro.html', 'facts': 'GeneralQ-dependentBEB total formula; detailedBED SDCS separate; nonrelativistic10keV ceiling; generic5..20% statement toabout1keV.'},
        {'url': 'https://physics.nist.gov/cgi-bin/Ionization/ion_data.php?differential=Y&id=HI&initial=&ision=I', 'facts': 'H1s B=U=13.6057eV,N1,Q.5668,DoubleIonNo.'},
        {'url': 'https://physics.nist.gov/cgi-bin/Ionization/table.pl?ionization=H2O', 'facts': 'Original four water valence orbitals, no core orbital.'},
        {'url': 'https://www.thinksrs.com/downloads/pdfs/applicationnotes/Residual%20Gas%20Analysis%20Basics.pdf', 'facts': 'Quantitative mixture inference needs fragmentation factors and species sensitivity; generic manufacturer formalism only.'}
    ],
    'retrieval_failures': [
        {'url': 'https://iopscience.iop.org/article/10.1088/1742-6596/635/1/012019', 'outcome': 'web tool internal error; retained originalPDF and previously visually inspected pages4/5 used'},
        {'url': 'https://discovery.ucl.ac.uk/id/eprint/10128988/1/5.0035315.pdf', 'outcome': 'web tool timeout; retainedSHA-bound originalPDF used'}
    ],
    'new_visual_source_check': {'pdf': str(thesis), 'page_one_based': 44, 'printed_page': 19, 'render': str(A / 'smith-dissertation-page44.png'), 'full_page_visually_reviewed': True},
    'source_limitations': ['Song is a recommended compilation with primary lineage, not a claim all original experiments were independently digitized.', 'OriginalTan1978 was not retrieved; branching is traced through retainedMurakami2012 and its explicit citation.', 'Hargreaves2012 has retained web factual binding, not a newly archived localPDF.', 'Rosén/Jensen recombination sources are research leads, no admitted table/kernel.']
})

print(json.dumps({'audit': str(A / 'audit.json'), 'findings': len(findings), 'input_files_bound': len(paths), 'v8_snapshot_unchanged': True}))
