#!/usr/bin/env python3
"""Record immutable inputs for a contract review; performs no physics tests."""
from pathlib import Path
import difflib
import hashlib
import json
import resource
import shutil
import subprocess

resource.setrlimit(resource.RLIMIT_AS, (2 * 1024**3, 2 * 1024**3))
resource.setrlimit(resource.RLIMIT_CPU, (180, 180))
HERE = Path(__file__).resolve().parent
REPO = Path('/home/treefor/projects/runtimes/hall-coupled-transport-20261029/source')
MODULE = REPO / 'development/hall-canonical-driver-20261029'
OWNER = Path('/home/treefor/sim-workspace/powerflow/outputs/pic-electrode-heat-source-v3-20261029')
SNAPSHOT = OWNER / 'source-snapshot-v3'
NATIVE = Path('/home/treefor/projects/runtimes/pic-reaction-transactions-20261029/native-full-time')
EXPECTED_COMMIT = 'c2f81ab87ddc295e7bdad808d33135fa781b3942'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def copy_input(path, relative):
    target = HERE / 'inputs' / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(path, target)
    return {'original': str(path), 'copy': str(target.relative_to(HERE)),
            'bytes': target.stat().st_size, 'sha256': digest(target)}


def function_bytes(path, signature):
    text = path.read_text()
    start = text.index(signature)
    opening = text.index('{', start)
    depth = 1
    for i in range(opening + 1, len(text)):
        depth += (text[i] == '{') - (text[i] == '}')
        if depth == 0:
            return text[start:i + 1].encode()
    raise ValueError('Unterminated function')


assert not (HERE / 'inputs').exists(), 'Use a new review directory for new inputs'
head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip()
assert head == EXPECTED_COMMIT
assert not subprocess.check_output(['git', 'status', '--porcelain'], cwd=REPO, text=True)
records = []
seal = json.loads((OWNER / 'source-seal-v3.json').read_text())
for entry in seal['files']:
    path = SNAPSHOT / entry['path']
    assert digest(path) == entry['sha256']
    assert path.stat().st_size == entry['bytes']
    records.append(copy_input(path, Path('source-v3') / entry['path']))
for name in ['source-seal-v3.json', 'native-source-stage-port-v3.json', 'READY-v3.json']:
    records.append(copy_input(OWNER / name, Path('source-owner') / name))

comparisons = []
native_pairs = []
for entry in json.loads((SNAPSHOT / 'hall-inputs.json').read_text()):
    relative = Path(entry['path'].removeprefix('hall-owner/'))
    sealed = SNAPSHOT / entry['path']
    assert digest(sealed) == entry['sha256']
    if str(relative).startswith('include/raikou/json') or str(relative).startswith('development/'):
        current = REPO / relative
    else:
        current = MODULE / 'bindings/native-v6.2/source' / relative
    records.append(copy_input(current, Path('current-hall-binding') / relative))
    comparisons.append({'path': str(relative), 'sealed_sha256': entry['sha256'],
                        'current_sha256': digest(current),
                        'matches': digest(current) == entry['sha256']})
    if relative.name in ['mfem_gap_pic.h', 'mfem_gap_pic.cpp']:
        native_pairs.append((sealed, current))
        diff = ''.join(difflib.unified_diff(sealed.read_text().splitlines(True),
                                          current.read_text().splitlines(True),
                                          fromfile='source-v3/hall-owner/' + str(relative),
                                          tofile='canonical-native-v6.2/' + str(relative)))
        (HERE / (relative.name + '.diff')).write_text(diff)

module_files = ['axial_thermal_v2.hpp', 'axial_joint_v2.hpp', 'kinetic_boundary.hpp',
                'kinetic_boundary.cpp', 'motion_receipt.hpp', 'native_arena_motion.cpp',
                'native_joint_binding.hpp', 'MOTION-API.txt', 'KINETIC-BOUNDARY-CONTRACT.txt']
for name in module_files:
    records.append(copy_input(MODULE / name, Path('current-canonical-port') / name))
for name in ['source/include/raikou/mfem_pic_full_time.h', 'MOTION-RECEIPT-SCHEMA.txt']:
    records.append(copy_input(NATIVE / name, Path('native-owner-draft') / name))
for name in ['surface_source_binding.hpp', 'binding_seals.hpp', 'source-kernel/include/surface_emission.hpp']:
    records.append(copy_input(MODULE / 'bindings/heat-v2' / name, Path('heat-v2-reference') / name))

old_cpp, new_cpp = next(pair for pair in native_pairs if pair[0].suffix == '.cpp')
functions = []
for signature in ['PicVector pic_boris(', 'double pic_kinetic_energy_J(']:
    old, new = function_bytes(old_cpp, signature), function_bytes(new_cpp, signature)
    assert old == new
    functions.append({'signature': signature, 'byte_identical': True,
                      'sha256': hashlib.sha256(old).hexdigest()})
assert sum(x['matches'] for x in comparisons) == 13
assert len(comparisons) == 15
owner_results = json.loads((SNAPSHOT / 'evidence/results.json').read_text())
audit = {
    'schema': 'hall-source-v3-contract-review-inputs-v1',
    'kind': 'static_contract_review_not_runtime_qualification',
    'hall_commit': head,
    'git_clean_at_capture': True,
    'sealed_source_v3_files_verified': len(seal['files']),
    'source_v3_snapshot_unchanged': True,
    'hall_dependency_entries': len(comparisons),
    'dependencies_matching_current_hall_and_canonical_native': 13,
    'dependencies_requiring_owner_binding_update': 2,
    'dependency_comparison': comparisons,
    'unchanged_native_functions': functions,
    'source_owner_reported_qualification': owner_results,
    'source_owner_checks_rerun_here': False,
    'full_native_runtime_link_performed_here': False,
    'population_or_count_policy_selected': False,
    'native_ids_allocated': False,
    'accepted_state_changed': False,
    'resources': {'CPU': 6, 'nice': 15, 'threads': 1, 'GPU': False},
    'inputs': records,
}
(HERE / 'bindings-audit.json').write_text(json.dumps(audit, indent=2) + '\n')
print(json.dumps({k: audit[k] for k in [
    'kind', 'hall_commit', 'sealed_source_v3_files_verified',
    'dependencies_matching_current_hall_and_canonical_native',
    'dependencies_requiring_owner_binding_update']}))
