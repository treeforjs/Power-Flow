"""Immutable stage provenance and fail-closed technical sealing.

These receipts establish which code/inputs/binaries were exercised. They do
not supply scientific approval for a changed criterion or model assumptions.
"""
import hashlib
import importlib.util
import json
from pathlib import Path

REQUIRED = ('inputs', 'tables', 'policy', 'build', 'core', 'sanitized', 'gpu', 'cost', 'bench', 'swarm')
BINARIES = ('test_core', 'test_core_sanitized', 'gpu_parity', 'swarm_mc')
SCHEMA = 'pfc-qualification-bound-stages-v1'


def need(condition, message):
    if not condition:
        raise ValueError('qualification integrity: ' + message)


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def capture(root, dependencies):
    root = Path(root)
    files = {}
    for directory in ('include', 'device', 'tests', 'tools', 'data/extracted', 'data/tables'):
        for p in sorted((root / directory).rglob('*')):
            if p.is_file() and '__pycache__' not in p.parts:
                files[str(p.relative_to(root))] = sha(p)
    need(files, 'empty source tree')
    files['run_checks.py'] = sha(root / 'run_checks.py')
    manifest_path = root / 'data/tables/chemistry-v3.2-manifest.json'
    manifest = json.loads(manifest_path.read_text())
    inputs = {}
    for name, expected in manifest['input_files_sha256'].items():
        actual = sha(root / name)
        need(actual == expected, 'manifest input changed: ' + name)
        inputs[name] = actual
    # The swarm stage also reads these source sets directly, outside its production manifest.
    for p in sorted((root / 'data/sources').glob('lxcat-*/Cross section.txt')):
        inputs[str(p.relative_to(root))] = sha(p)
    dep = {str(p): sha(p) for p in sorted(dependencies)}
    result = {'source_files': files, 'table_inputs': inputs, 'dependencies': dep}
    result['sha256'] = digest(result)
    return result


def binaries(root):
    return {name: sha(Path(root) / 'build/bin' / name) for name in BINARIES
            if (Path(root) / 'build/bin' / name).is_file()}


def begin_stage(root, dependencies, receipt, name):
    need(receipt.get('schema') == SCHEMA, 'legacy or unbound receipt requires new qualification')
    need(name in REQUIRED, 'unexpected stage ' + name)
    old = receipt.get('stages', [])
    need(name == 'inputs' or not any(s['stage'] == name for s in old), 'duplicate stage ' + name)
    before = capture(root, dependencies)
    if old:
        need(old[0]['binding']['before'] == before, 'code/inputs/dependencies differ from prior stage')
    return {'schema': 'pfc-stage-binding-v1', 'before': before, 'binaries_before': binaries(root)}


def finish_stage(root, dependencies, before):
    result = dict(before)
    result['after'] = capture(root, dependencies)
    need(result['before'] == result['after'], 'code/inputs changed during stage')
    result['binaries_after'] = binaries(root)
    return result


def verify_receipt(root, receipt, output_dir):
    """No writes. Legacy, incomplete, stale, or mixed qualification refuses."""
    root, output_dir = Path(root), Path(output_dir)
    need(receipt.get('schema') == SCHEMA, 'receipt lacks bound stage schema')
    stages = receipt.get('stages', [])
    names = [s.get('stage') for s in stages]
    need(set(names) == set(REQUIRED), 'missing/extra stages: ' + repr(sorted(set(REQUIRED) ^ set(names))))
    need(names and names[0] == 'inputs', 'inputs must precede qualification')
    for name in REQUIRED[1:]:
        need(names.count(name) == 1, 'duplicate stage ' + name)
    order = [name for name in names if name != 'inputs']
    need(order == list(REQUIRED[1:]), 'stage order differs from required qualification')
    for stage in stages:
        need(stage.get('passed') is True, 'failed stage: ' + str(stage.get('stage')))
        entries = stage.get('entries', [])
        need(entries and all(type(e.get('exit')) is int and e['exit'] == 0 for e in entries), 'missing/failed entry')
    initial = stages[0]['binding']['before']
    need(initial == capture(root, initial['dependencies']), 'current code/inputs differ from qualification')
    current_binaries = binaries(root)
    need(set(current_binaries) == set(BINARIES), 'missing required binaries')
    build = next(s for s in stages if s['stage'] == 'build')
    built = build['binding']['binaries_after']
    need(built == current_binaries, 'built binaries differ from current binaries')
    after_build = False
    for stage in stages:
        binding = stage.get('binding', {})
        need(binding.get('schema') == 'pfc-stage-binding-v1', 'unbound stage ' + stage['stage'])
        need(binding.get('before') == initial and binding.get('after') == initial, 'stale/mixed stage source identity')
        if stage['stage'] == 'build':
            after_build = True
            need(binding.get('binaries_after') == built, 'build outputs mismatch')
        elif after_build:
            need(binding.get('binaries_before') == built and binding.get('binaries_after') == built,
                 'tested binaries differ from build: ' + stage['stage'])
    table = next(s for s in stages if s['stage'] == 'tables')
    policy = next(s for s in stages if s['stage'] == 'policy')
    need(table['blob_sha256'] == sha(root/'data/tables/chemistry-v3.2.blob'), 'table blob mismatch')
    need(policy['policy_sha256'] == sha(root/'data/tables/chemistry-v3.2-coverage-policy.json'), 'policy mismatch')
    swarm = next(s for s in stages if s['stage'] == 'swarm')
    spec = importlib.util.spec_from_file_location('qualification_job_proof', root/'tools/swarm_job_integrity.py')
    proof_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(proof_module)
    contract = proof_module.V.Contract(root/'tools/swarm-contract.json', root)
    for directory in ('swarm', 'swarm-multiseed'):
        need(swarm.get('job_plan_sha256', {}).get(directory) == sha(output_dir/directory/'jobs.json'), 'job plan changed')
        jobs = json.loads((output_dir/directory/'jobs.json').read_text())
        completed = json.loads((output_dir/directory/'jobs-receipt.json').read_text())['jobs']
        need(jobs and len(jobs) == len(completed), 'incomplete job receipt: '+directory)
        names = [j['name'] for j in jobs]
        need(len(names) == len(set(names)) and set(names) == {j['name'] for j in completed}, 'duplicate/missing job identity')
        need(all(type(j.get('exit')) is int and j['exit'] == 0 and j.get('completion_verified') is True
                 for j in completed), 'job completion not verified')
        if directory == 'swarm-multiseed':
            expected = {f'{label}-{seed}' for label in ('v3','v32') for seed in range(20261101,20261107)}
            need(set(names) == expected, 'multiseed job plan incomplete')
        else:
            expected = {f'{label}-{k}' for label in ('v3','v3-iso','lxcat-muroranit','lxcat-triniti',
                                                   'lxcat-hayashi','lxcat-phelps','lxcat-morgan') for k in range(3)}
            need(set(names) == expected, 'legacy comparison job plan incomplete')
        for job in jobs:
            need(Path(job['argv'][0]).resolve() == (root/'build/bin/swarm_mc').resolve(), 'swarm binary path differs')
            need(Path(job['output']).resolve() == (output_dir/directory/(job['name']+'.json')).resolve(), 'job output ownership differs')
            label, suffix = job['name'].rsplit('-',1)
            if directory == 'swarm-multiseed':
                expected_blob = contract.data['sets'][label]['blob_sha256']
                expected_energies = contract.energies
                expected_seed = int(suffix)
            else:
                expected_blob = swarm['sets'][label]
                ens = [30,60,100,150,200,300,500,700,1000,1500,2000,3000]
                expected_energies = ens[4*int(suffix):4*(int(suffix)+1)]
                expected_seed = 20261031
            args = job['argv']
            need(sha(args[1]) == expected_blob and [float(e) for e in args[3].split(',')] == expected_energies and
                 float(args[4]) == contract.data['N_m3'] and float(args[5]) == contract.data['T_gas_K'] and
                 int(args[6]) == contract.data['electrons'] and int(args[7]) == expected_seed, 'job plan/configuration differs')
            need(proof_module.completed(job, initial['sha256']), 'missing/stale per-job proof: '+job['name'])
    artifact = output_dir / 'SWARM-REVALIDATION.json'
    need(swarm.get('revalidation_sha256') == sha(artifact), 'swarm result missing/changed')
    result, _ = proof_module.V.stable_json(artifact)
    need(result.get('complete') is True and result.get('passed') is True and result.get('judged') == 43
         and len(result.get('rows', [])) == 50, 'swarm result incomplete/failed')
    contract_hash = sha(root/'tools/swarm-contract.json')
    need(result.get('contract_sha256') == contract_hash, 'swarm contract identity differs')
    current_ref = proof_module.V.summarize(output_dir/'swarm-multiseed','v3',contract)
    current_new = proof_module.V.summarize(output_dir/'swarm-multiseed','v32',contract)
    proof_module.V.validate_comparison(result, current_new, current_ref, contract)
    need(swarm.get('scientific_criterion_status') == result.get('scientific_criterion_status') and
         bool(result.get('scientific_criterion_status')), 'scientific criterion status missing')
    return {'schema': 'pfc-technical-qualification-v1', 'binding_sha256': initial['sha256'],
            'binaries': built, 'required_stages': list(REQUIRED), 'swarm_result_sha256': sha(artifact),
            'production_qualified': False, 'physics_model_validated': False,
            'scientific_criterion_status': result['scientific_criterion_status'],
            'scope': 'complete execution and source integrity only; scientific criterion/model review remains separate'}
