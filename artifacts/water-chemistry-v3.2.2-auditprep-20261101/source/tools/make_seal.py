#!/usr/bin/env python3
"""Write SEAL.json of chemistry v3.2 from a passed qualification receipt.

usage: make_seal.py RECEIPT.json [COST.json]
Code identity: sha256 of the compact sorted JSON map {path: sha256} over include/, device/, tests/ (as v3.1). Files:
every file of the component except build/, __pycache__/ and SEAL.json. The receipt must have passed every stage, and
its blob and policy identities must equal the shipped files.
"""
import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path

R = Path(__file__).resolve().parent.parent
_qspec = importlib.util.spec_from_file_location('qualification_integrity', R / 'tools/qualification_integrity.py')
QI = importlib.util.module_from_spec(_qspec)
_qspec.loader.exec_module(QI)


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def main(receipt_path, cost_path=None):
    rc = json.load(open(receipt_path))
    qualified = QI.verify_receipt(R, rc, Path(receipt_path).resolve().parent)
    blob = sha(R / 'data/tables/chemistry-v3.2.blob')
    pol = sha(R / 'data/tables/chemistry-v3.2-coverage-policy.json')
    st = {s['stage']: s for s in rc['stages']}
    if st['tables']['blob_sha256'] != blob or st['policy']['policy_sha256'] != pol:
        raise SystemExit('receipt identities differ from the shipped tables')
    files = {}
    for p in sorted(R.rglob('*')):
        if not p.is_file() or 'build' in p.relative_to(R).parts or '__pycache__' in p.parts or p.name == 'SEAL.json':
            continue
        files[str(p.relative_to(R))] = sha(p)
    code = {k: v for k, v in files.items() if k.split('/')[0] in ('include', 'device', 'tests')}
    code_id = hashlib.sha256(json.dumps(code, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    seal = {
        'schema': 'pfc-chem-v3-seal-1', 'version': '3.2', 'created_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        'code_identity': code_id,
        'technical_qualification': qualified,
        'production_qualified': False,
        'physics_model_validated': False,
        'code_identity_rule': 'sha256 of the compact sorted JSON map {path: sha256} over include/, device/, tests/',
        'table_blob_sha256': blob, 'coverage_policy_sha256': pol, 'qualification_receipt': str(Path(receipt_path).resolve()),
        'predecessor': '/home/treefor/projects/runtimes/pic-water-chemistry-v3.1-20261031 (SEAL.json there)',
        'pinned_dependencies': st['inputs']['pinned_dependencies'],
        'build_isa_flags': st['build'].get('PF_BUILD_ISA_FLAGS', ''),
        'cost_per_collision': json.load(open(cost_path)) if cost_path else st.get('cost', {}).get('cost_per_collision'),
        'files': files,
    }
    QI.need(qualified == QI.verify_receipt(R, rc, Path(receipt_path).resolve().parent), 'inputs changed during sealing')
    temporary = R / 'SEAL.json.partial'
    with open(temporary, 'x') as stream:
        json.dump(seal, stream, indent=1, allow_nan=False)
        stream.write('\n')
        stream.flush()
        __import__('os').fsync(stream.fileno())
    temporary.replace(R / 'SEAL.json')
    print('SEAL.json: code', code_id, 'blob', blob, 'policy', pol, 'files', len(files))


if __name__ == '__main__':
    main(*sys.argv[1:3])
