"""Claim a bounded primary-source acquisition under the saved lease lock."""
import datetime
import fcntl
import hashlib
import json
from pathlib import Path
import subprocess

root = Path('/workspace/research_work')
repo = Path('/workspace/Power-Flow')
status_path = root / 'mitl-upstream-status.json'
pass_id = 'mitl-upstream-20261019T0953'
work = root / pass_id
base = '5019b86aabf0472037b0fa3262922c3d0dcd2e95'
source = Path('/workspace/research_sources/ss304_properties/new_mills2004/mills_2004_pages_text.txt')
scope = [
    'data/materials/npl_processability_1996_source_manifest.json',
    'docs/npl_processability_1996_source_assessment.md',
    str(work),
    '/workspace/research_sources/npl_processability_1996',
    '/workspace/handoffs/mitl-npl-processability-20261019T0953*',
]
task = ('Discover and acquire the previously unassessed McCormick/Brooks 1996 NPL MTS '
        'processability primary report through legitimate official routes; establish '
        'identity, access and alloy coverage before any grade-specific use.')
with (root / 'mitl-upstream-status.lock').open('a+') as lock:
    fcntl.flock(lock, fcntl.LOCK_EX)
    before = status_path.read_bytes()
    status = json.loads(before)
    now = datetime.datetime.now(datetime.timezone.utc)
    until = status.get('active_editing_lease_until_utc')
    if until and datetime.datetime.fromisoformat(until) > now:
        raise RuntimeError('Another editing lease is unexpired: ' + str(status.get('active_lease_owner')))
    if status.get('reserved_paths'):
        raise RuntimeError('Unreleased reserved paths require resolution before a new pass')
    assert subprocess.check_output(['git', 'branch', '--show-current'], cwd=repo, text=True).strip() == status['branch'] == 'research/flash-publication-foundation'
    assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip() == status['last_complete_commit'] == base
    assert not subprocess.check_output(['git', 'status', '--porcelain'], cwd=repo, text=True).strip()
    handoff = Path(status['handoff_document'])
    assert handoff.read_text().startswith('# Manufactured positive joint boundary quadrature')
    assert hashlib.sha256(Path(status['handoff']).read_bytes()).hexdigest() == status['handoff_sha256']
    assert 'McCormick and R. F. Brooks' in source.read_text()
    assert not work.exists()
    expires = now + datetime.timedelta(hours=2)
    work.mkdir()
    (work / 'status-before-lease.json').write_bytes(before)
    record = {
        'utc': now.isoformat(), 'event': 'lease_claimed', 'pass_id': pass_id,
        'baseline_commit': base, 'status_before_sha256': hashlib.sha256(before).hexdigest(),
        'prior_handoff': status['handoff'], 'prior_handoff_sha256': status['handoff_sha256'],
        'prior_handoff_document': str(handoff),
        'prior_handoff_document_sha256': hashlib.sha256(handoff.read_bytes()).hexdigest(),
        'task': task, 'saved_queue_priority': 5, 'reserved_paths': scope,
        'lease_until_utc': expires.isoformat(),
        'new_lead': {'citing_source_DOI': '10.2355/isijinternational.44.1661',
                     'reference': 8, 'printed_page': 1668, 'PDF_page': 8,
                     'cached_text': str(source), 'cached_text_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
                     'prior_acquisition_or_report_assessment_found': False,
                     'citing_use': 'Liquid 316 density; 304 coverage unverified'},
        'limits': 'Source identity/access/alloy assessment only. No property adapter, physical parameter adoption, numerical replay, solver or Quark changes; no blocked source routes retried. Up to two new discovery search calls and six official catalog/landing/document requests unless a concrete new location justifies a recorded extension.',
    }
    (work / 'lease.json').write_text(json.dumps(record, indent=2) + '\n')
    status.update(active_lease_owner=pass_id, active_editing_lease_until_utc=expires.isoformat(),
                  reserved_paths=scope, active_phase=task, last_update_utc=now.isoformat())
    temporary = status_path.with_suffix('.json.claim-new')
    assert not temporary.exists()
    temporary.write_text(json.dumps(status, indent=2) + '\n')
    temporary.replace(status_path)
    with (root / 'mitl-upstream-journal.jsonl').open('a') as stream:
        stream.write(json.dumps(record) + '\n')
    print(json.dumps(record, indent=2))
