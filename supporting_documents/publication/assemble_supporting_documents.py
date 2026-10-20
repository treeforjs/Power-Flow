"""Publish reviewed historical handoffs and label preserved incomplete source work."""
import datetime
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import re
import shutil
import tarfile

REPO = Path('/workspace/Power-Flow')
WORK = Path(__file__).resolve().parent
OUT = REPO / 'supporting_documents'
ORIGINALS = Path('/workspace/handoffs')
PARTIAL = Path('/workspace/research_work/mitl-upstream-20261019T0953')
BASE = '5019b86aabf0472037b0fa3262922c3d0dcd2e95'
SPECIAL = 'mitl-source-window-resolution-20261019T0150.tar.gz'
ATTRIBUTION = re.compile(bytes.fromhex('5c62283f3a4f70656e41497c436861744750547c436f6465787c436c617564657c416e7468726f7069637c4750542d5b302d395d295c62'), re.I)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def relative(path):
    return str(path.relative_to(REPO))


def pin(path):
    return {'path': relative(path), 'bytes': path.stat().st_size, 'sha256': sha(path)}


def add_bytes(archive, name, data):
    item = tarfile.TarInfo(name)
    item.size = len(data)
    item.mode = 0o644
    item.mtime = 0
    archive.addfile(item, io.BytesIO(data))


assert not OUT.exists()
OUT.mkdir()
archives = sorted(ORIGINALS.glob('*.tar.gz'))
assert len(archives) == 25
entries = []
for original in archives:
    assert original.stat().st_size < 90 * 1024 * 1024
    parent_sha = sha(original)
    bundle = original.name.removesuffix('.tar.gz')
    published_name = (bundle + '-publication.tar.gz') if original.name == SPECIAL else original.name
    destination = OUT / 'archives' / published_name
    destination.parent.mkdir(exist_ok=True)
    documents, metadata, omissions, omitted_dirs = [], [], [], []
    with tarfile.open(original, 'r:gz') as archive:
        members = archive.getmembers()
        seen = set()
        for item in members:
            path = PurePosixPath(item.name)
            assert not path.is_absolute() and '..' not in path.parts
            assert item.name not in seen and (item.isfile() or item.isdir())
            seen.add(item.name)
        roots = {PurePosixPath(item.name).parts[0] for item in members}
        assert len(roots) == 1
        archive_root = next(iter(roots))
        top = OUT / 'handoffs' / bundle
        top.mkdir(parents=True)
        verification = None
        for item in members:
            parts = PurePosixPath(item.name).parts
            if not item.isfile() or len(parts) != 2:
                continue
            if parts[-1].endswith('.md') or parts[-1] in {'verification.json', 'MANIFEST.sha256', 'manifest.sha256'}:
                data = archive.extractfile(item).read()
                assert not ATTRIBUTION.search(data), item.name
                target = top / parts[-1]
                target.write_bytes(data)
                record = {**pin(target), 'archive_member': item.name}
                (documents if target.suffix == '.md' else metadata).append(record)
                if parts[-1] == 'verification.json':
                    verification = json.loads(data)
        assert documents and verification is not None
        if original.name != SPECIAL:
            shutil.copyfile(original, destination)
            assert sha(destination) == parent_sha
        else:
            retained = []
            for item in members:
                is_cache = '__pycache__' in PurePosixPath(item.name).parts or item.name.endswith('.pyc')
                if is_cache:
                    if item.isfile():
                        data = archive.extractfile(item).read()
                        omissions.append({'member': item.name, 'bytes': len(data),
                                          'sha256': hashlib.sha256(data).hexdigest(),
                                          'reason': 'Generated compiled Python cache; source remains retained.'})
                    else:
                        omitted_dirs.append(item.name)
                else:
                    retained.append(item)
            assert len(omissions) == 28
            note = {
                'schema': 'mitl-historical-handoff-publication-derivative-v1',
                'original_archive': original.name, 'original_archive_sha256': parent_sha,
                'original_archive_bytes': original.stat().st_size,
                'changes': 'Only generated compiled Python caches and their empty directories are omitted.',
                'omitted_files': omissions, 'omitted_directories': omitted_dirs,
                'scientific_source_data_or_results_changed': False,
                'historical_verification': 'Original verification/manifest bytes are retained as historical records of the original archive.',
                'current_payload_manifest': 'PUBLICATION_MANIFEST.sha256',
                'scientific_tests_or_models_rerun': False,
            }
            rows = []
            with tarfile.open(destination, 'w:gz', format=tarfile.PAX_FORMAT) as published:
                for item in retained:
                    if item.isfile():
                        data = archive.extractfile(item).read()
                        published.addfile(item, io.BytesIO(data))
                        rows.append((hashlib.sha256(data).hexdigest(), item.name))
                    else:
                        published.addfile(item)
                data = (json.dumps(note, indent=2) + '\n').encode()
                note_name = archive_root + '/PUBLICATION.json'
                add_bytes(published, note_name, data)
                rows.append((hashlib.sha256(data).hexdigest(), note_name))
                manifest = ''.join(digest + '  ' + name + '\n' for digest, name in sorted(rows)).encode()
                add_bytes(published, archive_root + '/PUBLICATION_MANIFEST.sha256', manifest)
            save(top / 'PUBLICATION.json', note)
            metadata.append(pin(top / 'PUBLICATION.json'))
        entries.append({
            'bundle': bundle, 'original_archive': original.name,
            'original_sha256': parent_sha, 'original_bytes': original.stat().st_size,
            'publication_archive': pin(destination), 'publication_derivative': bool(omissions),
            'omitted_compiled_files': omissions, 'omitted_cache_directories': omitted_dirs,
            'readable_handoff_documents': documents, 'readable_verification_records': metadata,
            'source_commit': verification.get('source_commit'),
            'base_commit': verification.get('base_commit'),
            'verification_schema': verification.get('schema'),
            'historical_baseline_required': True,
        })

# Preserve project-generated, unfinished work separately from completed packages.
partial_out = OUT / 'in_progress' / PARTIAL.name
partial_out.mkdir(parents=True)
partial_files = []
for source in sorted(PARTIAL.rglob('*')):
    if not source.is_file() or '__pycache__' in source.parts:
        continue
    rel = source.relative_to(PARTIAL)
    copy = source.suffix == '.py' or str(rel) in {
        'commands.jsonl', 'scope.json', 'lease.json', 'citation-review/completion.json',
        'source-acquisition/catalog-links.json', 'delivery-preparation/contract.json',
        'delivery-preparation/draft-notes.json',
    }
    if copy:
        data = source.read_bytes()
        assert not ATTRIBUTION.search(data), str(source)
        target = partial_out / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
        partial_files.append({**pin(target), 'original_relative_path': str(rel),
                              'original_sha256': sha(source), 'status': 'incomplete_unexecuted_or_partially_executed'})

responses = []
for source in sorted((PARTIAL / 'source-acquisition').glob('*.json')):
    document = json.loads(source.read_text())
    if document.get('schema') != 'mitl-bounded-source-response-v1':
        continue
    response = document['response']
    responses.append({'input_record': source.name, 'input_record_sha256': sha(source),
                      'started_utc': document['started_utc'], 'finished_utc': document['finished_utc'],
                      'command': document['command'], 'exit_code': document['exit_code'],
                      'http_code': response.get('http_code'), 'effective_url': response.get('url_effective'),
                      'content_type': response.get('content_type'), 'files': document['files'],
                      'body_is_PDF': document['body_is_PDF']})
assert len(responses) == 4 and all(row['http_code'] == 200 for row in responses)
save(partial_out / 'acquisition_records.json', {
    'schema': 'mitl-incomplete-source-acquisition-records-v1', 'status': 'incomplete',
    'responses': responses, 'completed_primary_report_assessment': False,
    'projection': 'Retains request/response identity and hashes; connection diagnostics and third-party response bodies remain in the local source cache.',
})

discovery = []
for name in ('initial-discovery-response.json', 'official-discovery-response.json'):
    source = PARTIAL / 'source-acquisition' / name
    record = json.loads(source.read_text())
    discovery.append({'original_file': name, 'original_sha256': sha(source),
                      'utc': record['utc'], 'queries': record['args']['system2_search_query']})
save(partial_out / 'discovery_records.json', {'status': 'incomplete', 'discovery_requests': discovery})

raw_cache = Path('/workspace/research_sources/npl_processability_1996')
raw_omissions = [{'local_path': str(path), 'bytes': path.stat().st_size, 'sha256': sha(path),
                  'reason': 'Third-party source, OCR or source-page image retained locally; use the cited official source.'}
                 for path in sorted(raw_cache.rglob('*')) if path.is_file()]
save(partial_out / 'source_cache_index.json', {'raw_source_files_republished': False, 'records': raw_omissions})

(partial_out / 'README.md').write_text('''# Preserved incomplete source acquisition

Status: **incomplete**. This is the interrupted 19 October 2026 NPL-source pass,
not a completed or validated handoff. No property input was adopted, and no
scientific model or source integration was repeated.

The sought source was the 1996 McCormick/Brooks NPL MTS processability report
cited by Mills 2004, reference 8, specifically for liquid 316 density. Four
official HTTP requests completed. A related 1998 CMMT(A)140 methods PDF was
acquired; it is a different report. A catalog/cover identity discrepancy and
alloy coverage still require completion of the page review. The original 1996
report remains unlocated. The saved citation review is separate from that
unfinished new-source review.

The acquisition and review scripts, exact command records, request identities,
source hashes and unexecuted packet-helper drafts are preserved here. Absolute
paths describe the original workspace. These drafts are not a runnable release.
The packet-helper review identified two unresolved checks: validation of every
live command-record suffix, and independent required archive/evidence coverage
bound to the contract. Correct them before any packet helper is executed.

Third-party PDFs, source HTML, source-page images, OCR and search-result bodies
remain in the local source cache. Official URLs and their recorded sizes/hashes
are retained in acquisition_records.json and source_cache_index.json.

Next useful experimental source lead: Lamppa et al., *Assessment of Electrode
Contamination Mitigation at 0.5 MA Scale*, SAND-2021-12691, DOI
[10.2172/1825219](https://doi.org/10.2172/1825219),
[official full text](https://www.osti.gov/servlets/purl/1825219).
It is dissertation reference 58 and is linked there to the Mykonos MP3 platform.
Its full text has not yet been acquired or assessed in this workstream; no
same-sample geometry, calibration or surface parameter is inferred.
''')

snapshot = json.loads((WORK / 'status-before-lease.json').read_text())
for key in ('active_lease_owner', 'active_editing_lease_until_utc', 'reserved_paths', 'active_phase'):
    snapshot.pop(key, None)
snapshot['snapshot_qualification'] = 'Historical scientific checkpoint before this publication pass; operational leases omitted. The 0953 source acquisition remains incomplete.'
save(OUT / 'research_checkpoint.json', snapshot)

catalog = {
    'schema': 'mitl-supporting-document-publication-catalog-v1',
    'assembled_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'completed_code_baseline_commit': BASE, 'completed_archive_count': len(entries),
    'original_archive_bytes': sum(entry['original_bytes'] for entry in entries),
    'artifacts': entries,
    'preserved_incomplete_work': {'status': 'incomplete', 'pass_id': PARTIAL.name,
        'readme': relative(partial_out / 'README.md'), 'completed_handoff': False,
        'preserved_files': partial_files, 'acquisition_records': pin(partial_out / 'acquisition_records.json'),
        'source_cache_index': pin(partial_out / 'source_cache_index.json')},
    'source_document_policy': 'Project reports and generated verification figures are retained in the historical packages. Third-party full papers, licensed source, credentials, private implementations and build directories are not publication inputs; citations and source hashes remain available.',
    'verification_record': 'supporting_documents/publication/verification.json',
    'scientific_models_tests_or_recovery_reexecuted_by_assembly': False,
}
save(OUT / 'catalog.json', catalog)

lines = ['# Supporting documents and historical handoffs', '',
         'This catalog publishes all 25 completed handoff packages and their readable',
         'reports. The current source and research documents are in the repository.',
         'Each package preserves its own baseline, source versions, evidence and',
         'limitations; patch subsets are not standalone simulation packages.', '',
         '[Project handoff](../HANDOFF.md) · [Machine-readable catalog](catalog.json) ·',
         '[Publication verification](publication/verification.json)', '',
         '| Package | Readable handoff | Evidence archive |',
         '| --- | --- | --- |']
for entry in entries:
    document = next((item for item in entry['readable_handoff_documents'] if '/systems-' in item['path']), entry['readable_handoff_documents'][0])
    doc_link = str(Path(document['path']).relative_to('supporting_documents'))
    arc_link = str(Path(entry['publication_archive']['path']).relative_to('supporting_documents'))
    label = entry['bundle'].removeprefix('mitl-')
    lines.append(f'| {label} | [Report]({doc_link}) | [Download]({arc_link}) |')
lines += ['', 'Twenty-four archives are exact byte copies. The 0150 source-window package has',
          'a labelled publication derivative omitting 28 generated Python bytecode files.',
          'Its original archive hash and omission map are retained, along with the',
          'unchanged historical verification records and a new publication manifest.', '',
          'For fresh 0737 independent-reference replay, retain the separately included',
          '0710 archive containing its original seed checkpoints. All recorded baseline',
          'commits are available in this branch history.', '',
          f'[Preserved incomplete source work](in_progress/{PARTIAL.name}/README.md) is',
          'labelled separately. [The scientific checkpoint](research_checkpoint.json)',
          'retains the completed queue and evidence context. No unfinished pass is',
          'promoted to a completed calculation.', '']
(OUT / 'README.md').write_text('\n'.join(lines))

checksums = []
for entry in entries:
    record = entry['publication_archive']
    checksums.append(record['sha256'] + '  ' + str(Path(record['path']).relative_to('supporting_documents')))
(OUT / 'ARCHIVE_SHA256SUMS').write_text('\n'.join(checksums) + '\n')
publication = OUT / 'publication'
publication.mkdir()
shutil.copyfile(__file__, publication / 'assemble_supporting_documents.py')
report = {'schema': 'mitl-supporting-document-assembly-v1',
          'catalog': pin(OUT / 'catalog.json'), 'copied_archives': 24, 'publication_derivatives': 1,
          'total_archives': len(entries), 'publication_archive_bytes': sum(entry['publication_archive']['bytes'] for entry in entries),
          'compiled_files_omitted': 28, 'incomplete_source_pass_completed': False,
          'scientific_models_or_tests_executed': False}
save(WORK / 'supporting-document-assembly.json', report)
print(json.dumps(report, indent=2))
