"""Independently check publication copies and local navigation; execute no models."""
from __future__ import annotations

import argparse
import ast
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import sys
import tarfile
import traceback
import xml.etree.ElementTree as ET

from markdown_links import verify_documents

REPO = Path('/workspace/Power-Flow')
WORK = Path('/workspace/research_work/mitl-github-push-20261020T0601')
OUT = WORK / 'publication-verification'
SUPPORT = REPO / 'supporting_documents'
INCOMPLETE = Path('/workspace/research_work/mitl-upstream-20261019T0953')
SELF_REPORT = SUPPORT / 'publication/verification.json'


def require(condition, detail):
    if not condition:
        raise ValueError(detail)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def pairs(items):
    result = {}
    for key, value in items:
        require(key not in result, f'Duplicate JSON key: {key}')
        result[key] = value
    return result


def parsed(data):
    def invalid(value):
        raise ValueError('Nonfinite JSON token: '+value)
    return json.loads(data, object_pairs_hook=pairs, parse_constant=invalid)


def relative(root, name):
    part = PurePosixPath(name)
    require(not part.is_absolute() and '..' not in part.parts and '\\' not in name,
            'Unsafe relative path: '+name)
    result = root / name
    require(result.resolve().is_relative_to(root.resolve()), 'Escaping path: '+name)
    require(not result.is_symlink(), 'Unexpected symbolic link: '+name)
    return result


PINS = {}


def pin(path):
    path = Path(path)
    if path not in PINS:
        require(path.is_file() and not path.is_symlink(), 'Missing regular file: '+str(path))
        h = hashlib.sha256()
        length = 0
        with path.open('rb') as stream:
            for block in iter(lambda: stream.read(1048576), b''):
                h.update(block)
                length += len(block)
        PINS[path] = {'bytes': length, 'sha256': h.hexdigest()}
    return PINS[path]


def check_pin(path, item, *, prefix=''):
    actual = pin(path)
    require(actual['sha256'] == item[prefix+'sha256'], 'SHA mismatch: '+str(path))
    if prefix+'bytes' in item:
        require(actual['bytes'] == item[prefix+'bytes'], 'Size mismatch: '+str(path))
    return actual


def checksum_rows(data):
    result = {}
    require(data.endswith(b'\n'), 'Checksum list lacks final newline')
    for line in data.decode('utf-8').splitlines():
        match = re.fullmatch(r'([0-9a-f]{64})  (.+)', line)
        require(match is not None, 'Invalid checksum row')
        value, name = match.groups()
        require(name not in result, 'Duplicate checksum member: '+name)
        result[name] = value
    return result


def tar_members(archive):
    members = {}
    for member in archive.getmembers():
        name = member.name.rstrip('/')
        require(name not in members, 'Duplicate archive member: '+name)
        part = PurePosixPath(name)
        require(not part.is_absolute() and '..' not in part.parts and '\\' not in name,
                'Unsafe archive member: '+name)
        require(member.isfile() or member.isdir(), 'Unexpected archive member type: '+name)
        members[name] = member
    return members


def member_bytes(archive, member):
    require(member.isfile(), 'Expected regular archive member: '+member.name)
    stream = archive.extractfile(member)
    require(stream is not None, 'Unreadable archive member: '+member.name)
    return stream.read()


def verify_derivative(entry, original, published, original_members, published_members):
    omissions = {item['member']: item for item in entry['omitted_compiled_files']}
    require(len(omissions) == len(entry['omitted_compiled_files']) == 28, 'Wrong compiled omission count')
    omitted_dirs = set(entry['omitted_cache_directories'])
    require(len(omitted_dirs) == 3, 'Wrong cache-directory omission count')
    require(set(original_members)-set(published_members) == set(omissions)|omitted_dirs,
            'Derivative removed undeclared members')
    root = entry['bundle']
    note_name = root+'/PUBLICATION.json'
    manifest_name = root+'/PUBLICATION_MANIFEST.sha256'
    require(set(published_members)-set(original_members) == {note_name, manifest_name},
            'Unexpected derivative additions')
    for name, expected in omissions.items():
        require(name.endswith('.pyc') and '__pycache__' in PurePosixPath(name).parts,
                'Omission is not compiled cache: '+name)
        data = member_bytes(original, original_members[name])
        require(len(data) == expected['bytes'] and digest(data) == expected['sha256'],
                'Omitted member differs from recorded original: '+name)
    for name in omitted_dirs:
        require(original_members[name].isdir() and '__pycache__' in PurePosixPath(name).parts,
                'Unexpected omitted directory: '+name)
    retained = {}
    for name, member in published_members.items():
        if member.isfile():
            data = member_bytes(published, member)
            retained[name] = digest(data)
            if name in original_members:
                require(data == member_bytes(original, original_members[name]),
                        'Changed retained historical bytes: '+name)
    note = parsed(member_bytes(published, published_members[note_name]))
    require(note['schema'] == 'mitl-historical-handoff-publication-derivative-v1', 'Derivative schema')
    require(note['original_archive'] == entry['original_archive'] and
            note['original_archive_sha256'] == entry['original_sha256'] and
            note['original_archive_bytes'] == entry['original_bytes'], 'Derivative original association')
    require(note['omitted_files'] == entry['omitted_compiled_files'] and
            note['omitted_directories'] == entry['omitted_cache_directories'], 'Derivative omission association')
    require(note['scientific_source_data_or_results_changed'] is False and
            note['scientific_tests_or_models_rerun'] is False, 'Derivative qualification')
    manifest = checksum_rows(member_bytes(published, published_members[manifest_name]))
    require(manifest == {name: value for name, value in retained.items() if name != manifest_name},
            'Publication manifest does not cover exact retained regular payload')
    return {'retained_regular_members': len(retained)-2, 'omitted_compiled_files': len(omissions),
            'omitted_cache_directories': len(omitted_dirs), 'publication_manifest_entries': len(manifest)}


def verify_archives(catalog):
    rows, expected_sums, document_count, record_count = [], {}, 0, 0
    require(len(catalog['artifacts']) == catalog['completed_archive_count'] == 25, 'Expected 25 completed packets')
    bundles = [entry['bundle'] for entry in catalog['artifacts']]
    require(len(set(bundles)) == 25 and not any('0953' in name for name in bundles), 'Completed bundle identities')
    for entry in catalog['artifacts']:
        original_path = relative(Path('/workspace/handoffs'), entry['original_archive'])
        check_pin(original_path, {'sha256':entry['original_sha256'], 'bytes':entry['original_bytes']})
        public = entry['publication_archive']
        published_path = relative(REPO, public['path'])
        check_pin(published_path, public)
        expected_sums[str(published_path.relative_to(SUPPORT))] = public['sha256']
        derivative = entry['publication_derivative']
        require(type(derivative) is bool, 'Invalid derivative flag')
        if not derivative:
            require(pin(original_path) == pin(published_path), 'Archive is not exact copied bytes')
            require(not entry['omitted_compiled_files'] and not entry['omitted_cache_directories'],
                    'Exact copy declares omissions')
        with tarfile.open(original_path, 'r:gz') as original, tarfile.open(published_path, 'r:gz') as published:
            old_members, new_members = tar_members(original), tar_members(published)
            require(not any(name.endswith('.pyc') or '__pycache__' in PurePosixPath(name).parts
                            for name in new_members), 'Published compiled cache')
            derivative_record = verify_derivative(entry, original, published, old_members, new_members) if derivative else None
            for kind in ('readable_handoff_documents', 'readable_verification_records'):
                require(entry[kind], 'Missing readable top-level documents: '+entry['bundle'])
                for item in entry[kind]:
                    path = relative(REPO, item['path'])
                    check_pin(path, item)
                    member_name = item.get('archive_member', entry['bundle']+'/'+path.name)
                    source_archive, source_members = (original, old_members) if member_name in old_members else (published, new_members)
                    require(member_name in source_members and path.read_bytes() == member_bytes(source_archive, source_members[member_name]),
                            'Readable copy differs from its archive: '+str(path))
                    if kind == 'readable_handoff_documents':
                        document_count += 1
                    else:
                        record_count += 1
        rows.append({'bundle': entry['bundle'], 'original_archive_sha256': entry['original_sha256'],
                     'publication_archive_sha256': public['sha256'], 'publication_archive_bytes': public['bytes'],
                     'exact_archive_copy': not derivative, 'derivative': derivative_record})
    pin(SUPPORT/'ARCHIVE_SHA256SUMS')
    actual_sums = checksum_rows((SUPPORT/'ARCHIVE_SHA256SUMS').read_bytes())
    require(actual_sums == expected_sums, 'Archive checksum list coverage or values differ')
    require({str(path.relative_to(SUPPORT)) for path in (SUPPORT/'archives').iterdir() if path.is_file()} == set(expected_sums),
            'Archive directory contains uncatalogued or missing regular files')
    require(sum(item['exact_archive_copy'] for item in rows) == 24, 'Wrong exact-copy count')
    require(sum(item['publication_archive_bytes'] for item in rows) == 109653944, 'Publication archive byte total')
    return {'completed_archives': 25, 'exact_archive_copies': 24, 'publication_derivatives': 1,
            'readable_handoff_documents': document_count, 'readable_metadata_records': record_count,
            'archive_checksums': len(actual_sums), 'archives': rows}


def verify_incomplete(catalog):
    record = catalog['preserved_incomplete_work']
    require(record['status'] == 'incomplete' and record['completed_handoff'] is False and
            record['pass_id'] == INCOMPLETE.name, 'Incomplete pass promoted to completed')
    for item in record['preserved_files']:
        check_pin(relative(REPO, item['path']), item)
        source = relative(INCOMPLETE, item['original_relative_path'])
        require(pin(source)['sha256'] == item['original_sha256'] == item['sha256'], 'Incomplete copy differs')
        require(item['status'] == 'incomplete_unexecuted_or_partially_executed', 'Incomplete file status')
    readme_path = relative(REPO, record['readme'])
    readme = readme_path.read_text().lower().replace('**','')
    for phrase in ('status: incomplete', 'not a completed', 'not a runnable release',
                   'command-record suffix', 'archive/evidence coverage'):
        require(phrase in readme, 'Missing incomplete/readiness qualification: '+phrase)
    for key in ('acquisition_records','source_cache_index'):
        check_pin(relative(REPO, record[key]['path']), record[key])
    acquisition = parsed(relative(REPO, record['acquisition_records']['path']).read_bytes())
    source_index = parsed(relative(REPO, record['source_cache_index']['path']).read_bytes())
    require(acquisition['status'] == 'incomplete' and source_index['raw_source_files_republished'] is False,
            'Source acquisition qualification')
    expected_files = {relative(REPO,item['path']) for item in record['preserved_files']}
    expected_files.update({readme_path,relative(REPO,record['acquisition_records']['path']),
                           relative(REPO,record['source_cache_index']['path'])})
    require({path for path in readme_path.parent.rglob('*') if path.is_file()} == expected_files,
            'Incomplete subtree has undeclared or missing files')
    for path in readme_path.parent.rglob('*'):
        require(not path.is_symlink(), 'Incomplete subtree contains symlink')
        require(not path.is_file() or path.suffix.lower() not in {'.pdf','.html','.htm','.png','.jpg','.jpeg','.tsv','.body','.headers'},
                'Unexpected raw source payload in incomplete subtree')
    return {'preserved_files': len(record['preserved_files']), 'source_cache_metadata_records':len(source_index['records']),
            'raw_source_files_opened_by_verifier':0, 'completed_handoff':False,
            'unfinished_validation_gaps_preserved':True}


class DecodePolicy(ast.NodeTransformer):
    def visit_Call(self, node):
        node = self.generic_visit(node)
        if (isinstance(node.func, ast.Attribute) and node.func.attr == 'fromhex' and
                isinstance(node.func.value, ast.Name) and node.func.value.id == 'bytes' and
                len(node.args) == 1 and isinstance(node.args[0], ast.Constant) and not node.keywords):
            return ast.copy_location(ast.Constant(bytes.fromhex(node.args[0].value)), node)
        if (isinstance(node.func, ast.Attribute) and node.func.attr == 'decode' and len(node.args) == 1 and
                isinstance(node.args[0], ast.Constant) and node.args[0].value in ('ascii','utf-8') and not node.keywords):
            inner = node.func.value
            if isinstance(inner, ast.Constant) and isinstance(inner.value, bytes):
                return ast.copy_location(ast.Constant(inner.value.decode(node.args[0].value)), node)
        return node


def command_records(data):
    require(data.endswith(b'\n'), 'Command journal has incomplete trailing record')
    result = [parsed(line) for line in data.splitlines()]
    for record in result:
        require(isinstance(record, dict) and isinstance(record.get('command'), list) and
                record['command'] and all(isinstance(item,str) for item in record['command']) and
                type(record.get('exit_code')) is int and isinstance(record.get('cwd'),str) and
                re.fullmatch(r'[0-9a-f]{64}', record.get('log_sha256','')) is not None,
                'Command record schema')
    return result


def verify_additional_evidence():
    policy_path = SUPPORT/'publication/helper_policy_encoding.json'
    pin(policy_path)
    policy = parsed(policy_path.read_bytes())
    require(policy['scientific_source_data_results_or_equations_changed'] is False and
            policy['models_tests_or_source_acquisition_reexecuted'] is False, 'Helper projection qualification')
    projections = []
    for item in policy['files']:
        path = relative(REPO,item['published_path'])
        original = WORK/path.name
        require(pin(original)['sha256'] == item['original_source_sha256'] and
                pin(path)['sha256'] == item['published_source_sha256'], 'Helper source version binding')
        original_tree = ast.parse(original.read_text())
        published_tree = DecodePolicy().visit(ast.parse(path.read_text()))
        require(ast.dump(original_tree) == ast.dump(published_tree), 'Helper projection changed code beyond encoding')
        projections.append({'path':item['published_path'], 'original_sha256':item['original_source_sha256'],
                            'published_sha256':item['published_source_sha256'], 'equivalent_AST_after_literal_decoding':True})
    snapshot_path = SUPPORT/'publication/command_snapshot.jsonl'
    snapshot = snapshot_path.read_bytes()
    live = (WORK/'commands.jsonl').read_bytes()
    require(live.startswith(snapshot), 'Publication command snapshot is not exact byte prefix')
    snapshot_records, live_records = command_records(snapshot), command_records(live)
    retirement_root = SUPPORT/'publication/mhd_retirement'
    pin(retirement_root/'evidence_index.json')
    index = parsed((retirement_root/'evidence_index.json').read_bytes())
    require(index['focused_tests_passed'] == 7 and index['focused_test_runs'] == 1 and
            index['CLI_help_import_checks_passed'] == 3 and
            index['scientific_evolution_or_old_suites_replayed'] is False, 'Retirement evidence labels')
    for item in index['files']:
        check_pin(relative(REPO,item['path']),item)
    for record in index['commands']:
        require(record in live_records and record['exit_code'] == 0, 'Retirement command association')
        log = retirement_root/Path(record['log']).name
        require(pin(log)['sha256'] == record['log_sha256'], 'Retirement copied command log')
    require(len(index['commands']) == 4, 'Retirement command count')
    xml = ET.fromstring((retirement_root/'focused-results.xml').read_bytes())
    cases = list(xml.iter('testcase'))
    require(len(cases) == 7 and not list(xml.iter('failure')) and not list(xml.iter('error')), 'Retirement JUnit record')
    require(parsed((SUPPORT/'publication/assembly_record.json').read_bytes()) ==
            parsed((WORK/'supporting-document-assembly.json').read_bytes()), 'Assembly record copied bytes/contents')
    require((SUPPORT/'publication/assembly_record.json').read_bytes() ==
            (WORK/'supporting-document-assembly.json').read_bytes(), 'Assembly record is not exact copied bytes')
    copy_failure = SUPPORT/'publication/evidence-copy-initial-failure.log'
    original_failure = WORK/'preserve-maintenance-publication-evidence.log'
    require(pin(copy_failure) == pin(original_failure), 'Initial copy failure evidence differs')
    references_root = SUPPORT/'references/yacora_rotational'
    publication_record = references_root/'publication_record.json'
    pin(publication_record)
    references = parsed(publication_record.read_bytes())
    require(references['schema'] == 'mitl-yacora-literature-publication-v1' and
            references['full_papers_or_solver_data_republished'] is False, 'Literature publication scope')
    require(len(references['files']) == 9, 'Literature copied file count')
    for item in references['files']:
        check_pin(relative(references_root,item['name']),item)
    literature_verification = parsed((references_root/'verification.json').read_bytes())
    literature_summary = parsed((references_root/'summary.json').read_bytes())
    for item in references['files']:
        if item['name'] in {'summary.json','requests.jsonl','summarize_findings.py','read_official_pages.py'}:
            require(literature_verification['artifact_sha256'][item['name']] == item['sha256'],
                    'Literature verification copied-file association')
    require(literature_verification['status'] == 'passed' and
            literature_verification['primary_citations_checked'] == len(literature_summary['primary_citations']) == 3 and
            literature_verification['full_copyrighted_PDFs_downloaded_or_committed'] is False,
            'Literature evidence retained qualification')
    amendment = parsed((references_root/'helper-retention-amendment.json').read_bytes())
    require(amendment['schema'] == 'powerflow-yacora-source-helper-retention-amendment-v1' and
            amendment['status'] == 'metadata_lineage_clarified', 'Literature retention amendment identity')
    for key in ('preserved_verification','pre_retention_helper','current_helper','focused_difference','summary'):
        item = amendment[key]
        check_pin(relative(references_root,item['path']),item)
    require(amendment['preserved_verification']['artifact_helper_sha256'] ==
            amendment['current_helper']['sha256'] == literature_verification['artifact_sha256']['read_official_pages.py'],
            'Literature current/preserved verification helper binding')
    require((references_root/'verification.before_helper_lineage_amendment.json').read_bytes() ==
            (references_root/'verification.json').read_bytes(), 'Literature preserved current verification copy')
    require(amendment['earlier_verification_observation']['original_bytes_available_now'] is False and
            'reconstruction' in amendment['pre_retention_helper']['preservation_method'].lower() and
            'not an original' in amendment['pre_retention_helper']['preservation_method'].lower() and
            amendment['summary']['bibliographic_findings_changed'] is False,
            'Literature unavailable/reconstructed lineage qualifications')
    require(all(amendment[key] == 0 for key in ('new_source_retrievals','reader_or_aggregation_reruns','network_calls','scientific_calls')),
            'Literature metadata-only amendment')
    require({path.name for path in references_root.iterdir() if path.is_file()} ==
            {item['name'] for item in references['files']}|{'publication_record.json'},
            'Literature publication has undeclared files')
    return {'helper_policy_projections':projections,'retirement_evidence_files':len(index['files']),
            'recorded_retirement_test_cases':len(cases),'recorded_retirement_commands':len(index['commands']),
            'literature_metadata_copies':len(references['files']),'literature_recorded_primary_citations':3,
            'command_snapshot':{**pin(snapshot_path),'complete_records':len(snapshot_records),
                                'live_complete_records_at_verification':len(live_records),'exact_byte_prefix':True}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--catalog-sha256',required=True)
    args = parser.parse_args()
    destination = OUT/'report.json'
    require(not destination.exists(), 'Refusing to overwrite existing publication verification')
    result = {'schema':'mitl-independent-supporting-publication-verification-v1',
              'verified_utc':datetime.now(timezone.utc).isoformat(), 'status':'failed'}
    try:
        catalog_path = SUPPORT/'catalog.json'
        require(pin(catalog_path)['sha256'] == args.catalog_sha256, 'Trusted final catalog SHA mismatch')
        catalog = parsed(catalog_path.read_bytes())
        require(catalog['schema'] == 'mitl-supporting-document-publication-catalog-v1' and
                catalog['scientific_models_tests_or_recovery_reexecuted_by_assembly'] is False, 'Catalog labels')
        pin(SUPPORT/'publication/assembly_record.json')
        pin(WORK/'supporting-document-assembly.json')
        assembly = parsed((SUPPORT/'publication/assembly_record.json').read_bytes())
        require(assembly['catalog']['sha256'] == args.catalog_sha256 and
                assembly['catalog']['bytes'] == pin(catalog_path)['bytes'] and
                assembly['incomplete_source_pass_completed'] is False, 'Assembly catalog binding')
        result['archive_verification'] = verify_archives(catalog)
        result['incomplete_work_verification'] = verify_incomplete(catalog)
        result['publication_evidence'] = verify_additional_evidence()
        documents = sorted({REPO/'README.md',REPO/'HANDOFF.md',REPO/'docs/external_plasma_yacora_comparison.md',
                            *SUPPORT.rglob('*.md')})
        for path in documents:
            pin(path)
        links = verify_documents(REPO,documents,allowed_pending_outputs={SELF_REPORT})
        result['relative_links'] = {'documents_checked':len(documents),'relative_links_checked':len(links),
                                   'links':links,'planned_self_report_links':sum(item['planned_output'] for item in links)}
        require(catalog['verification_record'] == str(SELF_REPORT.relative_to(REPO)), 'Unexpected planned report target')
        result['status'] = 'passed_for_frozen_publication_inputs'
        result['all_input_checks_passed'] = True
        result['pending_output_closure'] = {'path':str(SELF_REPORT.relative_to(REPO)),
                                          'requirement':'Copy this exact real report, then verify its SHA and target existence.',
                                          'passing_placeholder_used':False}
        result['limits'] = ['Existing historical verification records were preserved, not scientifically rerun.',
                            'This verifies publication bytes, associations and local links, not physical truth.',
                            'External URLs were not fetched; source-cache contents were not opened.',
                            'The copied self-report link requires the separate exact-byte output closure.']
        result['scientific_models_or_tests_executed_by_verifier'] = False
        result['source_acquisition_or_external_writes_by_verifier'] = False
        result['repository_or_status_writes_by_verifier'] = False
    except Exception as error:
        result['all_input_checks_passed'] = False
        result['error'] = {'type':type(error).__name__,'message':str(error),'traceback':traceback.format_exc()}
    for path in (Path(__file__),OUT/'markdown_links.py',OUT/'plan.json'):
        pin(path)
    result['input_sha256'] = {str(path):record for path,record in sorted(PINS.items())}
    destination.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'report':str(destination),'sha256':digest(destination.read_bytes()),
                      'error':result.get('error')},indent=2))
    return 0 if result['all_input_checks_passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
