"""Check the real publication report/source copies without repeating archive work."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

WORK = Path('/workspace/research_work/mitl-github-push-20261020T0601')
LOCAL = WORK/'publication-verification'
REPO = Path('/workspace/Power-Flow')
PUBLIC = REPO/'supporting_documents/publication'


def pin(path):
    data = path.read_bytes()
    return {'sha256':hashlib.sha256(data).hexdigest(),'bytes':len(data)}


def require(value,message):
    if not value:
        raise ValueError(message)


def main():
    output = LOCAL/'output-closure.json'
    require(not output.exists(),'Refusing to overwrite output closure')
    final = LOCAL/'report.final.json'
    require(pin(final)['sha256'] == '1477c000600428878bece2056928be8a28b0a7f4b9272ae9f7942799eb9c7d7d',
            'Unexpected final input report')
    report = json.loads(final.read_text())
    require(report['all_input_checks_passed'] is True and report['status'] == 'passed_for_frozen_publication_inputs',
            'Input verification did not pass')
    copies = []
    for local_name,public_name in (('report.final.json','verification.json'),
                                  ('verify_publication.py','verify_publication.py'),
                                  ('markdown_links.py','markdown_links.py'),('plan.json','plan.json'),
                                  ('verify_publication.initial.py','verify_publication.initial.py'),
                                  ('report.initial-failure.json','report.initial-failure.json')):
        source,target = LOCAL/local_name,PUBLIC/public_name
        require(target.is_file() and not target.is_symlink() and source.read_bytes() == target.read_bytes(),
                'Publication copy mismatch: '+public_name)
        copies.append({'source':str(source),'published':str(target.relative_to(REPO)),**pin(target)})
    pending = [item for item in report['relative_links']['links'] if item['planned_output']]
    require(len(pending) == 2,'Unexpected planned output-link count')
    for item in pending:
        require(item['resolved_target'] == 'supporting_documents/publication/verification.json' and not item['fragment'],
                'Unexpected planned output target')
        document = REPO/item['document']
        require(pin(document) == report['input_sha256'][str(document)],'Link-containing document changed')
        require((REPO/item['resolved_target']).is_file(),'Planned output target remains missing')
    journal = (WORK/'commands.jsonl').read_bytes()
    require(journal.endswith(b'\n'),'Incomplete command journal')
    entries = [json.loads(line) for line in journal.splitlines()]
    commands = []
    for label,expected_exit in (('publication-artifact-readback',1),('publication-metadata-link-continuation',0)):
        matches = [row for row in entries if Path(row.get('log','')).stem == label]
        require(len(matches) == 1 and matches[0]['exit_code'] == expected_exit,'Verification command lineage mismatch')
        require(pin(Path(matches[0]['log']))['sha256'] == matches[0]['log_sha256'],'Verification command log mismatch')
        commands.append(matches[0])
    result = {'schema':'mitl-supporting-publication-output-closure-v1','status':'passed',
              'completed_utc':datetime.now(timezone.utc).isoformat(),'all_checks_passed':True,
              'input_verification_report_sha256':pin(final)['sha256'],'exact_copies':copies,
              'previously_pending_relative_links_resolved':pending,'pending_output_links_remaining':0,
              'verification_commands':commands,'source_sha256':{Path(__file__).name:pin(Path(__file__))['sha256']},
              'full_archive_readbacks':1,'focused_metadata_continuations':1,'focused_output_closures':1,
              'archives_scientific_models_or_tests_reexecuted_in_closure':False,
              'repository_status_source_cache_or_external_writes_by_closure':False}
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'report':str(output),**pin(output)},indent=2))


if __name__ == '__main__':
    main()
