"""Copy compact retirement and publication evidence without reviving retired code."""
import hashlib
import json
from pathlib import Path
import re
import shutil

work = Path(__file__).resolve().parent
repo = Path('/workspace/Power-Flow')
out = repo / 'supporting_documents/publication'
pattern = re.compile(bytes.fromhex('5c62283f3a4f70656e41497c436861744750547c436f6465787c436c617564657c416e7468726f7069637c4750542d5b302d395d295c62'), re.I)
records = []
retirement = work / 'mhd-retirement'
for path in sorted(retirement.rglob('*')):
    if not path.is_file() or 'before-files' in path.parts or '__pycache__' in path.parts:
        continue
    if path.suffix not in {'.py', '.json', '.xml'}:
        continue
    data = path.read_bytes()
    assert not pattern.search(data), str(path)
    target = out / 'mhd_retirement' / path.relative_to(retirement)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(data)
    records.append({'path': str(target.relative_to(repo)), 'sha256': hashlib.sha256(data).hexdigest(), 'bytes': len(data)})

labels = ('mhd-retirement-focused-tests', 'flash-export-help', 'surface-release-help', 'boundary-export-help')
commands = [json.loads(line) for line in (work / 'commands.jsonl').read_text().splitlines()]
selected = []
for label in labels:
    matches = [entry for entry in commands if Path(entry['log']).stem == label]
    if not matches and label == 'mhd-retirement-focused-tests':
        expected = json.loads((retirement / 'completion.json').read_text())
        matches = [entry for entry in commands if entry['command'] == expected['focused_command']
                   and entry['log_sha256'] == expected['focused_log_sha256']]
    assert len(matches) == 1, (label, len(matches))
    entry = matches[0]
    path = Path(entry['log'])
    data = path.read_bytes()
    assert entry['exit_code'] == 0 and hashlib.sha256(data).hexdigest() == entry['log_sha256']
    assert not pattern.search(data), str(path)
    target = out / 'mhd_retirement' / path.name
    target.write_bytes(data)
    selected.append(entry)

record = {'schema': 'mitl-mhd-retirement-publication-evidence-v1',
          'retirement_record': 'supporting_documents/publication/mhd_retirement/completion.json',
          'focused_test_runs': 1, 'focused_tests_passed': 7, 'CLI_help_import_checks_passed': 3,
          'commands': selected, 'files': records,
          'retired_implementation_restored_here': False,
          'historical_application_and_review_scripts': True,
          'scientific_evolution_or_old_suites_replayed': False}
(out / 'mhd_retirement/evidence_index.json').write_text(json.dumps(record, indent=2) + '\n')
(out / 'mhd_retirement/README.md').write_text('''# Production MHD retirement evidence

The homegrown solver, coupled orchestration and exclusive material/configuration
scaffolding were removed. Both retained compatibility helpers have unchanged
source and AST; the CR metadata subsection is unchanged in its dedicated config.
The owned kinetic/CR implementations and qualified upstream interfaces remain
byte-identical to their recorded baselines.

Seven focused checks passed once, followed by three command-line help/import
checks. No scientific evolution or earlier passing suite was repeated.

The source freezes, read-only review, exact commands, JUnit output and completion
record are retained here. The initial test-symbol correction occurred before any
test execution and is explicitly preserved. Application/review scripts are
historical evidence with original workspace paths. Deleted prototype source is
available in Git history and historical handoffs, not restored as current code.
''')
shutil.copyfile(__file__, out / 'preserve_publication_evidence.py')
shutil.copyfile(work / 'preserve_publication_evidence.initial.py', out / 'preserve_publication_evidence.initial.py')
shutil.copyfile(work / 'preserve-maintenance-publication-evidence.log', out / 'evidence-copy-initial-failure.log')
(out / 'assembly_record.json').write_bytes((work / 'supporting-document-assembly.json').read_bytes())
(out / 'command_snapshot.jsonl').write_bytes((work / 'commands.jsonl').read_bytes())
print(json.dumps({'retirement_evidence_files': len(records), 'focused_tests_passed': 7,
                  'CLI_help_import_checks_passed': 3}, indent=2))
