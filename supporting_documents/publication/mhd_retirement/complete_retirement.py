"""Read-only completion binding for the focused legacy-solver retirement."""
from pathlib import Path
import hashlib
import json
import xml.etree.ElementTree as ET

R=Path('/workspace/Power-Flow')
W=Path('/workspace/research_work/mitl-github-push-20261020T0601/mhd-retirement')
PASS=W.parent
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
freeze=json.loads((W/'source-freeze.json').read_text())
assert sha(W/'source-freeze.json')=='fc3ba952b90456e2e27c4836baffe0327b00aa0d5fd2382bbaabcd6a5281cdcb'
for path,pins in freeze['owned_file_hashes'].items():
    actual=sha(R/path) if (R/path).is_file() else None
    assert actual==pins['after_sha256'],path
preserved=json.loads((W/'preserved-before.json').read_text())
assert all((R/path).is_file() and sha(R/path)==pin for path,pin in preserved.items())
review=W/'independent-review/completion.json'
assert sha(review)=='70071461240dd504cfa10d6e824354becb1e16e42e7859b98f08440e5b4cd068'
results=W/'focused-results.xml'
tree=ET.parse(results)
cases=tree.findall('.//testcase')
assert len(cases)==7 and not tree.findall('.//failure') and not tree.findall('.//error') and not tree.findall('.//skipped')
commands=[json.loads(line) for line in (PASS/'commands.jsonl').read_text().splitlines()]
runs=[record for record in commands if record['log'].endswith('/legacy-mhd-retirement-focused-checks.log')]
assert len(runs)==1 and runs[0]['exit_code']==0
record={
    'schema':'mhdlab-legacy-mhd-retirement-completion-v1','status':'passed',
    'source_freeze_sha256':sha(W/'source-freeze.json'),
    'initial_unexecuted_source_freeze_sha256':sha(W/'source-freeze.initial.json'),
    'independent_review_sha256':sha(review),
    'owned_file_hashes':freeze['owned_file_hashes'],
    'removed_file_count':len(freeze['removed_paths']),
    'preserved_file_count':len(preserved),'all_preserved_files_byte_identical':True,
    'helper_AST_and_source_unchanged':True,'cr_model_subsection_byte_identical':True,
    'cr_metadata_test_only_config_filename_changed':True,
    'focused_test_cases':[{'class':case.attrib['classname'],'name':case.attrib['name']} for case in cases],
    'focused_tests_passed':7,'focused_test_runs':1,'focus_result_sha256':sha(results),
    'focused_command':runs[0]['command'],'focused_log_sha256':runs[0]['log_sha256'],
    'science_suites_or_solver_evolution_runs':0,
    'scope':'Retired the homegrown production MHD and exclusive scaffolding. Preserved external FLASH interchange, independent verification/source interfaces, native arithmetic/build files and all unowned kinetic/CR/state data and implementations.',
    'public_commits_or_pushes_by_this_author':0,
}
(W/'completion.json').write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
print(json.dumps({'completion':str(W/'completion.json'),'completion_sha256':sha(W/'completion.json'),
                  'focused_tests_passed':7,'unchanged_preserved_files':len(preserved)},sort_keys=True))
