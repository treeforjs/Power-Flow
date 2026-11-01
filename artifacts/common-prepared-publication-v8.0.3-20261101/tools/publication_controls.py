#!/usr/bin/python -I
from pathlib import Path
import hashlib,importlib.util,json,subprocess,sys
R=Path(__file__).resolve().parents[1]
s=importlib.util.spec_from_file_location('runner',R/'run_v8.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
label,tag,expect=sys.argv[1:]
assert expect in ('preimage','refusal') and m.inside_frozen()
adm=m.admission('publication-controls')
def run(*args):subprocess.run([sys.executable,'-I',str(R/'run_v8.py'),*args],check=True)
run('build',label)
run('run',label,'publication-guard',tag+'-allocation')
run('run',label,'unchanged-source-replacement',tag+'-replacement')
p=m.OUT/(tag+'-replacement/unchanged-source-replacement.json');r=json.loads(p.read_text())
if expect=='preimage':assert r['can_publish_after_replacement'] and r['published'] and r['actual_vs_recorded_mismatch'] and r['hall_advanced'] and r['foreign_thermal_unchanged']
else:assert not r['can_publish_after_replacement'] and not r['published'] and not r['actual_vs_recorded_mismatch'] and not r['hall_advanced'] and r['foreign_thermal_unchanged']
q=m.OUT/(tag+'-allocation/publication-guard.json');a=json.loads(q.read_text());assert a['passed'] and a['whole_publication_allocations']==0 and a['source_digests_match'] and a['arrival_digest_matches']
(m.OUT/(tag+'-result.json')).write_text(json.dumps({'passed':True,'admission':adm,'source_identity_replacement_expected':expect,'replacement_result_sha256':m.sha(p),'allocation_result_sha256':m.sha(q),'harness_sha256':m.sha(Path(__file__)),'physical_run':False},indent=2)+'\n')
