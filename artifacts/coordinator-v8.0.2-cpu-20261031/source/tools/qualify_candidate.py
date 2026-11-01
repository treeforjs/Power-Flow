#!/usr/bin/python -I
"""Build and qualification, only inside a frozen one-core SimQueue job."""
import hashlib,importlib.util,json,os,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
spec=importlib.util.spec_from_file_location('runner',ROOT/'run_v8.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def call(argv):
    subprocess.run([str(x) for x in argv],check=True,cwd=ROOT)

def main(label):
    if not os.getenv('SIMQ_JOB_ID') or not m.inside_frozen():raise RuntimeError('frozen SimQueue execution required')
    admission=m.admission('hotpath-qualification')
    if admission['cpu_affinity'][0] not in (6,7) and os.getenv('SIMQ_BORROW')!='1':raise RuntimeError('Power-Flow or explicitly borrowed core required')
    out=m.OUT/('hotpaths-'+label);out.mkdir(parents=True,exist_ok=False)
    before=m.own_sources();pins=m.sealed_pins()
    exe=out/'test_cpu_hotpaths'
    cmd=['/usr/bin/c++',*m.COMMON,'-Wall','-Wextra','-Wpedantic','-Werror',*m.includes(),
         '-MD','-MF',out/'test_cpu_hotpaths.d',ROOT/'tests/test_cpu_hotpaths.cpp','-o',exe]
    records=[m.execute('compile-hotpaths',cmd,out,600)]
    if records[-1]['returncode']:raise RuntimeError(Path(records[-1]['log']).read_text())
    binary=m.sha(exe)
    records.append(m.execute('hotpaths',[exe,m.DATA],out,120))
    unchanged=m.own_sources()==before and m.sealed_pins()==pins and m.sha(exe)==binary
    result={'schema':'coordinator-cpu-hotpath-qualification-v1','admission':admission,'commands':records,
       'binary_sha256':binary,'source_sha256':before,'sealed':pins,'inputs_unchanged':unchanged,
       'passed':unchanged and all(c['returncode']==0 for c in records)}
    (out/'qualification.json').write_text(json.dumps(result,indent=2)+'\n')
    if not result['passed']:raise RuntimeError('hotpath checks failed: '+records[-1]['log'])
    print(Path(records[-1]['log']).read_text(),flush=True)
    call([sys.executable,'-I',ROOT/'run_v8.py','build',label])
    call([sys.executable,'-I',ROOT/'run_v8.py','run',label,'checks',label+'-checks'])
    call([sys.executable,'-I',ROOT/'tools/audit_regressions.py',label,label+'-continuation','replay'])
    call([sys.executable,'-I',ROOT/'tools/audit_run_integrity.py',m.OUT/(label+'-run-integrity-results.json')])

if __name__=='__main__':main(*sys.argv[1:])
