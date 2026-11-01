#!/usr/bin/python -I
"""Instrument coordinator and rebuilt source/Hall TUs; pinned native archive is unchanged."""
import importlib.util,json,os,sys
from pathlib import Path
R=Path(__file__).resolve().parents[1]
s=importlib.util.spec_from_file_location('runner',R/'run_v8.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
label=sys.argv[1]
assert m.inside_frozen() and int(os.getenv('SIMQ_MEM','0'))>=8*1024**3
adm=m.admission('publication-sanitized')
manifest={l.split(':')[0]:int(l.split()[1])*1024 for l in Path('/proc/meminfo').read_text().splitlines() if l.split()[1].isdigit()}
assert manifest['MemAvailable']>=int(os.environ['SIMQ_MEM'])+2*1024**3
m.COMMON=[x for x in m.COMMON if not x.startswith('-O')]+['-O1','-g','-fno-omit-frame-pointer','-fsanitize=address,undefined','-fno-sanitize-recover=all']
execute=m.execute
def instrument(name,cmd,*args,**kwargs):
 if name=='link':cmd=[*cmd,'-fsanitize=address,undefined']
 return execute(name,cmd,*args,**kwargs)
m.execute=instrument
if m.build(label,m.DEFAULT_BINDINGS):raise RuntimeError('instrumented coordinator build failed')
records=[]
for mode in ('checks','publication-guard','unchanged-source-replacement','retained-source-lifetime'):
 tag=label+'-'+mode
 rc=m.run(label,mode,tag,[]);records.append({'mode':mode,'returncode':rc})
 if rc:raise RuntimeError('instrumented '+mode+' failed')
p=m.OUT/(label+'-unchanged-source-replacement/unchanged-source-replacement.json');r=json.loads(p.read_text())
assert not r['can_publish_after_replacement'] and not r['published'] and not r['actual_vs_recorded_mismatch'] and not r['hall_advanced']
(m.OUT/(label+'-sanitized-result.json')).write_text(json.dumps({'passed':True,'admission':adm,'harness_sha256':m.sha(Path(__file__)),'records':records,'scope':'rebuilt coordinator/source/Hall translation units; existing pinned native archive and numerical libraries are not rebuilt here','memory_before':manifest},indent=2)+'\n')
