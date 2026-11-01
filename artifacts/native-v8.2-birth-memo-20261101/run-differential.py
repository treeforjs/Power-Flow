"""Bounded actual-archive comparison. SMT occupancy invalidates timing conclusions,
but never suppresses physical parity/refusal validation."""
import sys
from tooling import *
resource_guard()
label=sys.argv[1] if len(sys.argv)>1 else 'differential'
n=int(sys.argv[2]) if len(sys.argv)>2 else 10000
steps=int(sys.argv[3]) if len(sys.argv)>3 else 8
out=OUT/label; r=Runner(out,'differential')
build=json.loads((OUT/'build-differential/build-result.json').read_text())
inputs=sealed_inputs()
def stats():
    rows={}
    for line in Path('/proc/stat').read_text().splitlines():
        t=line.split()
        if t and t[0] in ['cpu6','cpu7','cpu14','cpu15']: rows[t[0]]=list(map(int,t[1:9]))
    return rows
def processes():
    return subprocess.check_output(['ps','-eLo','pid,tid,psr,pcpu,stat,comm,args','--sort=-pcpu'],text=True)
before=stats(); (out/'processes-before.txt').write_text(processes())
reports={}; receipts=[]
for side in ['baseline','candidate','candidate','baseline']:
    run=side+'-'+str(sum(x.startswith(side) for x in reports))
    r.run([P/'build'/(side+'-birth-memo-diff'),out/run,n,steps],run+'.log',timeout=240,rlimit=3*1024**3)
    reports[run]=json.loads((out/run/'report.json').read_text())
    receipts.append({'run':run,'physical_sha256':sha256(out/run/'physical.bin'),
                     'physical_bytes':(out/run/'physical.bin').stat().st_size,
                     'birth_prepare_seconds':reports[run]['birth_prepare_seconds'],
                     'full_step_seconds':reports[run]['full_step_seconds']})
after=stats(); (out/'processes-after.txt').write_text(processes())
busy={}
for cpu,first in before.items():
    delta=[b-a for a,b in zip(first,after[cpu])]; total=sum(delta)
    busy[cpu]=None if total<=0 else (total-delta[3]-delta[4])/total
active=next(iter(os.sched_getaffinity(0)))
sibling_text=Path(f'/sys/devices/system/cpu/cpu{active}/topology/thread_siblings_list').read_text().strip()
siblings=set()
for part in sibling_text.split(','):
    if '-' in part:
        lo,hi=map(int,part.split('-'));siblings.update(range(lo,hi+1))
    else:siblings.add(int(part))
siblings.discard(active)
quiet=all(busy.get('cpu'+str(cpu),1) is not None and busy.get('cpu'+str(cpu),1)<=.01 for cpu in siblings)
clean=quiet and active in (6,7) and os.environ.get('PF_TIMING_QUALIFIED')=='1'
ref=reports['baseline-0']; differences=[]; seals={}
for key,report in reports.items():
    for field in ['physical_step_sha256','refusals','bulk_carriers','steps','births','sunk','newborn_sink_steps','final_carriers','checks','failures']:
        if report[field]!=ref[field]:differences.append(key+': '+field)
    if receipts[list(reports).index(key)]['physical_sha256']!=receipts[0]['physical_sha256']:
        differences.append(key+': physical.bin')
    seals[key]={field:sum(a[field]!=b[field] for a,b in zip(ref['release_provenance_seals'],report['release_provenance_seals']))
                for field in ref['release_provenance_seals'][0]}
need(sealed_inputs()==inputs,'dependencies changed during differential')
need(sha256(P/'build/libpic_native_v8.a')==build['candidate_library'],'candidate archive changed')
need(sha256(RT/'native-full-time-candidate008-auditfix/build/libpic_native_v8.a')==build['baseline_library'],'baseline archive changed')
result={'passed':not differences,'differences':differences,'build':build,'runs':receipts,
        'refusal_controls':ref['refusals'],'source_identity':{k:v['source_identity'] for k,v in reports.items()},
        'source_identity_bound_seal_differences':seals,'cpu':active,'thread_siblings_list':sibling_text,
        'cpu_busy_fraction':busy,'sibling_quiet':quiet,'timing_uncontaminated':clean,
        'timing_note':'Timing is qualified only for an explicitly designated owner-core timing job when its SMT sibling stays at <=1% busy; physical parity remains valid regardless. Borrowed functional qualification is not timing evidence. This manufactured zero-field birth/sink fixture is not a production forecast.',
        'commands':r.records}
(out/'differential-result.json').write_text(json.dumps(result,indent=2)+'\n')
need(not differences,'actual archive differential differs: '+str(differences))
print('actual-archive physical bytes/refusals match; timing uncontaminated:',clean)
