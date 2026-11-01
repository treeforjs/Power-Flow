#!/usr/bin/python -I
"""Final-build parity and interleaved CPU timing in a normal one-core SimQueue job.

The reference is qualified v8.0.1/c2, never the older optimization study binary.
Timing evidence excludes a sample if the actual SMT sibling was busy >1% (or
>0.03 CPU seconds for very short samples). It never stops another worker.
"""
import array,hashlib,importlib.util,json,os,resource,statistics,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
BASE=ROOT.parent/'pic-water-common-commit-v8.0.1-20261031'
BASE_OUT=Path('/home/treefor/projects/.workspace/outputs/pic-water-common-commit-v8.0.1-20261031')
MESH=ROOT.parent/'resolution-series-prep-20261031/meshes/mesh-200um.json'

def module(name,file):
    spec=importlib.util.spec_from_file_location(name,file);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
m=module('candidate_run',ROOT/'run_v8.py')
compare=module('compare_v8',ROOT/'tools/compare_v8.py')

def stats(cpus):
    rows={}
    for line in Path('/proc/stat').read_text().splitlines():
        key,*fields=line.split()
        if key in {'cpu'+str(c) for c in cpus}:rows[key]=list(map(int,fields[:8]))
    return rows

def workers(cpus):
    # Evidence only. A broad allowed affinity is not by itself running load.
    result=[]
    for proc in Path('/proc').glob('[0-9]*'):
        try:
            for task in (proc/'task').iterdir():
                raw=(task/'stat').read_text();fields=raw[raw.rfind(')')+2:].split()
                allowed=os.sched_getaffinity(int(task.name));last_cpu=int(fields[36])
                if last_cpu in cpus or (len(allowed)==1 and allowed&set(cpus)):
                    result.append({'pid':int(proc.name),'tid':int(task.name),'state':fields[0],
                       'last_cpu':last_cpu,'allowed':sorted(allowed),'utime_ticks':int(fields[11]),'stime_ticks':int(fields[12]),
                       'cmd':(proc/'cmdline').read_bytes().replace(b'\0',b' ').decode(errors='replace')[:400]})
        except (FileNotFoundError,ProcessLookupError,PermissionError):pass
    return result

def memory_snapshot():
    info={k:int(v.split()[0])*1024 for k,v in (line.split(':',1) for line in Path('/proc/meminfo').read_text().splitlines())}
    tasks=[]
    for proc in Path('/proc').glob('[0-9]*'):
        try:
            status=dict(line.split(':',1) for line in (proc/'status').read_text().splitlines())
            rss=int(status.get('VmRSS','0 kB').split()[0])*1024
            if rss:tasks.append({'pid':int(proc.name),'name':status.get('Name','').strip(),'rss_bytes':rss})
        except (FileNotFoundError,PermissionError,ProcessLookupError):pass
    return {'mem_available_bytes':info['MemAvailable'],'mem_total_bytes':info['MemTotal'],
            'swap_total_bytes':info['SwapTotal'],'swap_free_bytes':info['SwapFree'],
            'largest_processes':sorted(tasks,key=lambda p:p['rss_bytes'],reverse=True)[:12]}

def table(path,kind):
    cols=json.loads((path/(kind+'-schema.json')).read_text())['columns']
    data=array.array('d');data.frombytes((path/(kind+'.bin')).read_bytes())
    if len(data)%len(cols):raise RuntimeError('partial timing row')
    return [dict(zip(cols,data[i:i+len(cols)])) for i in range(0,len(data),len(cols))]

def main(label,tag,mode='parity'):
    if mode not in ('parity','timing'):raise RuntimeError('parity or timing required')
    if not os.getenv('SIMQ_JOB_ID') or not m.inside_frozen():raise RuntimeError('frozen SimQueue execution required')
    adm=m.admission('candidate-'+mode)
    cpu=adm['cpu_affinity'][0]
    if mode=='timing' and (cpu not in (6,7) or os.getenv('SIMQ_BORROW')=='1'):raise RuntimeError('timing requires an owner Power-Flow core')
    if cpu not in (6,7) and os.getenv('SIMQ_BORROW')!='1':raise RuntimeError('Power-Flow or explicitly borrowed core required')
    sibling_list=Path('/sys/devices/system/cpu/cpu%d/topology/thread_siblings_list'%cpu).read_text().strip()
    siblings=[int(c) for c in sibling_list.split(',') if int(c)!=cpu]
    output=m.OUT/tag;output.mkdir(parents=True,exist_ok=False)
    mesh_sha=m.sha(MESH);own_before=m.own_sources();pins=m.sealed_pins()
    br=json.loads((m.OUT/('build-'+label)/'build-receipt.json').read_text())
    baseline_br=json.loads((BASE_OUT/'build-c2/build-receipt.json').read_text())
    if baseline_br['binary_sha256']!='3e5d52b2eb93d15fdfbaa79a7a3541911fa987905a9bbb55c14f9411d7e86072':raise RuntimeError('qualified baseline binary pin differs')
    results={'schema':'coordinator-cpu-comparison-v1','mode':mode,'admission':adm,'smt_siblings':siblings,
       'baseline_binary_sha256':baseline_br['binary_sha256'],'candidate_binary_sha256':br['binary_sha256'],
       'mesh_sha256':mesh_sha,'runs':[],'comparisons':[],'memory_admissions':[],'timing_policy':'1% sibling busy time or 0.03 CPU s absolute floor; 3 paired repetitions; no timing acceptance from a contaminated sample'}
    def save():
        (output/'comparison.json').write_text(json.dumps(results,indent=2)+'\n')
    if mode=='timing':
        pre=stats(siblings);time.sleep(2);post=stats(siblings)
        busy={str(c):sum(post['cpu'+str(c)][i]-pre['cpu'+str(c)][i] for i in (0,1,2,5,6,7))/os.sysconf('SC_CLK_TCK') for c in siblings}
        results['preflight_sibling_busy_cpu_s']=busy
        results['preflight_sibling_tasks']=workers(siblings)
        if any(v>.03 for v in busy.values()):
            results['passed']=False
            results['timing_not_qualified']='SMT sibling already busy at preflight; no numerical timing run started'
            save()
            raise RuntimeError(results['timing_not_qualified'])
    def run(name,which,options):
        memory=memory_snapshot()
        cap=int(os.environ['SIMQ_MEM'])
        reserve=4*1024**3
        admitted=0<cap<=4*1024**3 and memory['mem_available_bytes']>=cap+reserve
        results['memory_admissions'].append({'name':name,'snapshot':memory,'scheduler_cap_bytes':cap,
           'physical_reserve_bytes':reserve,'required_mem_available_bytes':cap+reserve,'admitted':admitted})
        save()
        if not admitted:raise RuntimeError('physical memory admission refused before '+name+'; need 4GiB job cap and cap+4GiB MemAvailable')
        root,build=(BASE,'c2') if which=='baseline' else (ROOT,label)
        out=output/name
        command=[sys.executable,'-I',root/'run_v8.py','run',build,'run',out,'--hall-owner','owner',*options]
        before=stats(siblings);active_before=workers(siblings);wall=time.monotonic();ru=resource.getrusage(resource.RUSAGE_CHILDREN)
        p=subprocess.run(list(map(str,command)),cwd=root,capture_output=True,text=True)
        wall=time.monotonic()-wall;after=stats(siblings);ru_after=resource.getrusage(resource.RUSAGE_CHILDREN)
        if not out.is_dir():raise RuntimeError(p.stdout+p.stderr)
        (out/'wrapper.log').write_text(p.stdout+p.stderr)
        receipt=json.loads((out/'run-receipt.json').read_text())
        if not all(receipt.get(k) is True for k in ('sealed_unchanged_after_run','own_sources_unchanged_after_run','binary_unchanged_after_run')):raise RuntimeError('run integrity failed')
        summary=json.loads((out/'run-summary.json').read_text())
        hz=os.sysconf('SC_CLK_TCK');busy={}
        for c in siblings:
            a=before['cpu'+str(c)];b=after['cpu'+str(c)];delta=[y-x for x,y in zip(a,b)]
            # user/nice/system/irq/softirq/steal; idle+iowait are not executing sibling work.
            busy[str(c)]=sum(delta[i] for i in (0,1,2,5,6,7))/hz
        quiet=all(v<=max(.03,.01*wall) for v in busy.values())
        rows=table(out,'ledger');timings=table(out,'coordinator-timings')
        record={'name':name,'which':which,'directory':str(out),'returncode':p.returncode,'stop':summary['stop'],
           'steps':summary['steps_completed'],'final_time_s':summary['final_time_s'],'wall_s':wall,
           'child_cpu_s':ru_after.ru_utime+ru_after.ru_stime-ru.ru_utime-ru.ru_stime,
           'memory_after':memory_snapshot(),'sibling_busy_cpu_s':busy,'quiet_sibling':quiet,'sibling_tasks_before':active_before,'sibling_tasks_after':workers(siblings),
           'ledger_timing_rows':[{k:v for k,v in r.items() if k.startswith('t_') or k in ('step','rss_kib','carriers_electron','carriers_ion','carriers_h2o')} for r in rows],
           'coordinator_timing_rows':timings}
        results['runs'].append(record);save();print(name,p.returncode,summary['steps_completed'],summary['stop'],'quiet='+str(quiet),flush=True)
        return record
    def pair(name,a,b):
        if (a['returncode'],a['steps'],a['stop'],a['final_time_s'])!=(b['returncode'],b['steps'],b['stop'],b['final_time_s']):raise RuntimeError(name+': accepted/refused trajectory differs')
        if a['steps']<=0:raise RuntimeError(name+': no accepted steps for parity')
        result=compare.identity(a['directory'],b['directory'],a['steps'])
        result['name']=name;results['comparisons'].append(result);save()
        if not result['bit_identical']:raise RuntimeError(name+': deterministic parity failed')
    if mode=='parity':
        checks=compare.checks(BASE_OUT/'c2-checks',m.OUT/(label+'-checks'))
        results['checks']=checks
        if not checks['reference_cases_first_same_names_and_outcomes'] or checks['evidence_differences'] or checks['passed_candidate']!=53:raise RuntimeError('53-check evidence differs')
        fixtures=[('qual20',['--fixture','qual','--steps','20','--frame-every','10']),
          ('dense-adaptive',['--fixture','dense','--blob-weight','8e10','--dt-mode','adaptive','--sync-every','5e-14','--window','5e-13','--progress-every','50']),
          ('load100k',['--fixture','load','--load-carriers','100000','--steps','3','--mesh',str(MESH),'--progress-every','1']),
          ('load1m',['--fixture','load','--load-carriers','1000000','--steps','2','--mesh',str(MESH),'--progress-every','1']),
          ('pilot2000',['--fixture','pilot','--steps','2000','--frame-every','1000','--progress-every','1000'])]
        for name,args in fixtures:
            a=run(name+'-baseline','baseline',args);b=run(name+'-candidate','candidate',args);pair(name,a,b)
    else:
        args=['--fixture','load','--load-carriers','1000000','--steps','3','--mesh',str(MESH),'--progress-every','1']
        pairs=[]
        for i,order in enumerate((('baseline','candidate'),('candidate','baseline'),('baseline','candidate'))):
            runs={which:run('rep%d-%s'%(i+1,which),which,args) for which in order}
            pair('rep%d'%(i+1),runs['baseline'],runs['candidate']);pairs.append(runs)
        admitted=[r for r in pairs if all(x['quiet_sibling'] and x['returncode']==0 and x['steps']==3 for x in r.values())]
        results['timing_pairs_admitted']=len(admitted)
        if len(admitted)==3:
            vals={which:[statistics.median(r[which]['coordinator_timing_rows'][j]['t_step_total_ms'] for j in (1,2)) for r in admitted] for which in ('baseline','candidate')}
            med={which:statistics.median(v) for which,v in vals.items()}
            results['qualified_timing']={'scope':'manufactured 1e6-carrier 200um load fixture; steps 2 and 3, one CPU core, all CPU backends',
              'per_repetition_median_step_ms':vals,'median_step_ms':med,'baseline_over_candidate':med['baseline']/med['candidate'],
              'candidate_over_baseline':med['candidate']/med['baseline']}
        else:results['timing_not_qualified']='SMT sibling contamination or run refusal; repeat when quiet, no speed claim'
    results['inputs_unchanged']=m.own_sources()==own_before and m.sealed_pins()==pins and m.sha(MESH)==mesh_sha
    results['passed']=results['inputs_unchanged'] and all(c['bit_identical'] for c in results['comparisons'])
    save()
    if not results['passed']:raise RuntimeError('final integrity/parity failure')

if __name__=='__main__':main(*sys.argv[1:])
