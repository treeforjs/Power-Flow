#!/usr/bin/env python3
"""Bounded static controls. Invoke under the documented CPU7 exclusive lease."""
import hashlib,json,os,pathlib,subprocess,sys,time

PREFIX=pathlib.Path('/home/treefor/projects/runtimes/pic-electric-exterior-20261029')
OUT=pathlib.Path('/home/treefor/sim-workspace/powerflow/outputs/pic-electric-exterior-20261029')
FIXTURE=pathlib.Path('/home/treefor/sim-workspace/powerflow/outputs/pic-profile/cuda-gather-investigation/actual-midpoint-fixture')
QUALIFIED=pathlib.Path('/home/treefor/sim-workspace/powerflow/outputs/pic-profile/em-field-consistency-audit/static-domain-preparation/cpu-qualified')
FROZEN='/home/treefor/projects/runtimes/arch-base-20261028/run-frozen'

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def save(name,value):
    (OUT/name).write_text(json.dumps(value,indent=2)+'\n')

def main():
    assert os.sched_getaffinity(0)=={7},'CPU7 affinity required'
    assert os.getpriority(os.PRIO_PROCESS,0)>=15,'nice15 required'
    records=[]
    def run(argv):
        begin=time.monotonic()
        result=subprocess.run(argv,check=False,text=True,capture_output=True,timeout=25)
        records.append(dict(argv=argv,exit_code=result.returncode,elapsed_s=time.monotonic()-begin,stdout=result.stdout,stderr=result.stderr))
        save('cpu-control-commands.json',dict(CPU=[7],nice=15,one_thread=True,AS_limit_bytes=4294967296,commands=records))
        assert result.returncode==0,records[-1]
        return result.stdout
    cases=[('cpu-preflight',9,0,'uniform-control-sealed'),('cpu-geometry',13,32,'graded-13mm-32-sealed'),('cpu-geometry',17,40,'graded-17mm-40-sealed'),('cpu-geometry',25,48,'graded-25mm-48-sealed'),('cpu-geometry',40,64,'graded-40mm-64-sealed'),('cpu-geometry',80,80,'graded-80mm-80-sealed'),('cpu-geometry',80,120,'graded-80mm-120-sealed'),('cpu-geometry',80,160,'graded-80mm-160-sealed'),('cpu-geometry',160,96,'graded-160mm-96-sealed')]
    cases += [('cpu-geometry',0,n,f'nested-105-{n}-sealed') for n in [40,60,80,96,110,120]]
    if '--finish-existing' not in sys.argv:
        manufactured=run([FROZEN,str(PREFIX/'build/exterior_manufactured')])
        save('manufactured-controls.json',dict(records=[json.loads(line) for line in manufactured.splitlines()],passed=True))
        for mode,width,margin,name in cases:
            run([FROZEN,str(PREFIX/'build/electric_exterior_cpu'),mode,str(FIXTURE),str(OUT/name),str(width),str(margin)])
    else:
        prior=json.loads((OUT/'cpu-control-commands.json').read_text())
        assert len(prior['commands'])==1+len(cases) and all(x['exit_code']==0 for x in prior['commands'])
    pairs={}
    for name in ['full_operator_I.bin','full_operator_J.bin','full_operator_V.bin','constrained_operator_I.bin','constrained_operator_J.bin','constrained_operator_V.bin','charge_rhs','phi_rhs','psi_rhs']:
        old=QUALIFIED/name
        pairs[name]={'new_sha256':sha(OUT/'uniform-control-sealed'/name),'qualified_sha256':sha(old)}
        pairs[name]['equal']=pairs[name]['new_sha256']==pairs[name]['qualified_sha256']
    save('uniform-CSR-RHS-parity.json',dict(passed=all(r['equal'] for r in pairs.values()),comparison='Exact bytes; no tolerances',pairs=pairs))
    assert all(r['equal'] for r in pairs.values()),'Uniform CSR/RHS bytes changed'
    summary=[]
    for _,width,margin,name in cases:
        folder=OUT/name
        result=json.loads((folder/'result.json').read_text());layout=json.loads((folder/'work-and-layout.json').read_text());resources=json.loads((folder/'resources.json').read_text())
        assert result['passed'] and result['nodes']<=1010000 and resources['max_rss_KiB']*1024<4294967296
        if margin:
            transfer=json.loads((folder/'interpolation-control.json').read_text())
            assert transfer['passed'] and transfer['max_E_difference_V_m']==0 and transfer['max_B_difference_T']==0
        summary.append(dict(case=name,width_mm=width,margin_cells=margin,nodes=result['nodes'],cells=result['cells'],nnz=result['nnz'],RSS_KiB=resources['max_rss_KiB'],growth=layout['outer_growth_x'],last_outer_cell_m=layout['last_outer_cell_width_m'],passed=True))
    nested=[OUT/f'nested-105-{n}-sealed'/'axes_x.bin' for n in [40,60,80,96,110,120]]
    import array
    previous=None;previous_n=None
    for path,n in zip(nested,[40,60,80,96,110,120]):
        values=array.array('d');values.frombytes(path.read_bytes())
        if previous is not None:assert values[n-previous_n:n-previous_n+len(previous)]==previous,'Nested interior coordinates changed'
        previous=values;previous_n=n
    save('nested-axis-parity.json',dict(passed=True,shared_interior_axis_coordinates_bitwise_equal=True,fixed_growth=1.05,margin_cells=[40,60,80,96,110,120]))
    save('cpu-controls.json',dict(passed=True,cases=summary,GPU_executed=False,actual_fixture_Poisson_solves=0,PIC_evolution=False))
    print('CPU static exterior controls PASS; actual-fixture solves=0; GPU executions=0')

if __name__=='__main__':main()
