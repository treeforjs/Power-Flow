#!/usr/bin/env python3
import hashlib,json,os,pathlib,shlex,subprocess,sys,time

if not __debug__:raise RuntimeError('optimized guard execution refused')
R=pathlib.Path(__file__).resolve().parents[1]
OUT=R.parent.parent/'.workspace/outputs'/R.name
BASE=R.parent/'pic-arrival-owner-v1.1-20261030'
CANONICAL=R.parent/'pic-arrival-owner-v1-20261030/include'
NATIVE=R.parent/'pic-reaction-transactions-20261029/native-full-time-candidate008-auditfix/source/include'
FROZEN=R.parent/'arch-base-20261028'
mode,label=sys.argv[1:];assert mode in ('normal','sanitized')
assert os.getenv('SIMQ_JOB_ID','').isdigit() and os.getenv('SIMQ_LAUNCH_TOKEN')
unit=os.getenv('SIMQ_UNIT','');assert unit and unit in pathlib.Path('/proc/self/cgroup').read_text()
cores=sorted(os.sched_getaffinity(0));assert len(cores)==1 and os.getenv('SIMQ_CPU_LIST')==str(cores[0])
assert os.getenv('SIMQ_PARTITION')=='powerflow'
assert cores[0] in (6,7) or (cores[0] in (4,5) and os.getenv('SIMQ_BORROW')=='1')
assert all(os.getenv(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'))
assert os.getenv('CUDA_VISIBLE_DEVICES')==''
mem=int(os.getenv('SIMQ_MEM','0'));assert mem>0
available=next(int(line.split()[1])*1024 for line in pathlib.Path('/proc/meminfo').read_text().splitlines() if line.startswith('MemAvailable:'))
assert available>=mem+2*1024**3
st=pathlib.Path('/').stat();f=(FROZEN/'root').stat();assert (st.st_dev,st.st_ino)==(f.st_dev,f.st_ino),'actual frozen root required'
assert pathlib.Path('/proc/sys/kernel/osrelease').read_text().strip()==(FROZEN/'kernel.release').read_text().strip()
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def physical(p):
    p=pathlib.Path(p)
    return FROZEN/'root'/p.relative_to('/') if str(p).startswith(('/usr/','/lib/')) else p
def inputs():
    files=[p for d in ('include','tests','tools') for p in (R/d).rglob('*') if p.is_file()]
    files += list((BASE/'include').rglob('*.hpp'))+list(CANONICAL.rglob('*.hpp'))+list(CANONICAL.rglob('*.h'))
    files += [NATIVE/'raikou/exact_particle_count.h',FROZEN/'run-frozen',pathlib.Path('/usr/bin/c++'),pathlib.Path('/usr/lib/libcrypto.so.3')]
    return {str(physical(p)):sha(physical(p)) for p in sorted(files)}
before=inputs();D=OUT/label;D.mkdir(parents=True,exist_ok=False)
receipt={'mode':mode,'job':os.environ['SIMQ_JOB_ID'],'affinity':cores,'frozen_root':[st.st_dev,st.st_ino],'input_hashes':before,'runs':[],'passed':False,'allocation':{k:v for k,v in os.environ.items() if k.startswith('SIMQ_') and k!='SIMQ_LAUNCH_TOKEN'}}
def run(cmd,name):
    t=time.monotonic()
    with (D/(name+'.log')).open('w') as log:r=subprocess.run(list(map(str,cmd)),stdout=log,stderr=subprocess.STDOUT)
    receipt['runs'].append({'command':list(map(str,cmd)),'returncode':r.returncode,'wall_s':time.monotonic()-t,'log':str(D/(name+'.log'))})
    print(name,r.returncode,flush=True)
    if r.returncode:raise RuntimeError((D/(name+'.log')).read_text()[-5000:])
try:
    deps={};binaries={}
    for name,root in [('baseline',BASE),('prepared',R)]:
        cmd=['/usr/bin/c++','-std=c++20','-Wall','-Wextra','-Werror','-ffp-contract=off','-fno-fast-math','-I',root/'include','-I',CANONICAL,'-I',NATIVE,'-DPF_ARRIVAL_EXACT_COUNT_HEADER="raikou/exact_particle_count.h"']
        cmd+=['-O1','-g','-fno-omit-frame-pointer','-fsanitize=address,undefined','-fno-sanitize-recover=all'] if mode=='sanitized' else ['-O2']
        if name=='prepared':cmd+=['-DPF_TEST_PREPARED=1']
        binary=D/name
        objects=[]
        for tu in ['test_publication','allocation_probe']:
            obj=D/(name+'-'+tu+'.o');dep=D/(name+'-'+tu+'.d')
            run(cmd+['-MD','-MF',dep,'-c',R/'tests'/(tu+'.cpp'),'-o',obj],'build-'+name+'-'+tu)
            objects.append(obj)
            for path in shlex.split(dep.read_text().split(':',1)[1].replace('\\\n',' ')):
                p=pathlib.Path(path);p=p if p.is_absolute() else R/p;p=physical(p);deps[str(p)]=sha(p)
        run(cmd+objects+['-lcrypto','-o',binary],'link-'+name);binaries[str(binary)]=sha(binary)
        run([binary],name)
    histories=[next(line for line in (D/(name+'.log')).read_text().splitlines() if line.startswith('history ')) for name in ('baseline','prepared')]
    assert len(set(histories))==1,'actual baseline/checkpoint/ledger parity differs'
    assert inputs()==before and all(sha(p)==h for p,h in deps.items()) and all(sha(p)==h for p,h in binaries.items())
    receipt.update({'passed':True,'inputs_unchanged':True,'baseline_parity':histories[0],'binary_hashes':binaries,'compiler_dependency_hashes':deps})
finally:(D/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
