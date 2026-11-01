#!/usr/bin/python -I
import hashlib,json,os,pathlib,subprocess,sys,time
if not __debug__: raise RuntimeError('optimized Python guard execution refused')
R=pathlib.Path(__file__).resolve().parents[1]
OUT=pathlib.Path('/home/treefor/projects/.workspace/outputs/gas-target-owner-v1.4-binding-20261101')
FORMS=R.parent/'pic-reaction-transactions-20261029/native-full-time-candidate008-auditfix/forms/include/pf_v7'
GPU=R.parent/'hall-gpu-resident-v1-20261030'
NV=R.parent/'pic-cuda-gather-investigation-20261029/tools/nvrtc/nvidia/cuda_nvrtc'
CUDA=R.parent/'mfem-cuda-benchmark-20261029/cuda-libs/nvidia/cuda_runtime/include'
MODE=sys.argv[1];assert MODE in ('normal','sanitized','device');OUT.mkdir(parents=True,exist_ok=True);D=OUT/sys.argv[2];D.mkdir(exist_ok=False)
assert os.getenv('SIMQ_JOB_ID','').isdigit(),'SimQueue required'
unit=os.getenv('SIMQ_UNIT','');affinity=sorted(os.sched_getaffinity(0))
assert unit and unit in pathlib.Path('/proc/self/cgroup').read_text() and os.getenv('SIMQ_LAUNCH_TOKEN')
assert len(affinity)==1 and os.getenv('SIMQ_CPU_LIST')==str(affinity[0])
assert os.getenv('SIMQ_PARTITION')=='powerflow'
assert affinity[0] in (6,7) or (affinity[0] in (4,5) and os.getenv('SIMQ_BORROW')=='1')
assert all(os.getenv(x)=='1' for x in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'))
mem=int(os.getenv('SIMQ_MEM','0'));assert mem>0
meminfo={line.split(':')[0]:int(line.split()[1])*1024 for line in pathlib.Path('/proc/meminfo').read_text().splitlines() if line.split()[1].isdigit()}
assert meminfo['MemAvailable']>=mem+2*1024**3,'physical memory reserve unavailable'
if MODE=='device':
    assert 'GPU-a7b9b458-167f-cd57-fec9-3dea0995b6e9' in os.getenv('SIMQ_GPUS','') and int(os.getenv('SIMQ_GPU_MEM','0'))>0,'scheduled shared GPU allocation required'

FROZEN=R.parent/'arch-base-20261028'
root_stat=pathlib.Path('/').stat();frozen_root_stat=(FROZEN/'root').stat()
root_identity=(root_stat.st_dev,root_stat.st_ino)
assert root_identity==(frozen_root_stat.st_dev,frozen_root_stat.st_ino),'outer run-frozen filesystem required'
assert pathlib.Path('/proc/sys/kernel/osrelease').read_text().strip()==(FROZEN/'kernel.release').read_text().strip()
for executable in ('usr/bin/c++','usr/bin/python'):
    assert hashlib.sha256(pathlib.Path('/'+executable).read_bytes()).digest()==hashlib.sha256((FROZEN/'root'/executable).read_bytes()).digest(),'frozen executable required'

def hashes():
    paths=[p for d in ('include','src','device','tests','tools','baseline','reference-v1.1') for p in (R/d).rglob('*') if p.is_file()]
    paths += [FORMS/'forms_core.h',GPU/'common/hall_gpu_common.hpp',GPU/'common/hall_gpu_common.cpp']
    return {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}
before=hashes();receipt={'frozen_root_identity':root_identity,'frozen_root_verified':True,'mode':MODE,'job':os.environ['SIMQ_JOB_ID'],'affinity':sorted(os.sched_getaffinity(0)),'input_hashes':before,'passed':False,'runs':[],'allocation':{k:v for k,v in os.environ.items() if k.startswith('SIMQ_') and k!='SIMQ_LAUNCH_TOKEN'},'physical_memory_before':meminfo,'frozen_launcher_sha256':hashlib.sha256((FROZEN/'run-frozen').read_bytes()).hexdigest()}
def run(cmd,name,env=None):
    t=time.monotonic()
    with (D/(name+'.log')).open('w') as f:r=subprocess.run(list(map(str,cmd)),stdout=f,stderr=subprocess.STDOUT,env=env)
    receipt['runs'].append({'command':list(map(str,cmd)),'rc':r.returncode,'wall_s':time.monotonic()-t,'log':str(D/(name+'.log'))})
    print(name,r.returncode,flush=True)
    if r.returncode:raise RuntimeError(name+' failed: '+(D/(name+'.log')).read_text()[-4000:])
try:
    cmd=['/usr/bin/c++','-std=c++20','-Wall','-Wextra','-Werror','-ffp-contract=off','-fno-fast-math','-I',R/'include','-I',FORMS]
    if MODE=='sanitized':cmd+=['-O1','-g','-fno-omit-frame-pointer','-fsanitize=address,undefined','-fno-sanitize-recover=all']
    else:cmd+=['-O2']
    cmd += [R/'tests/test_guards.cpp',R/'tests/allocation_probe.cpp',R/'src/gas_owner.cpp']
    if MODE=='device':
        cmd+=[R/'tests/device_begin_fault.cpp','-DPF_GAS_BEGIN_FAULT_WRAPPER','-Wl,--wrap=_ZN9pf_gas_v110DeviceLane10begin_stepEm','-Wno-c++26-extensions','-DPF_GAS_WITH_DEVICE','-DPF_GAS_FORMS_CORE="'+str(FORMS/'forms_core.h')+'"','-I',GPU,'-isystem',CUDA,'-isystem',NV/'include',R/'src/gas_device.cpp',GPU/'common/hall_gpu_common.cpp','-ldl','-lcrypto']
    cmd+=['-o',D/'test_guards']
    run(cmd,'build');receipt['binary_sha256']=hashlib.sha256((D/'test_guards').read_bytes()).hexdigest()
    run([D/'test_guards'],'host')
    if MODE=='normal':
        for name,base in (('baseline',R/'baseline'),('fixed',R)):
            run(['/usr/bin/c++','-std=c++20','-O2','-ffp-contract=off','-fno-fast-math','-I',base/'include','-I',FORMS,R/'tests/test_preimages.cpp',base/'src/gas_owner.cpp','-o',D/('preimage-'+name)],'build-preimage-'+name)
            run([D/('preimage-'+name),name],'reproduce-'+name)
        for name,base in (('baseline',R/'reference-v1.1'),('fixed',R)):
            run(['/usr/bin/c++','-std=c++20','-O2','-ffp-contract=off','-fno-fast-math','-I',base/'include','-I',FORMS,R/'tests/test_motion_preimages.cpp',base/'src/gas_owner.cpp','-o',D/('motion-preimage-'+name)],'build-motion-preimage-'+name)
            run([D/('motion-preimage-'+name),name],'motion-reproduce-'+name)

    if MODE=='device':
        run([D/'test_guards','device'],'device');run([D/'test_guards','verified'],'verified')
        states=[next(l for l in (D/(n+'.log')).read_text().splitlines() if l.startswith('state ')) for n in ('host','device','verified')]
        assert len(set(states))==1,'state parity differs'
        motion_states=[next(l for l in (D/(n+'.log')).read_text().splitlines() if l.startswith('motion_state ')) for n in ('host','device','verified')]
        assert len(set(motion_states))==1,'transport state parity differs'
        publication_states=[next(l for l in (D/(n+'.log')).read_text().splitlines() if l.startswith('publication_state ')) for n in ('host','device','verified')]
        assert len(set(publication_states))==1,'publication state parity differs'
        binding_states=[next(l for l in (D/(n+'.log')).read_text().splitlines() if l.startswith('binding_state ')) for n in ('host','device','verified')]
        assert len(set(binding_states))==1,'binding-path state parity differs'
        receipt['state_parity']=True
    receipt['binary_unchanged']=hashlib.sha256((D/'test_guards').read_bytes()).hexdigest()==receipt['binary_sha256']
    assert receipt['binary_unchanged']
    receipt['inputs_unchanged']=hashes()==before
    assert receipt['inputs_unchanged']
    receipt['passed']=True
finally:
    (D/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
