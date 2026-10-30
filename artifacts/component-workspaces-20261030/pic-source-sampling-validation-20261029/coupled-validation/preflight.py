#!/usr/bin/env python3
"""Build and exercise CPU-only validation controls. No GPU/coupled launch path."""
import argparse
import csv
import fcntl
import hashlib
import json
import os
import shutil
import struct
import subprocess
import sys
import time
from pathlib import Path
import metrics
import budget

ROOT=Path(__file__).resolve().parent
OUT=Path('<workspace>/outputs/pic-profile/source-sampling-validation/coupled-preflight')
PROFILE=OUT.parent.parent
ORIGINAL=PROFILE/'particle-audit/copied-checkpoint'
FROZEN=Path('/home/treefor/projects/runtimes/mfem-gpu-pic-20261029')
PROPOSAL=OUT.parent/'coupled-validation-proposal.json'
THERMAL=Path('<workspace>/outputs/heating-history-full.bin')
VOLTAGE=Path('<workspace>/outputs/gap-pulse-midplane110kV.csv')
FIELD=Path('/home/treefor/projects/Power-Flow/runs/restore_field_mpi8_domain40mm')
CHECKPOINT='checkpoint_400040.bin'
CHECK_SHA='957da7e1419c8ca09d844265312bef1e850a74fc93d144ca557422aad39c752a'
META_SHA='93f1c94b61f50efe9452ac63e4a38edf355f19d2ffdeeb542a14f8f08c0fe2ba'

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def load(path):
    def pairs(items):
        out={}
        for k,v in items:
            if k in out:raise ValueError('Duplicate JSON key')
            out[k]=v
        return out
    return json.loads(Path(path).read_text(),object_pairs_hook=pairs,parse_constant=lambda x:(_ for _ in ()).throw(ValueError('Nonfinite JSON')))

def owned(path):
    path=Path(path).resolve();path.relative_to(OUT.resolve());return path

def write(path,value):
    owned(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')

def command(args,log,expected=0,env=None,timeout=90):
    p=subprocess.run([str(x) for x in args],capture_output=True,text=True,env=env,timeout=timeout)
    owned(log).write_text(p.stdout+p.stderr)
    if p.returncode!=expected:raise RuntimeError(f'Command returned {p.returncode}, expected {expected}: {args[0]}; see {log}')
    return p.stdout

def build(where):
    command([sys.executable,ROOT/'prepare.py'],where/'prepare.log')
    base=Path('/home/treefor/projects/runtimes/mfem-4.10-pic-20261029')
    includes=['-I'+str(ROOT/'frozen/include'),'-I'+str(ROOT/'frozen/src'),'-I'+str(ROOT/'frozen/examples'),'-I'+str(ROOT)]
    for p in [base/'install-mpi/include',base/'tpl/hypre/include',base/'tpl/metis/include']:includes+=['-isystem',str(p)]
    flags=['/usr/sbin/c++','-std=c++20','-O3','-DNDEBUG','-Wall','-Wextra','-Werror=return-type']
    command(flags+includes+['-c',ROOT/'mfem_gap_pic_candidate.cpp','-o',ROOT/'pic.o'],where/'build-native.log')
    libraries=[base/'install-mpi/lib/libmfem.a',base/'tpl/hypre/lib/libHYPRE.a',base/'tpl/metis/lib/libmetis.a',base/'tpl/GKlib/lib/libGKlib.a','/usr/lib/libmpi.so','-Wl,-rpath,/usr/lib']
    command(flags+includes+[ROOT/'trial.cpp',ROOT/'pic.o','-o',ROOT/'coupled_trial']+libraries,where/'build-driver.log')
    command(flags+['-I'+str(ROOT),ROOT/'source_preflight.cpp','-o',ROOT/'source_preflight'],where/'build-source.log')

def proof():
    manifest=load(ROOT/'candidate-manifest.json')
    if manifest['gpu_execution_enabled'] or manifest['coupled_execution_enabled']:raise ValueError('Execution fence changed')
    for name,want in manifest['frozen_source_sha256'].items():
        if digest(name)!=want:raise ValueError('Frozen source drift')
    for name,want in manifest['generated_files_sha256'].items():
        if digest(ROOT/name)!=want:raise ValueError('Owned generated source drift')
    original=(FROZEN/'source/src/mfem_gap_pic.cpp').read_text()
    transformed=(ROOT/'mfem_gap_pic_candidate.cpp').read_text()
    restored=transformed.replace('#include "gpu_disabled.hpp"','#include "cudss_poisson.hpp"').replace('\n            if(absorption_observer_)absorption_observer_(p,hit);','')
    if restored!=original:raise ValueError('Native physics/guard code changed beyond declared hooks')
    candidate=(ROOT/'frozen/include/raikou/mfem_gap_pic.h').read_text()
    candidate=candidate.replace('\n    void set_absorption_observer(std::function<void(const GapParticle&,int)> f) { absorption_observer_=std::move(f); }\n','')
    candidate=candidate.replace('\n    std::function<void(const GapParticle&,int)> absorption_observer_;','')
    if candidate!=(FROZEN/'source/include/raikou/mfem_gap_pic.h').read_text():raise ValueError('Native header changed beyond passive observer')
    if digest(ROOT/'frozen/src/mfem_pic_distributed.hpp')!=digest(FROZEN/'source/src/mfem_pic_distributed.hpp'):raise ValueError('Poisson guard changed')
    return manifest

def field_support(start,end):
    index=FIELD/'em_field_snapshots.csv'
    with index.open() as f:rows=list(csv.reader(f))[1:]
    times=[float(r[0]) for r in rows]
    if not times or any(b<=a for a,b in zip(times,times[1:])) or start<times[0] or end>times[-1]:raise ValueError('Field clock does not cover trial')
    selected={index,FIELD/'em_mesh.mesh',FIELD/'em_axes.csv'}
    import bisect
    first=max(0,bisect.bisect_right(times,start)-1);last=min(len(times)-1,bisect.bisect_left(times,end))
    for row in rows[first:last+1]:
        for name in row[1:3]:
            p=(FIELD/name).resolve();p.relative_to(FIELD.resolve())
            if not p.is_file():raise ValueError('Missing field snapshot')
            selected.add(p)
    return {str(p):digest(p) for p in sorted(selected)}

def descriptors(proposal,inputs,identity):
    jobs=[]
    for stage in proposal['stages']:
        seeds=[0] if stage['name']=='resource_pilot' else stage['screen_seed_indices']+stage['additional_seed_indices']
        for seed in seeds:
            for count in [64,32,16]:
                jobs.append(dict(stage=stage['name'],seed=seed,count=count,reference_count=64,
                    mode='adaptive' if stage['name']=='native_adaptive_confirmation' else 'reference' if count==64 else 'replay',
                    input_directory=str(inputs[count]),input_identity=identity,
                    checkpoint_sha256=CHECK_SHA,restart_sha256=META_SHA,
                    max_accepted_steps=stage['accepted_steps_per_run'],wall_cap_s=stage['wall_seconds_per_run'],
                    duration_ps=stage['duration_ps'],output_cap_bytes=256*1024*1024,
                    host_RSS_cap_bytes=4*1024**3,
                    gpu_execution_enabled=False,coupled_execution_enabled=False,
                    conservative_ceiling_s=None,ceiling_qualified=False,
                    reference_trace_sha256=None,reference_trace_qualified=False,
                    execution_status='disabled_not_launched'))
    if len(jobs)!=39:raise ValueError('Campaign budget changed')
    return jobs

def run(args):
    began=time.monotonic()
    if os.sched_getaffinity(0)!={6} or os.getpriority(os.PRIO_PROCESS,0)!=15:raise RuntimeError('CPU6/nice15 required')
    if not Path('/proc/self/fd/9').exists():raise RuntimeError('Use build-and-preflight.sh with inherited testing lease')
    lock=PROFILE/'testing-core.lock'
    if os.fstat(9).st_ino!=lock.stat().st_ino or os.fstat(9).st_dev!=lock.stat().st_dev:raise RuntimeError('Wrong lease file descriptor')
    fcntl.flock(9,fcntl.LOCK_EX|fcntl.LOCK_NB)
    out=owned(args.output);out.mkdir(exist_ok=True)
    i=1
    while (out/f'attempt-{i:03d}').exists():i+=1
    where=out/f'attempt-{i:03d}';where.mkdir()
    write(where/'progress.json',dict(phase='preparing',coupled_execution_enabled=False))
    protected=[ORIGINAL/CHECKPOINT,ORIGINAL/'restart.meta',THERMAL,VOLTAGE,FROZEN/'qualified/heating-history-pic-gpu',PROPOSAL]
    before={str(p):digest(p) for p in protected}
    if before[str(ORIGINAL/CHECKPOINT)]!=CHECK_SHA or before[str(ORIGINAL/'restart.meta')]!=META_SHA:raise ValueError('Copied checkpoint identity changed')
    seals={r['path']:r['sha256'] for r in load(FROZEN/'sealed-files.json')}
    for p in protected:
        if str(p) in seals and before[str(p)]!=seals[str(p)]:raise ValueError('Protected seal mismatch')
    proposal=load(PROPOSAL)
    if digest(proposal['proposal_path'])!=proposal['proposal_sha256']:raise ValueError('Proposal text/spec mismatch')
    if proposal['counts']!=[64,32,16] or proposal['resources']['total_max_runs']!=39:raise ValueError('Proposal scope changed')
    if args.build:build(where)
    manifest=proof()
    linkages={}
    for name in ['coupled_trial','source_preflight']:
        text=command(['ldd',ROOT/name],where/f'linkage-{name}.txt')
        if any(x in text.lower() for x in ['cuda','cudss','nvidia','nvml','nvrtc','cublas']):raise ValueError('GPU library dependency detected')
        linkages[name]=text
    advertised=json.loads(command([ROOT/'coupled_trial','--describe'],where/'describe.json'))
    if advertised['gpu_execution_enabled'] or advertised['coupled_execution_enabled']:raise ValueError('Execution enabled')
    inputs={}
    for count in [64,32,16]:
        dest=where/'inputs'/str(count);dest.mkdir(parents=True)
        for name in [CHECKPOINT,'restart.meta']:
            shutil.copyfile(ORIGINAL/name,dest/name)
            if digest(dest/name)!=before[str(ORIGINAL/name)]:raise ValueError('Restart clone differs')
        inputs[count]=dest
        d=where/f'native-{count}';d.mkdir()
        command([ROOT/'coupled_trial','--preflight',dest,d,THERMAL],where/f'native-{count}.log')
        parsed=load(d/'native-preflight.json')
        if not parsed['passed'] or parsed['field_solves'] or parsed['trajectory_steps'] or parsed['MPI_initialized']:raise ValueError('Preflight executed coupled work')
        if (d/'native-particles.bin').read_bytes()!=(ORIGINAL/CHECKPOINT).read_bytes()[392:]:raise ValueError('Initial population was altered')
    source=where/'source';source.mkdir()
    command([ROOT/'source_preflight',inputs[64],OUT.parent/'events.bin',source],where/'source.log')
    emitted=load(source/'source-preflight.json')
    if not emitted['passed'] or emitted['coupled_trace_qualified']:raise ValueError('Source-only trace promoted to coupled evidence')
    fixtures=where/'negative-fixtures';fixtures.mkdir();negative=[]
    # Malformed restart files are confined to owned fixture directories.
    for kind in ['truncated','trailing','duplicate-id','clock-mismatch','unsafe-name']:
        p=fixtures/kind;p.mkdir();shutil.copyfile(inputs[64]/'restart.meta',p/'restart.meta');data=bytearray((inputs[64]/CHECKPOINT).read_bytes())
        if kind=='truncated':data=data[:-1]
        elif kind=='trailing':data+=b'x'
        elif kind=='duplicate-id':data[472:480]=data[392:400]
        elif kind=='clock-mismatch':
            text=(p/'restart.meta').read_text();parts=text.split();parts[-2]=str(float(parts[-2])+1e-9);(p/'restart.meta').write_text(' '.join(parts)+'\n')
        elif kind=='unsafe-name':(p/'restart.meta').write_text('../forbidden.bin 400040 27\n')
        (p/CHECKPOINT).write_bytes(data)
        dest=p/'output';dest.mkdir()
        command([ROOT/'coupled_trial','--preflight',p,dest,THERMAL],where/f'negative-{kind}.log',expected=78)
        if any(dest.iterdir()):raise ValueError('Rejected restart wrote result data')
        negative.append(kind)
    inspected=json.loads(command([ROOT/'coupled_trial','--inspect-trace',source/'manufactured-trace.bin'],where/'trace-inspection.json'))
    if inspected['coupled_verified']:raise ValueError('Manufactured trace qualified itself')
    bad=fixtures/'truncated-trace.bin';bad.write_bytes((source/'manufactured-trace.bin').read_bytes()[:-1])
    command([ROOT/'coupled_trial','--inspect-trace',bad],where/'negative-trace.log',expected=78);negative.append('truncated-trace')
    env=os.environ.copy();env.update(RAIKOU_PIC_FIELD_BACKEND='cudss',CUDA_VISIBLE_DEVICES='0',VALIDATION_ENABLE_GPU='1')
    for tag,extra in [('run',['--run']),('flag',['--run','--enable-gpu']),('config',['--run','--config','ignored.json'])]:
        command([ROOT/'coupled_trial',*extra],where/f'denied-{tag}.log',expected=78,env=env)
    denied_config=fixtures/'gpu-enabled.json';write(denied_config,dict(gpu_execution_enabled=True))
    try:check_execution_request(load(denied_config))
    except RuntimeError:negative.append('enabled-config-rejected')
    else:raise ValueError('GPU config enabled execution')
    outside=ROOT/'path-rejection-fixture';outside.mkdir(exist_ok=False)
    try:
        command([ROOT/'coupled_trial','--preflight',inputs[64],outside,THERMAL],where/'denied-output-path.log',expected=78)
        if any(outside.iterdir()):raise ValueError('Out-of-scope output written')
    finally:outside.rmdir()
    negative.append('output-path-escape-rejected')
    # Source-identity drift must fail before any numerical operation.
    current=ROOT/'native_driver_parts.hpp';original_bytes=current.read_bytes()
    try:
        current.write_bytes(original_bytes.replace(b'.05*pic_me',b'.10*pic_me',1))
        try:proof()
        except ValueError:negative.append('owned-guard-tamper-rejected')
        else:raise ValueError('Changed guard accepted')
    finally:current.write_bytes(original_bytes)
    metric_checks=metrics.self_test(where);write(where/'metric-checks.json',metric_checks)
    budget_checks=budget.self_test();write(where/'campaign-budget-checks.json',budget_checks)
    t0=proposal['common_start_time_ns']*1e-9
    field_hashes=field_support(t0,t0+100e-12)
    with VOLTAGE.open() as f:voltage=list(csv.reader(f))[1:]
    if float(voltage[0][0])>t0 or float(voltage[-1][1])<t0+100e-12:raise ValueError('Voltage clock coverage missing')
    identity=hashlib.sha256(json.dumps(dict(inputs=before,fields=field_hashes,source=manifest['frozen_source_sha256']),sort_keys=True).encode()).hexdigest()
    jobs=descriptors(proposal,inputs,identity);write(where/'campaign-plan.json',dict(execution_enabled=False,jobs=jobs))
    after={str(p):digest(p) for p in protected}
    if before!=after:raise ValueError('Protected input/executable changed')
    report=dict(status='CPU_preflight_passed_GPU_and_coupled_execution_disabled',counts=[64,32,16],
        native_restart_copies_identical=True,initial_macros=28098,next_particle_id=69633,
        source_preflight=emitted,negative_checks=negative,metric_checks=metric_checks,budget_checks=budget_checks,
        MPI_initialized=False,field_solves=0,trajectory_steps=0,GPU_execution_enabled=False,
        coupled_execution_enabled=False,convergence_established=False,protected_before_after_sha256=before,
        selected_field_input_sha256=field_hashes,CPU=6,nice=15,threads=1,local_jobs=1,
        testing_lock_held=True,wall_s=time.monotonic()-began,
        limitations=['No real field/trajectory evolution or coupled trace was qualified','Manufactured guard B bounds do not qualify an actual step ceiling','GPU backend is an unconditional throwing stub','GPU/host memory caps require a later authorized pilot','Confidence kernels were tested on manufactured scalars; no physical campaign result exists'])
    write(where/'report.json',report)
    write(where/'progress.json',dict(phase='complete',GPU_execution_enabled=False,coupled_execution_enabled=False))
    write(out/'latest.json',dict(result_directory=str(where),report=str(where/'report.json')))
    write(where/'seal.json',dict(files_sha256={str(p.relative_to(where)):digest(p) for p in sorted(where.rglob('*')) if p.is_file()},
                               runtime_sha256={p.name:digest(p) for p in ROOT.iterdir() if p.is_file()}))
    print(json.dumps(dict(status=report['status'],result_directory=str(where),wall_s=report['wall_s']),indent=2))

def check_execution_request(config):
    raise RuntimeError('Coupled/GPU execution disabled; no runtime configuration can enable it')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--build',action='store_true');parser.add_argument('--output',type=Path,default=OUT)
    args=parser.parse_args()
    try:run(args)
    except Exception as e:print(str(e),file=sys.stderr);raise
