#!/usr/bin/env python3
"""Build and check the isolated scalar C++ adapter under the shared CPU7 lease."""
from pathlib import Path
import csv
import ctypes as C
from decimal import Decimal as D, localcontext
import fcntl
import hashlib
import json
import math
import os
import random
import shutil
import subprocess
import sys
import time

SOURCE = Path(__file__).resolve().parent
OUTPUTS = Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs')
LOCK = OUTPUTS/'pic-profile/testing-core7.lock'
EPS = sys.float_info.epsilon

class Values(C.Structure):
    _fields_ = [(n,C.c_double) for n in (
        'temperature_k','h_j_kg','u_model_j_kg','cp_sensible_j_kg_k',
        'apparent_heat_capacity_j_kg_k','k_w_m_k','resistivity_ohm_m',
        'sigma_s_m','liquid_fraction','energy_j_m3_reference','capacity_j_m3_reference_k')]
    _fields_ += [('phase',C.c_int),('synthesized_high_temperature',C.c_int)]

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def hold():
    p = Path('/proc/278902'); raw = (p/'stat').read_text(); f = raw[raw.rfind(')')+2:].split()
    r = {'pid':278902,'state':f[0],'start_ticks':f[19],'exe':str((p/'exe').resolve()),
         'cgroup_events':dict(line.split() for line in Path('/sys/fs/cgroup/powerflow-pic-gpu-reserved/cgroup.events').read_text().splitlines())}
    assert r['state']=='T' and r['start_ticks']=='5288339'
    assert r['exe']=='/home/treefor/projects/runtimes/mfem-gpu-pic-20261029/qualified/heating-history-pic-gpu'
    assert r['cgroup_events']=={'populated':'1','frozen':'1'}
    return r

def main(report):
    os.sched_setaffinity(0,{7})
    current = os.getpriority(os.PRIO_PROCESS,0)
    if current < 15: os.nice(15-current)
    assert os.getpriority(os.PRIO_PROCESS,0)==15
    for n in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS']:os.environ[n]='1'
    os.environ['CUDA_VISIBLE_DEVICES']=''
    report.mkdir(parents=True,exist_ok=True)
    with LOCK.open('a') as lease:
        fcntl.flock(lease,fcntl.LOCK_EX)
        return audit(report)

def audit(report):
    start=time.time(); before=hold()
    source_before={str(p.relative_to(SOURCE)):sha(p) for p in sorted(SOURCE.rglob('*')) if p.is_file()}
    inputs=json.loads((report/'input-manifest.json').read_text())
    for x in inputs: assert sha(x['shared_path'])==x['sha256']==sha(x['snapshot_path'])
    build=SOURCE.parent/'build';build.mkdir(exist_ok=True)
    compiler=Path(shutil.which('c++')).resolve()
    command=[str(compiler),'-std=c++17','-O2','-fno-fast-math','-Wall','-Wextra','-Werror','-pedantic','-fPIC','-shared',str(SOURCE/'ss304_c_api.cpp'),'-o',str(build/'libss304_thermal.so')]
    proc=subprocess.run(command,capture_output=True,text=True)
    (report/'build.stdout.txt').write_text(proc.stdout);(report/'build.stderr.txt').write_text(proc.stderr)
    assert proc.returncode==0,proc.stderr
    lib=C.CDLL(str(build/'libss304_thermal.so'))
    lib.ss304_error.restype=C.c_char_p
    lib.ss304_evaluate.argtypes=[C.c_int,C.c_double,C.c_double,C.POINTER(Values)]
    lib.ss304_inverse.argtypes=[C.c_int,C.c_double,C.c_double,C.c_int,C.POINTER(C.c_double)]
    lib.ss304_deposit.argtypes=[C.c_int,C.c_double,C.c_double,C.c_double,C.c_double,C.POINTER(Values),C.POINTER(C.c_double)]
    counts={};metrics={};worst={}
    def check(ok,group):
        if not ok:raise AssertionError(group)
        counts[group]=counts.get(group,0)+1
    def close(a,b,group,rtol=64*EPS,atol=0):
        error=abs(a-b);limit=max(atol,rtol*max(abs(a),abs(b)))
        check(error<=limit,group)
        ratio=error/limit if limit else 0
        if ratio>metrics.get(group,0):metrics[group]=ratio;worst[group]={'a':a,'b':b,'absolute_error':error,'limit':limit}
    def ev(t,mode=1,rho=8000):
        v=Values();rc=lib.ss304_evaluate(mode,rho,t,C.byref(v))
        if rc:raise ValueError(lib.ss304_error().decode())
        return v
    def inverse(e,mode=1,rho=8000,vol=0):
        out=C.c_double(-123);rc=lib.ss304_inverse(mode,rho,e,vol,C.byref(out))
        if rc:raise ValueError(lib.ss304_error().decode())
        return out.value
    def rejected(call,group):
        try:call()
        except ValueError:check(True,group)
        else:check(False,group)
    def deposit(u,mass,q,mode=1):
        v=Values();r=C.c_double();rc=lib.ss304_deposit(mode,8000,u,mass,q,C.byref(v),C.byref(r))
        if rc:raise ValueError(lib.ss304_error().decode())
        return v,r.value
    scenario=json.loads((SOURCE/'data/scenario.json').read_text())
    check(sha(SOURCE/'data/ho_chu_304_recovered.csv')==scenario['transport']['table_sha256'],'scenario_source_table_binding')
    with localcontext() as context:
        context.prec=60
        def cp(t):return D(472)+D('.136')*t-D(2820000)/(t*t)
        def hs(t):return (t-D(300))*(D(472)+D('.068')*(t+D(300))-D(2820000)/(t*D(300)))
        cs=cp(D(1672));width=D(55);latent=D(160)*D(1727)
        def reference(t):
            if t<D(1672):return hs(t),cp(t),cp(t),D(0)
            if t<D(1727):
                dt=t-D(1672);f=dt/width;c=cs+(D(800)-cs)*f
                return hs(D(1672))+cs*dt+(D(800)-cs)*dt*dt/(2*width)+latent*f,c,c+latent/width,f
            h=hs(D(1672))+cs*width+(D(800)-cs)*width/2+latent+D(800)*(t-D(1727))
            return h,D(800),D(800),D(1)
        anchors=list(csv.DictReader((SOURCE/'data/ho_chu_304_recovered.csv').open()))
        for row in anchors:
            v=ev(float(row['temperature_k']))
            close(v.resistivity_ohm_m,float(D(row['resistivity_1e_minus8_ohm_m'])*D('1e-8')),'source_resistivity_anchors')
            close(v.k_w_m_k,float(row['thermal_conductivity_w_m_k']),'source_k_anchors')
        for a,b in zip(anchors,anchors[1:]):
            t=(float(a['temperature_k'])+float(b['temperature_k']))/2;v=ev(t)
            target=(D(a['resistivity_1e_minus8_ohm_m'])+D(b['resistivity_1e_minus8_ohm_m']))*D('5e-9')
            close(v.sigma_s_m,float(1/target),'sigma_reciprocal_of_interpolated_resistivity')
        rng=random.Random(20261029)
        points=sorted(set([300.,math.nextafter(300,math.inf),1500.,1672.,1727.,1800.]+[300+1500*rng.random() for _ in range(2000)] +
                          [math.nextafter(x,y) for x in [1500.,1672.,1727.] for y in [-math.inf,math.inf]]))
        previous=None
        for t in points:
            v=ev(t);h,c,ca,f=reference(D.from_float(t))
            for name,target in [('h_j_kg',h),('cp_sensible_j_kg_k',c),('apparent_heat_capacity_j_kg_k',ca),('liquid_fraction',f)]:
                close(getattr(v,name),float(target),'decimal_'+name,atol=1e-28)
            check(v.u_model_j_kg==v.h_j_kg,'explicit_incompressible_energy_identity')
            close(v.energy_j_m3_reference,8000*v.u_model_j_kg,'reference_volume_energy')
            close(v.capacity_j_m3_reference_k,8000*v.apparent_heat_capacity_j_kg_k,'reference_volume_capacity')
            check(all(math.isfinite(getattr(v,n)) for n,_ in Values._fields_),'finite_outputs')
            check(v.k_w_m_k>0 and v.sigma_s_m>0 and v.apparent_heat_capacity_j_kg_k>0,'positive_transport_capacity')
            check(v.synthesized_high_temperature==int(t>1500),'high_temperature_label')
            check(v.phase==(0 if t<1672 else 1 if t<1727 else 2),'phase_label')
            close(inverse(v.h_j_kg),t,'specific_energy_inverse',rtol=4*EPS)
            close(inverse(v.energy_j_m3_reference,vol=1),t,'reference_volume_inverse',rtol=4*EPS)
            if previous and t-previous[0]>1e-8:check(v.h_j_kg>previous[1],'strict_energy_monotonicity')
            previous=(t,v.h_j_kg)
            if t<=1500:close(ev(t,mode=0).h_j_kg,v.h_j_kg,'solid_only_same_law')
        # Finite-volume energy identity: exact Simpson integration of the linear
        # mushy Cp isolates the latent contribution, independently of inverse.
        vm=ev((1672+1727)/2);vs=ev(1672);vl=ev(1727)
        sensible=55*(vs.cp_sensible_j_kg_k+4*vm.cp_sensible_j_kg_k+vl.cp_sensible_j_kg_k)/6
        close((vl.h_j_kg-vs.h_j_kg)-sensible,276320,'fusion_counted_once')
        for t in [350,1000,1600,1680,1700,1720,1750,1790]:
            dt=.01;v=ev(t);dh=(ev(t+dt).h_j_kg-ev(t-dt).h_j_kg)/(2*dt)
            close(dh,v.apparent_heat_capacity_j_kg_k,'dH_dT_equals_apparent_capacity',rtol=1e-8)
        for t in [1672.,1727.]:
            center=ev(t).h_j_kg
            for direction in [-math.inf,math.inf]:
                adjacent=ev(math.nextafter(t,direction)).h_j_kg
                close(adjacent,center,'phase_energy_continuity',rtol=2e-14)
        # Integration over an arbitrary solid interval agrees with a high-order
        # numerical quadrature of source Cp, not merely evaluation of h twice.
        for a,b in [(300.,1500.),(1500.,1672.)]:
            n=1000;step=(b-a)/n
            integral=step/3*math.fsum((1 if i in [0,n] else 4 if i%2 else 2)*ev(a+i*step).cp_sensible_j_kg_k for i in range(n+1))
            close(integral,ev(b).h_j_kg-ev(a).h_j_kg,'integrated_solid_cp_energy',rtol=1e-11)
        for lower,upper in [(350,1200),(1200,1700),(1700,1750),(1750,400),(1700,1680)]:
            old,new=ev(lower),ev(upper);mass=2.3e-8;q=mass*(new.u_model_j_kg-old.u_model_j_kg)
            v,residual=deposit(old.u_model_j_kg,mass,q)
            close(v.temperature_k,float(upper),'signed_heat_transaction',rtol=4*EPS)
            close(mass*(v.u_model_j_kg-old.u_model_j_kg)-q,residual,'reported_energy_residual',atol=1e-18)
            check(abs(residual)<=2e-17,'heat_energy_balance_roundoff')
        for t in [math.nextafter(300,-math.inf),math.nextafter(1800,math.inf),-1,0,math.inf,-math.inf,math.nan]:rejected(lambda t=t:ev(t),'temperature_domain_rejection')
        rejected(lambda:ev(math.nextafter(1500,math.inf),mode=0),'solid_mode_extrapolation_rejection')
        for mode in [0,1]:
            top=ev(1500 if mode==0 else 1800,mode=mode)
            for vol,maximum in [(0,top.h_j_kg),(1,top.energy_j_m3_reference)]:
                for e in [math.nextafter(0,-math.inf),math.nextafter(maximum,math.inf),math.inf,-math.inf,math.nan]:rejected(lambda e=e:inverse(e,mode=mode,vol=vol),'energy_domain_rejection')
        for rho in [0,-1,math.inf,math.nan,sys.float_info.max]:rejected(lambda rho=rho:ev(500,rho=rho),'density_rejection')
        for mode in [-1,2]:rejected(lambda mode=mode:ev(500,mode=mode),'qualification_rejection')
        for mass,q in [(0,1),(-1,1),(math.inf,1),(1,math.nan),(1,math.inf)]:rejected(lambda:deposit(0,mass,q),'invalid_deposition_rejection')
        rejected(lambda:deposit(ev(1800).u_model_j_kg,1,1e-300),'upper_endpoint_positive_heat_rejection')
        rejected(lambda:deposit(0,1,-1e-300),'lower_endpoint_cooling_rejection')
        rejected(lambda:deposit(ev(1700).u_model_j_kg,1,1e99),'overshoot_no_clamp')
        tiny,residual=deposit(0,1,1e-300)
        check(tiny.temperature_k==300 and residual==-1e-300,'subresolution_heat_explicit_residual')
        for t in [300.,1500.,1672.,1700.,1727.,1800.]:
            old=ev(t);unchanged,residual=deposit(old.u_model_j_kg,1e-8,0)
            check(unchanged.temperature_k==t and residual==0,'zero_deposit_identity')
        # C ABI must leave caller outputs untouched on failed trial.
        v=Values();C.memset(C.byref(v),0x2a,C.sizeof(v));raw=bytes(v)
        check(lib.ss304_evaluate(1,8000,1900,C.byref(v))!=0 and bytes(v)==raw,'failed_trial_output_atomicity')
        r=C.c_double(123);check(lib.ss304_inverse(1,8000,-1,0,C.byref(r))!=0 and r.value==123,'failed_inverse_output_atomicity')
        check(lib.ss304_deposit(1,8000,0,1,-1,C.byref(v),C.byref(r))!=0 and bytes(v)==raw and r.value==123,'failed_deposit_output_atomicity')
        check(lib.ss304_evaluate(1,8000,500,None)!=0,'null_output_rejection')
        check(lib.ss304_inverse(1,8000,0,3,C.byref(r))!=0,'energy_measure_rejection')
        # Source latent/cp choices are approximations, not extra defaults from EOS.
        check(ev(1727).phase==2 and ev(1124.3).phase==0,'eos_parameter_is_not_melting_temperature')
        snapshots=[]
        for t in [300,773.15,1000,1500,1600,1672,1700,1727,1800]:
            v=ev(t);snapshots.append({n:getattr(v,n) for n,_ in Values._fields_})
    after=hold()
    for x in inputs:assert sha(x['shared_path'])==x['sha256']==sha(x['snapshot_path'])
    assert source_before=={str(p.relative_to(SOURCE)):sha(p) for p in sorted(SOURCE.rglob('*')) if p.is_file()}
    linked=subprocess.check_output(['ldd',str(build/'libss304_thermal.so')],text=True)
    assert not any(name in linked.lower() for name in ['cuda','cublas','cudss','mfem'])
    result={'status':'passed_conditional_nominal_thermal_property_checks_not_coupled_physics_validation',
            'checks':sum(counts.values()),'check_groups':counts,'worst_fraction_of_numeric_limit':metrics,'worst_numeric_pairs':worst,
            'source_and_model_points':snapshots,'compiler':str(compiler),'compiler_sha256':sha(compiler),
            'compiler_version':subprocess.check_output([str(compiler),'--version'],text=True).splitlines()[0],
            'build_command':command,'library_sha256':sha(build/'libss304_thermal.so'),
            'linked_dependencies':linked,'source_files_sha256':source_before,'source_unchanged_during_checks':True,
            'cpu_affinity':sorted(os.sched_getaffinity(0)),'nice':os.getpriority(os.PRIO_PROCESS,0),'threads':1,
            'lease':str(LOCK),'production_before':before,'production_after':after,'shared_inputs_unchanged':True,
            'gpu_or_production_signals':False,'wall_s':time.time()-start,
            'numeric_check_tolerances_are_not_source_uncertainty_or_changed_PIC_guards':True,
            'limits':['fixed_reference_volume_Cp_as_Cv_caloric_approximation','pooled_304_304L_nominal_transport','source_synthesis_above1500K','assumed_mushy_fraction_and_sensible_Cp_bridge','no_vapor_or_state_above1800K','no_mechanics_or_773K_guard_relaxation','not_foil_calibration_or_coupled_breakdown_validation']}
    (report/'verification.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ['status','checks','wall_s','cpu_affinity','nice','library_sha256']}))

if __name__=='__main__':main(Path(sys.argv[1]).resolve())
