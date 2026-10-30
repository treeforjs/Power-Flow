#!/usr/bin/env python3
"""Qualify conditional event data/kinematics on one leased CPU; no evolution."""
import fcntl
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]
OUT=Path('/home/treefor/sim-workspace/powerflow/outputs/h2o-collision-data-20261029')
FROZEN=Path('/home/treefor/projects/runtimes/arch-base-20261028/run-frozen')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--dataset',type=Path,default=ROOT/'data/collisions/h2o-gas-20261029-events.json');args=parser.parse_args()
    if os.sched_getaffinity(0)!={7} or os.getpriority(os.PRIO_PROCESS,0)!=15:raise RuntimeError('CPU7/nice15 required')
    lease=(OUT.parent/'pic-profile/testing-core7.lock').open('a')
    try:fcntl.flock(lease,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:raise SystemExit(75)
    if any(os.environ.get(k)!='1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']):raise RuntimeError('One thread required')
    i=1
    while (OUT/f'event-attempt-{i:03d}').exists():i+=1
    out=OUT/f'event-attempt-{i:03d}';out.mkdir();began=time.monotonic()
    def run(args,name,timeout=45):
        r=subprocess.run([str(FROZEN),*map(str,args)],capture_output=True,text=True,timeout=timeout)
        (out/name).write_text(r.stdout+r.stderr)
        if r.returncode:raise RuntimeError(f'Command failed: see {out/name}')
        return r.stdout
    base=ROOT/'data/collisions/h2o-gas-20261029.json';before=sha(base)
    assert before==json.loads(base.with_name('h2o-gas-20261029.lock.json').read_text())['dataset_sha256']
    path=args.dataset.resolve();data=json.loads(path.read_text())
    assert sha(path)==json.loads(path.with_suffix('.lock.json').read_text())['dataset_sha256']
    assert sha(OUT/'source-documents/ram-2010-water-dea.pdf')==data['sources']['ram2010']['pdf_sha256']
    if 'melton1972' in data['sources']:assert sha(OUT/'source-documents/melton-1972.pdf')==data['sources']['melton1972']['pdf_sha256']
    if 'murakami2012' in data['sources']:assert sha(OUT/'source-documents/murakami-2012-water.pdf')==data['sources']['murakami2012']['pdf_sha256']
    bindings=OUT/'source-documents/nist-OH-numerical-bindings.json'
    ram=next(m for m in data['conditional_event_models'] if m['id']=='Hminus_Ram65_factorized_window_v1')
    assert json.loads(bindings.read_text())['G_cm1_coefficients']==ram['OH_X_G_cm1_coefficients']
    flags=['/usr/bin/c++','-std=c++20','-O2','-Wall','-Wextra','-Wpedantic','-Werror','-I'+str(ROOT/'include')]
    run(['/usr/bin/c++','--version'],'frozen-compiler.txt')
    results={}
    cases=[('events','tests/test_collision_events.cpp'),('data','tests/test_collision_data.cpp')]
    if any(m['id']=='OHminus_ground_isotropic_CM_v1' for m in data['conditional_event_models']):cases.append(('two_body','tests/test_collision_two_body.cpp'))
    if any(m['id']=='A1B1_fixed_loss_dipole_ground_v1' for m in data['conditional_event_models']):cases.append(('dissociation','tests/test_collision_dissociation.cpp'))
    if any(m['id']=='outer3_BEB_Tan_parent_v1' for m in data['conditional_event_models']):cases.append(('outer_valence','tests/test_collision_outer_valence.cpp'))
    for name,source in cases:
        exe=out/('test-'+name)
        run(flags+[ROOT/source,'-o',exe],'build-'+name+'.log')
        linkage=run(['ldd',exe],'linkage-'+name+'.txt')
        assert not any(x in linkage.lower() for x in ['cuda','cudss','nvml','nvidia','libmpi'])
        results[name]=json.loads(run([exe,path],name+'-tests.json'))
        if name=='data':results['v1_regression']=json.loads(run([exe,base],'v1-regression.json'))
    assert all(r['passed'] for r in results.values()) and before==sha(base)
    binding_paths=[bindings]
    if len(cases)>3:binding_paths.append(OUT/'source-documents/hargreaves-2012-web-binding.json')
    if len(cases)>4:binding_paths.append(OUT/'source-documents/nist-fragment-ionization-bindings.json')
    summary=dict(status='conditional_event_kinematics_qualified_no_evolution',dataset=str(path),dataset_sha256=sha(path),
        v1_dataset_unchanged_sha256=before,tests=results,CPU=7,nice=15,threads=1,lease_held=True,
        GPU_used=False,PIC_or_chemistry_evolution=False,wall_s=time.monotonic()-began,
        frozen_userspace=True,frozen_runtime=str(FROZEN),frozen_runtime_receipts_sha256={str(p):sha(p) for p in [FROZEN,FROZEN.parent/'kernel.release',FROZEN.parent/'root-files.sha256']},
        source_bindings_sha256={str(p):sha(p) for p in binding_paths},dataset_rebuilt_or_modified=False,
        conditional_models=[dict(id=m['id'],support_eV=m['data_support_eV'],assumptions=m.get('assumptions',[])) for m in data['conditional_event_models']],
        limitations=['conditional near-anchor Hminus model only6.25..6.75eV','factorized laboratory angular/vibrational distributions; unresolved rotation and high-v tail','OHminus ground/isotropic endpoint has no measured energy-angle calibration','A1B1 fixed-loss/dipole/ground endpoint is conditional; broad spectrum/internal states unresolved','no complete8keV exclusive mechanism'])
    dump(out/'validation.json',summary);dump(OUT/'event-latest.json',dict(result_directory=str(out),dataset_sha256=sha(path)))
    sources=[base,path,ROOT/'include/raikou/collision_data.h',ROOT/'include/raikou/h2o_conditional_parent.h',ROOT/'include/raikou/h2o_dea_kernel.h',ROOT/'tests/test_collision_events.cpp',ROOT/'tools/build_h2o_event_extensions.py',Path(__file__)]
    if len(cases)>2:sources += [ROOT/'include/raikou/h2o_two_body_attachment.h',ROOT/'tests/test_collision_two_body.cpp',ROOT/'tools/build_h2o_attachment_extensions.py']
    if len(cases)>3:sources += [ROOT/'include/raikou/h2o_neutral_dissociation.h',ROOT/'include/raikou/h2o_common_mass_basis.h',ROOT/'tests/test_collision_dissociation.cpp',ROOT/'tools/build_h2o_dissociation_extensions.py']
    if len(cases)>4:sources += [ROOT/'include/raikou/h2o_outer_valence.h',ROOT/'tests/test_collision_outer_valence.cpp',ROOT/'tools/build_h2o_outer_valence_extensions.py']
    dump(out/'seal.json',dict(source_sha256={str(p.relative_to(ROOT)):sha(p) for p in sources},files_sha256={p.name:sha(p) for p in out.iterdir() if p.is_file()}))
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
