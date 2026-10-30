#!/usr/bin/env python3
"""Seal prepared static-domain source, bound fixture and CPU evidence; no GPU."""
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import tarfile

root=Path(__file__).resolve().parents[1]
out=Path('/home/treefor/sim-workspace/powerflow/outputs/pic-profile/em-field-consistency-audit/static-domain-preparation')
fixture=out.parents[1]/'cuda-gather-investigation'/'actual-midpoint-fixture'
fixture=Path('/home/treefor/sim-workspace/powerflow/outputs/pic-profile/cuda-gather-investigation/actual-midpoint-fixture')
control=out/'cpu-qualified'
gpu_base=Path('/home/treefor/projects/runtimes/mfem-cuda-benchmark-20261029/cuda-libs/nvidia')
library_path=':'.join(str(gpu_base/p) for p in ['cu12/lib','cublas/lib','cuda_runtime/lib'])
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1048576),b''):h.update(block)
    return h.hexdigest()
expected=json.loads((fixture/'sha256.json').read_text())
assert all(sha(fixture/n)==h for n,h in expected.items())
result=json.loads((control/'result.json').read_text());gates=json.loads((control/'control-gates.json').read_text())
assert result['passed'] and result['field_solves']==0 and result['factorizations']==0
assert result['max_E_difference_V_m']==result['max_B_difference_T']==0
assert gates['source_gate_passed'] and gates['nodal_gate_passed']
files=set(p for p in fixture.iterdir() if p.is_file())
files.update([root/'CMakeLists.txt',root/'README.txt',root/'build/static_domain_cpu',root/'build/static_domain_gpu'])
files.update(p for folder in ['source','tools'] for p in (root/folder).iterdir() if p.is_file())
files.add(Path('/home/treefor/projects/runtimes/mfem-gpu-pic-20261029/source/src/mfem_gap_pic.cpp'))
files.add(Path('/home/treefor/projects/runtimes/mfem-gpu-pic-20261029/source/src/cudss_poisson.hpp'))
files.add(Path('/home/treefor/projects/runtimes/mfem-gpu-pic-20261029/source/src/mfem_pic_distributed.hpp'))
files.add(Path('/home/treefor/sim-workspace/powerflow/outputs/pic-profile/cpu-independent-review/checkpoint-validation/cli_pass.json'))
for target in ['static_domain_cpu','static_domain_gpu']:
    dependency_text=subprocess.check_output(['ninja','-C',str(root/'build'),'-t','deps',f'CMakeFiles/{target}.dir/source/static_domain.cpp.o'],text=True)
    for line in dependency_text.splitlines():
        if line.startswith('    /'):files.add(Path(line.strip()))
    dynamic=subprocess.check_output(['ldd',str(root/'build'/target)],text=True,env={**os.environ,'LD_LIBRARY_PATH':library_path})
    assert 'not found' not in dynamic,dynamic
    (out/f'{target}-dynamic-dependencies.txt').write_text(dynamic)
    for match in re.findall(r'(/[^\s]+)',dynamic):
        path=Path(match)
        if path.is_file():files.add(path)
files.add(Path('/usr/bin/c++'))
frontend=subprocess.check_output(['/usr/bin/c++','-print-prog-name=cc1plus'],text=True).strip()
if Path(frontend).is_file():files.add(Path(frontend))
for path in [
 '/home/treefor/projects/runtimes/mfem-4.10-pic-20261029/install-mpi/lib/libmfem.a',
 '/home/treefor/projects/runtimes/mfem-4.10-pic-20261029/tpl/hypre/lib/libHYPRE.a',
 '/home/treefor/projects/runtimes/mfem-4.10-pic-20261029/tpl/metis/lib/libmetis.a',
 '/home/treefor/projects/runtimes/mfem-4.10-pic-20261029/tpl/GKlib/lib/libGKlib.a']:
    files.add(Path(path))
bindings=dict(schema='static-domain-input-bindings-v1',fixture=str(fixture),cpu_control=str(control),
    sha256={str(p):sha(p) for p in sorted(files)})
(root/'input-bindings.json').write_text(json.dumps(bindings,indent=2)+'\n')
cases=[]
for mm,nodes in [(9,130321),(13,271441),(17,463761),(25,1002001)]:
    cases.append(dict(full_width_mm=mm,nodes=nodes,cell_size_um=25,
        output=str(out/f'parent-gpu-{mm}mm'),
        argv=[str(root/'build/static_domain_gpu'),'gpu',str(fixture),str(out/f'parent-gpu-{mm}mm'),str(mm)],
        prerequisite='Bound CPU prerequisite and production hold' if mm==9 else '9mm GPU control passed all captured phi/psi/E/B/source/residual gates with the same input bindings'))
request=dict(schema='parent-static-domain-gpu-request-v1',status='prepared_not_executed',
    runtime=str(root),input_bindings=str(root/'input-bindings.json'),input_bindings_sha256=sha(root/'input-bindings.json'),
    verify_argv=['python',str(root/'tools/verify_inputs.py')],cpu_control=str(control),
    cpu_control_result=result,cpu_control_gates=gates,
    scope='Static field-only sensitivity, not a continuous restart, particle history or EEDF convergence run',
    production='Keep the full run and cgroup processes frozen; no resume or production edits',
    execution_owner='parent_only',environment={'PIC_STATIC_DOMAIN_PARENT_GPU':'1','LD_LIBRARY_PATH':library_path},
    serial_locks=['/home/treefor/sim-workspace/powerflow/outputs/pic-profile/testing-core.lock','/home/treefor/sim-workspace/powerflow/outputs/pic-profile/testing-gpu.lock'],
    cpu_affinity=[6],nice=15,jobs=1,cases=cases,
    limits=dict(max_nodes=1010000,max_particles=46726,host_RSS_bytes=4*1024**3,GPU_used_memory_increase_bytes=8*1024**3,wall_seconds_per_case=900,output_bytes_per_case=64*1024**2,factorizations_per_case=1,RHS_solves_per_case=2,enforcement='Parent monitors and aborts diagnostic case on a cap breach; limits are ceilings, not measured predictions'),
    validation_order=['Verify bound hashes and production hold before dispatch','Run only9mm GPU control first','Require result.passed and all control-gates plus same particle/input hashes','Only then consider13/17/25mm sequentially','Verify bound inputs unchanged after every case','Report field-only sensitivities, absolute errors and resource use; no accepted trajectory or EEDF claim'],
    graded_mesh=dict(status='deferred',reason='No graded source/gather mapping independently qualified; uniform25um only'),
    comparison_policy=dict(E_abs_V_m=1e-6,B_abs_T=1e-13,field_relative=5e-12,nodal_global_relative=1e-10,nodal_roundoff_factor=256,source_roundoff_factor=1024,residual='max(1e-30,2e-13*||b||2)',iterative_refinements=2),
    output_layout=dict(particle_fields='64B record:FP64 E[3],B[3],phi,psi; original particle order/IDs bound by particles.bin hash',nodal_fields='Raw FP64 MFEM nodal phi/psi, uniform mesh/DOF ordering from source',not_restart_checkpoint=True))
(out/'parent-request.json').write_text(json.dumps(request,indent=2)+'\n')
(out/'cpu-control-summary.json').write_text(json.dumps(dict(status='compiled_CPU_qualified_GPU_pending',source_prefix=str(root),CPU=result,gates=gates,resources=json.loads((control/'resources.json').read_text()),bound_files=len(files),cuDSS_header_unchanged=sha(root/'source/cudss_poisson.hpp')==sha('/home/treefor/projects/runtimes/mfem-gpu-pic-20261029/source/src/cudss_poisson.hpp'),no_GPU_execution=True),indent=2)+'\n')
(out/'cpu-preflight-v1/invalidated.json').write_text(json.dumps(dict(status='invalidated_scaffold_receipt',reason='Superseded:raw charge comparison previously used an MFEM-mutated RHS and unit RHS lifetime was not retained. Final cpu-qualified fixes both; do not use v1.'),indent=2)+'\n')
with tarfile.open(out/'source-and-contracts.tar.gz','w:gz') as archive:
    for name in ['CMakeLists.txt','README.txt','input-bindings.json','source','tools']:
        archive.add(root/name,arcname='pic-static-domain-20261029/'+name)
manifest={str(p.relative_to(out)):dict(bytes=p.stat().st_size,sha256=sha(p)) for p in out.rglob('*') if p.is_file() and p.name!='output-manifest.json'}
(out/'output-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(dict(status='compiled_CPU_qualified_GPU_pending',bound_files=len(files),parent_request=str(out/'parent-request.json'),source_package=str(out/'source-and-contracts.tar.gz'))))
