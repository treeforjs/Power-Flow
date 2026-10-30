#!/usr/bin/env python3
"""Seal source, dependencies, CPU receipts and a static-only GPU request."""
import datetime,difflib,hashlib,json,pathlib,shutil
P=pathlib.Path('/home/treefor/projects/runtimes/pic-electric-exterior-20261029')
W=pathlib.Path('/home/treefor/sim-workspace/powerflow/outputs')
O=W/'pic-electric-exterior-20261029'
F=W/'pic-profile/cuda-gather-investigation/actual-midpoint-fixture'
OLD=pathlib.Path('/home/treefor/projects/runtimes/pic-static-domain-20261029')
G=pathlib.Path('/home/treefor/projects/runtimes/mfem-cuda-benchmark-20261029/cuda-libs/nvidia')
FRESH=W/'cathode-emission-audit-20261029/static-domain/9mm'
REVIEW=W/'cathode-emission-audit-20261029/static-domain-control-review/additive-control-disposition.json'
EXTENDED=W/'cathode-emission-audit-20261029/static-domain-extended'

def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def record(p):return {'path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)}
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def bind(paths):return [record(p) for p in sorted(set(paths))]

assert json.loads((O/'cpu-controls.json').read_text())['passed']
assert json.loads((O/'uniform-CSR-RHS-parity.json').read_text())['passed']
assert json.loads((O/'nested-axis-parity.json').read_text())['passed']
assert sha(P/'source/cudss_poisson.hpp')==sha(OLD/'source/cudss_poisson.hpp')
source_files=[P/'CMakeLists.txt',P/'README.txt',P/'gauss-accounting.txt',*list((P/'source').glob('*')), *list((P/'tools').glob('*.py'))]
capture=pathlib.Path('/home/treefor/projects/runtimes/pic-cuda-gather-investigation-20261029/source')
external=[capture/n for n in ['batch_snapshot.cpp','batch_snapshot.hpp','snapshot_io.hpp','gather_push_core.hpp']]
external += [pathlib.Path('/home/treefor/projects/runtimes/mfem-gpu-pic-20261029/source/include/raikou/mfem_gap_pic.h'),pathlib.Path('/home/treefor/projects/runtimes/arch-base-20261028/run-frozen'),pathlib.Path('/home/treefor/projects/runtimes/mfem-4.10-pic-20261029/install-mpi/lib/libmfem.a')]
libraries=[G/'cu12/lib/libcudss.so.0',G/'cuda_runtime/lib/libcudart.so.12',G/'cublas/lib/libcublas.so.12',G/'cublas/lib/libcublasLt.so.12']
executables=[P/'build/electric_exterior_cpu',P/'build/electric_exterior_gpu',P/'build/exterior_manufactured']
prior=[REVIEW,EXTENDED/'independent-review.json',EXTENDED/'parent-continuation-contract.json',FRESH/'result.json',FRESH/'control-gates.json',*(FRESH/n for n in ['phi_nodes.bin','psi_nodes.bin','particle_fields.bin'])]
manifest={'schema':1,'created_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'sources':bind(source_files),'external_compile_inputs':bind(external),'libraries':bind(libraries),'executables':bind(executables),'fixture':bind(list(F.iterdir())),'preserved_prior_evidence':bind(prior),'GPU_executed':False,'actual_fixture_Poisson_solves':0,'PIC_evolution':False,'unchanged_two_refinement_solver_header':True}
save(P/'bindings.json',manifest);shutil.copy2(P/'bindings.json',O/'bindings.json')
bundle=O/'source';bundle.mkdir(exist_ok=True)
for p in source_files:
 relative=p.relative_to(P);dest=bundle/relative;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest)
diff=''.join(difflib.unified_diff((OLD/'source/static_domain.cpp').read_text().splitlines(keepends=True),(P/'source/static_domain.cpp').read_text().splitlines(keepends=True),fromfile=str(OLD/'source/static_domain.cpp'),tofile=str(P/'source/static_domain.cpp')))
(O/'static-domain-source.diff').write_text(diff)
preserved=O/'preserved-initial-ledger-failure';preserved.mkdir(exist_ok=True)
for p in (O/'uniform-control').glob('*.json'):shutil.copy2(p,preserved/p.name)
save(preserved/'disposition.json',{'status':'failed_initial_added_ledger_budget_preserved','failed_budget_C':1.2809575120063138e-20,'unchanged_capture_defect_C':2.7670845851123882e-20,'cause':'new budget omitted complete FEM residual and rounded full-operator constant-mode action; see gauss-accounting.txt','original_source_solver_nodal_field_thresholds_changed':False,'initial_development_source_not_sealed':True})

env={'PIC_STATIC_DOMAIN_PARENT_GPU':'1','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','BLIS_NUM_THREADS':'1','LD_LIBRARY_PATH':':'.join(str(G/x) for x in ['cu12/lib','cublas/lib','cuda_runtime/lib'])}
def command(width,n,name):
 return ['/home/treefor/projects/runtimes/arch-base-20261028/run-frozen','env',*[f'{k}={v}' for k,v in env.items()],str(P/'build/electric_exterior_gpu'),'gpu',str(F),str(O/'parent-static-runs'/name),str(width),str(n)]
request={
 'schema':1,'status':'CPU_qualified_parent_source_reviewed_static_GPU_request_not_executed',
 'bound_executable':record(P/'build/electric_exterior_gpu'),'bindings':record(P/'bindings.json'),
 'scope':'Captured fixed-state electrostatic phi/unit-bias psi solves and original-position field gathering only; no GapPic, Boris, source injection or trajectories',
 'resources':{'CPU':[7],'nice':15,'threads':1,'exclusive_lease':str(W/'pic-profile/testing-core7.lock'),'max_nodes':1010000,'host_RSS_limit_bytes':4294967296,'per_case_timeout_s':120,'one_case_at_a_time':True,'GPU_memory_measured':False,'CPU_RSS_measured_max_KiB':348940,'dispatch_condition':'Use the existing candidate-only resource monitor; abort this candidate on memory/time/solver/source/charge/work failure. Do not impose a 4GiB virtual-address limit on CUDA UVA.'},
 'uniform_control':{
  'argv':command(9,0,'uniform9'),
  'exact_CPU_CSR_RHS_receipt':record(O/'uniform-CSR-RHS-parity.json'),
  'unchanged_solver_header':record(P/'source/cudss_poisson.hpp'),
  'original_strict_comparison_must_remain_recorded':True,
  'known_prior_strict_failures':2733,
  'expected_exit_if_only_original_strict_field_gate_fails':3,
  'continuation_requires':'Fresh phi_nodes.bin,psi_nodes.bin and particle_fields.bin must each match the previously accepted fresh-solve control bitwise; all source,charge,nodal,work and original residual checks must independently pass; record the new original strict failure unchanged and explicitly apply the prior additive continuation recipe. No blanket exit3 acceptance.',
  'accepted_fresh_control':bind([FRESH/n for n in ['phi_nodes.bin','psi_nodes.bin','particle_fields.bin']]),
  'prior_additive_disposition':record(REVIEW),
  'prior_extended_review':record(EXTENDED/'independent-review.json'),
  'additional_condition':'If any fresh-control binary differs, stop for a new causal audit; prior acceptance does not cover unexplained differences.'
 },
 'later_cases_require_separate_explicit_continuation_after_uniform_control':True,
 'stages':[
  {'purpose':'Match prior uniform13/17/25mm fixed-state results and measure graded exterior discretization, without claiming equality','cases':[{'argv':command(w,n,f'graded-{w}-{n}'),'width_mm':w,'margin':n} for w,n in [(13,32),(17,40),(25,48)]]},
  {'purpose':'Move far boundary on a bitwise nested fixed1.05 mesh; every existing interior cell stays identical','cases':[{'argv':command(0,n,f'nested-105-{n}'),'fixed_growth':1.05,'margin':n} for n in [40,60,80,96,110,120]]},
  {'purpose':'Separate outer interpolation/discretization refinement at fixed80mm full extent','cases':[{'argv':command(80,n,f'width80-{n}'),'width_mm':80,'margin':n} for n in [80,120,160]]}
 ],
 'mandatory_guards':['FP64 original epsilon0/operator and same physical finite-foil Dirichlet ownership','unchanged point charge/depth and original source charge conservation','independent physical source assembly check and bitwise core-source mapping','original phi and psi max(1e-30,2e-13*norm(rhs)) residual guards','default cuDSS FP64 with exactly two iterative-refinement steps','source-potential work identity on actual source-work scale','driven energy identity and complete free-residual/null-mode Gauss accounting, with first ledger failure retained','finite particle fields, unchanged particle bytes, factorization1 and solves2','B bitwise invariant across all cases; no shifted coordinates or extrapolated native B'],
 'required_reporting':['Actual axis extents, node/nnz counts and peak host/GPU memory','All/above-upper-back/gap/near-edge field changes at the same captured particle positions; component maxima and weighted norms','Separate E perpendicular/parallel relative to unchanged native B; do not advance approximate guiding centres','Compare extent refinement and fixed-extent outer refinement separately','Report actual raw charge/self-field source norms; driven-bias residual alone cannot establish tiny space-charge accuracy','No declaration of physical/infinite-domain convergence from finite successive differences; no production/PIC continuation'],
 'GPU_dispatch_in_this_preparation':False,'physics_or_solver_threshold_changes':False,'production_actions':False}
save(P/'gpu-request.json',request);shutil.copy2(P/'gpu-request.json',O/'gpu-request.json')
save(O/'handoff.json',{'status':'READY_CPU_controls_and_reviewed_static_request','prefix':str(P),'executables':bind(executables),'source_bundle':str(bundle),'CPU_controls':str(O/'cpu-controls.json'),'uniform_exact_parity':str(O/'uniform-CSR-RHS-parity.json'),'manufactured':str(O/'manufactured-controls.json'),'nested_axis_parity':str(O/'nested-axis-parity.json'),'GPU_request':str(O/'gpu-request.json'),'bindings':str(O/'bindings.json'),'GPU_executed':False,'actual_fixture_solves':0,'PIC_evolution':False,'production_actions':False,'evolving_driver_integration_qualified':False})
(O/'READY').write_text('CPU-qualified static exterior; pinned reviewed GPU request only; no GPU or PIC execution.\n')
print('Pinned source/handoff READY; GPU not executed')
