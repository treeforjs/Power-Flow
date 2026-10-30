"""Read-only pins and independent scalar/uniform-control checks; no GPU API."""
import array,hashlib,json,math,pathlib,struct,sys
O=pathlib.Path(__file__).resolve().parents[1]
P=pathlib.Path('/home/treefor/projects/runtimes/pic-electric-exterior-20261029')
F=O.parent/'pic-profile/cuda-gather-investigation/actual-midpoint-fixture'
EPS=sys.float_info.epsilon
ELECTRON_CHARGE=1.602176634e-19
FIELDS=struct.Struct('<8d');PARTICLES=struct.Struct('<Q9d')

def sha(p):
 h=hashlib.sha256()
 with pathlib.Path(p).open('rb') as f:
  for block in iter(lambda:f.read(1048576),b''):h.update(block)
 return h.hexdigest()
def record(p):
 p=pathlib.Path(p);return dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p))
def finite(value):
 if isinstance(value,float):assert math.isfinite(value),'Nonfinite JSON value'
 elif isinstance(value,dict):
  for item in value.values():finite(item)
 elif isinstance(value,list):
  for item in value:finite(item)
def read_json(p):
 result=json.loads(pathlib.Path(p).read_text());finite(result);return result
def write_json(p,v):
 finite(v);pathlib.Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def verify_record(row):
 p=pathlib.Path(row['path']);assert p.is_file() and p.stat().st_size==row['bytes'] and sha(p)==row['sha256'],f'Pin mismatch: {p}'
def doubles(p,count=None):
 p=pathlib.Path(p);assert p.stat().st_size%8==0,f'Truncated vector: {p}'
 if count is not None:assert p.stat().st_size==8*count,f'Vector size: {p}'
 a=array.array('d');a.frombytes(p.read_bytes());assert sys.byteorder=='little'
 assert all(math.isfinite(v) for v in a),f'Nonfinite vector: {p}'
 return a
def ints(p):
 a=array.array('i');a.frombytes(pathlib.Path(p).read_bytes());assert a.itemsize==4 and sys.byteorder=='little';return a
V2_BINDING_AT_IMPORT=sha(O/'dispatcher-v2-bindings.json')
def verify_pins():
 assert sha(O/'dispatcher-v2-bindings.json')==V2_BINDING_AT_IMPORT,'Version2 manifest changed during dispatch'
 request=read_json(O/'gpu-request.json');verify_record(request['bindings']);verify_record(request['bound_executable'])
 bindings=read_json(request['bindings']['path']);count=0
 for key in ['sources','external_compile_inputs','libraries','executables','fixture','preserved_prior_evidence']:
  for item in bindings[key]:verify_record(item);count+=1
 preparation=read_json(O/'dispatcher-bindings.json')
 assert preparation['gpu_request_sha256']==sha(O/'gpu-request.json')
 for item in preparation['tools_and_guards']:verify_record(item);count+=1
 extra=read_json(O/'dispatcher-v2-bindings.json')
 assert extra['gpu_request_sha256']==sha(O/'gpu-request.json') and extra['original_dispatcher_bindings_sha256']==sha(O/'dispatcher-bindings.json')
 for row in extra['tools_and_evidence']:verify_record(row);count+=1
 return dict(passed=True,verified_records=count,gpu_request_sha256=sha(O/'gpu-request.json'),candidate_bindings_sha256=sha(P/'bindings.json'),dispatcher_version=2,dispatcher_v2_bindings_sha256=V2_BINDING_AT_IMPORT)
def check_result(folder,uniform=False):
 folder=pathlib.Path(folder);r=read_json(folder/'result.json');g=read_json(folder/'control-gates.json');w=read_json(folder/'work-and-layout.json');resources=read_json(folder/'resources.json')
 assert r['mode']=='gpu' and r['nodes']<=1010000 and r['foil_dirichlet_nodes']==3618 and r['particles']==46726
 assert r['particle_state_unchanged'] and r['no_particle_advance'] and r['same_operator_for_two_RHS']
 assert r['factorizations']==1 and r['field_solves']==2
 assert resources['allowed_logical_cpus']==[7] and resources['max_rss_KiB']*1024<=4*1024**3
 assert r['charge_error_C']==abs(r['direct_charge_C']-r['deposited_charge_C'])
 assert r['charge_error_C']<=1e-11*max(abs(r['direct_charge_C']),ELECTRON_CHARGE)
 assert g['independent_source_checked'] and g['source_gate_passed'] and 0<=r['independent_source_max_abs_error_C_per_m']<=g['source_roundoff_limit_C_per_m']
 for label in ['phi','psi']:
  assert 0<=r[label+'_residual_C_per_m']<=r[label+'_residual_limit_C_per_m']
 assert g['original_E_abs_V_m']==1e-6 and g['original_B_abs_T']==1e-13 and g['original_field_relative']==5e-12
 assert w['solution_available'] and w['work_checks_performed'] and w['work_checks_passed'] and w['field_energy_J']>=0
 assert abs(w['source_potential_work_J']-w['nodal_source_work_J'])<=w['source_work_limit_J']
 assert abs(w['electrostatic_energy_identity_defect_J'])<=w['energy_identity_limit_J']
 assert abs(w['gauss_defect_C'])<=w['gauss_limit_C'] and w['finite_vector_gauss_identity_error_C']<=w['finite_sum_roundoff_bound_C']
 computed=(abs(w['gauss_source_deposition_term_C'])+abs(w['full_operator_null_action_sum_C'])+w['free_residual_Cauchy_Schwarz_bound_C']+w['finite_sum_roundoff_bound_C'])/(1-EPS)
 assert abs(computed-w['gauss_limit_C'])<=4*math.ulp(max(computed,w['gauss_limit_C']))
 assert w['physical_particle_coordinates_unshifted'] and w['native_coordinates_unshifted'] and w['natural_Neumann_exterior'] and not w['Robin_or_sheath_smoothing']
 assert g['captured_nodal_comparison_performed']==uniform and g['capture_field_comparison_performed']==uniform and g['nodal_gate_passed']
 if not uniform:assert r['passed'] and r['captured_field_failures']==0
 doubles(folder/'phi_nodes.bin',r['nodes']);doubles(folder/'psi_nodes.bin',r['nodes'])
 assert (folder/'particle_fields.bin').stat().st_size==r['particles']*FIELDS.size
 with (folder/'particle_fields.bin').open('rb') as stream:
  for b in iter(lambda:stream.read(FIELDS.size),b''):assert len(b)==FIELDS.size and all(math.isfinite(v) for v in FIELDS.unpack(b))
 return dict(result=r,control_gates=g,work=w,resources=resources,independent_scalar_checks_passed=True)
def uniform_audit(folder,returncode):
 folder=pathlib.Path(folder);request=read_json(O/'gpu-request.json');receipt=check_result(folder,True);r=receipt['result'];g=receipt['control_gates']
 prior={pathlib.Path(x['path']).name:x for x in request['uniform_control']['accepted_fresh_control']}
 comparisons={}
 for name,item in prior.items():
  verify_record(item);comparisons[name]=dict(new=record(folder/name),previous=item,bitwise_equal=sha(folder/name)==item['sha256']);assert comparisons[name]['bitwise_equal'],f'Fresh uniform output differs: {name}'
 old=read_json(pathlib.Path(prior['phi_nodes.bin']['path']).parent/'result.json')
 assert not old['passed'] and old['captured_field_failures']==2733
 assert returncode==3 and not r['passed'] and r['captured_field_failures']==old['captured_field_failures']
 assert r['max_E_difference_V_m']==old['max_E_difference_V_m'] and r['max_B_difference_T']==0
 control=O/'uniform-control-sealed';I=ints(control/'constrained_operator_I.bin');J=ints(control/'constrained_operator_J.bin');V=doubles(control/'constrained_operator_V.bin')
 assert len(I)==r['nodes']+1 and I[0]==0 and I[-1]==len(J)==len(V)
 residuals={};nodal={}
 for label in ['phi','psi']:
  x=doubles(folder/(label+'_nodes.bin'),r['nodes']);b=doubles(control/(label+'_rhs'),r['nodes']);capture=doubles(F/(label+'.bin'),r['nodes'])
  squared=math.fsum((math.fsum(V[k]*x[J[k]] for k in range(I[i],I[i+1]))-b[i])**2 for i in range(len(b)))
  absolute=math.sqrt(squared);limit=max(1e-30,2e-13*math.sqrt(math.fsum(v*v for v in b)));assert absolute<=limit
  delta=max(abs(a-bb) for a,bb in zip(x,capture));bound=(1e-10+256*EPS)*max(abs(v) for v in capture);assert delta<=bound
  residuals[label]=dict(independent_CSR_fsum_residual_C_per_m=absolute,original_limit_C_per_m=limit,passed=True)
  nodal[label]=dict(max_difference=delta,original_global_limit=bound,passed=True)
 # Recheck source work from actual saved point values and exact raw charge.
 x=doubles(folder/'phi_nodes.bin');q=doubles(control/'charge_rhs');depth=read_json(F/'manifest.json')['depth']
 nodal_work=depth*math.fsum(a*b for a,b in zip(x,q));terms=[]
 with (F/'particles.bin').open('rb') as particles,(folder/'particle_fields.bin').open('rb') as fields:
  for _ in range(r['particles']):
   p=PARTICLES.unpack(particles.read(PARTICLES.size));field=FIELDS.unpack(fields.read(FIELDS.size));terms.append(p[7]*p[9]*field[6])
 point_work=math.fsum(terms);work_limit=1024*EPS*max(math.fsum(abs(v) for v in terms),1e-30);assert abs(point_work-nodal_work)<=work_limit
 receipt.update(status='uniform_additive_match_original_strict_failure_preserved',passed_additive_control=True,original_strict_passed=False,original_strict_failures=r['captured_field_failures'],binary_comparisons=comparisons,independent_residuals=residuals,independent_nodal_checks=nodal,independent_source_work=dict(point_work_J=point_work,nodal_work_J=nodal_work,limit_J=work_limit,passed=True))
 return receipt
