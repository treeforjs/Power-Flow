"""Private source-bound candidate table; no production or model promotion."""
import copy,hashlib,importlib.util,json,math,sys
from decimal import Decimal as D,Context,ROUND_HALF_EVEN,localcontext,InvalidOperation,DivisionByZero,Overflow
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE=Path('/home/treefor/projects/runtimes/pic-water-chemistry-v3.1.1-20261031')
REF=Path('/home/treefor/projects/runtimes/water-fragment-source-inventory-20261101/references')
CTX=Context(prec=70,rounding=ROUND_HALF_EVEN,Emin=-999999,Emax=999999,flags=[],traps=[InvalidOperation,DivisionByZero,Overflow])
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def module(name,path):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m
def write(path,data):Path(path).write_text(json.dumps(data,sort_keys=True,indent=2)+'\n')
def main():
 pins=json.loads((ROOT/'INPUTS.json').read_text())
 for p,h in pins['dependencies'].items():assert sha(p)==h,p
 original=json.loads((BASE/'data/tables/chemistry-v3.1-spec.json').read_text())
 registry=json.loads((BASE/'data/tables/chemistry-v3.1-manifest.json').read_text());sp={s['name']:s for s in registry['species']}
 builder=module('candidate_base_builder',BASE/'tools/build_tables.py')
 h2=module('candidate_h2',Path('/home/treefor/projects/runtimes/h2-dissociation-primary-v1.1-20261101/src/h2_primary.py'))
 oh=module('candidate_oh',Path('/home/treefor/projects/runtimes/oh-primary-supplement-v1.2-20261101/src/oh_primary.py'))
 hdata=h2._parse_source(h2._read_pinned(REF/'MCCC-el-H2-b3Su_DE.X1Sg_vi=0.txt',h2.SOURCE_SHA256,65536))
 ochannel=oh.load_channel(REF/'PSSTab264csupp_e-OH-data.xlsx','excitation_C')
 assert ochannel.label=='X 2\\Pi --> A 4\\Sigma-' and ochannel.value_unit=='bohr^2' and ochannel.energy_unit=='eV'
 with localcontext(CTX):
  m=D.from_float(sp['e']['mass_eV']);a0=D.from_float(registry['constants']['bohr_m'])
  qh=(sp['H']['ref_eV']+sp['H']['ref_eV'])-sp['H2']['ref_eV'];mh=D.from_float(sp['H2']['mass_eV']);Q=D.from_float(qh)
  ht=Q*(2*(m+mh)+Q)/(2*mh)
  hrows=[(e,s) for e,s in hdata.rows if e>ht]
  assert hrows[0][0]==D('4.5') and hrows[-1][0]==D(150) and all(s>0 for e,s in hrows)
  mo=D.from_float(sp['OH']['mass_eV']);loss=D.from_float(7.36);ot=loss*(2*(m+mo)+loss)/(2*mo)
  orows=[(e,s) for e,s in ochannel.rows if e>=D('7.3605')]
  assert orows[0][0]>ot and orows[-1][0]==D(10) and all(s>0 for e,s in orows)
  def channel(cid,rows,doi,locator):return {'id':cid,'energy_eV':[float(e) for e,s in rows],'values':[float(s*a0*a0) for e,s in rows],'ordinate_scale_to_SI':1.0,'source':doi,'source_locator':locator}
  bank={'id':'primary-h2-oh-fragment-candidate-v1','channels':[channel('H2:b3Su:integral-rate-only',hrows,'10.1016/j.adt.2020.101361','v_i=0 DE; exact source knots from first positive point above registry threshold; no differential closure'),channel('OH:a4Sigma:vertical-iso',orows,'10.1088/1361-6595/ab364c','Excitation C; fixed 7.36 eV vertical loss and ground O+H diagnostic closure')]}
 write(ROOT/'data/private-primary-rates.json',bank)
 spec=copy.deepcopy(original);spec['id']='water-fragment-table-candidate-v1.2-based-on-v3.1.1'
 bottom,top=spec['energy_bottom_eV'],spec['energy_top_eV']
 def missing(target,label,kind,products,loss=0):
  p={'id':f'e + {target} {label} (no event-qualified data)','kind':kind,'tag':'model_derived','angular':'ISO','products':products,'known_zero_below_eV':0.0,'placeholder':True,'sigma':{'type':'piecewise','segments':[{'source':{'type':'none','lo':bottom,'hi':top,'locator':label},'class':'UNRESOLVED','tag':'model_derived'}]},'notes':'Explicit missing channel, never silently identified with known zero. No numerical upper bound or event admission.'}
  if kind in ('DISSOC','EXCITE'):p['loss_eV']=loss
  return p
 def products(target):return [['e','SURVIVOR'],[target,'GAS']]
 qo=sp['H']['ref_eV']+sp['O']['ref_eV']-sp['OH']['ref_eV']
 op={'id':'e + OH -> e + H + O via a4Sigma (diagnostic vertical-ISO model)','kind':'DISSOC','loss_eV':7.36,'kinetic_fraction':1.0,'sink':'NONE','angular':'ISO','tag':'model_derived','products':[['e','SURVIVOR'],['H','RESERVOIR'],['O','RESERVOIR']],
     'sigma':{'type':'piecewise','segments':[{'source':{'type':'raikou','file':str(ROOT/'data/private-primary-rates.json'),'channel':'OH:a4Sigma:vertical-iso'},'lo':float(orows[0][0]),'hi':10.0,'class':'MODEL','tag':'model_derived'}]},
     'notes':'Calculated source integral; diagnostic differential closure uses fixed 7.36 eV vertical loss, isotropic electron pair-CM and isotropic two-body ground O+H breakup. No unspecified internal energy. Source values below the first actual knot above finite-target threshold are not admitted. Other OH channels remain unresolved.'}
 oprocess=[op,missing('OH','elastic/momentum transfer','ELASTIC',products('OH')),
           missing('OH','rotational/vibrational/other electronic excitation','EXCITE',products('OH')),
           missing('OH','other neutral dissociation','DISSOC',[['e','SURVIVOR'],['H','RESERVOIR'],['O','RESERVOIR']],qo),
           missing('OH','single and dissociative ionization','IONIZE',[['e','SURVIVOR'],['e','BIRTH'],['OH+','BIRTH']]),
           missing('OH','dissociative attachment','CAPTURE',[['H','RESERVOIR'],['O-','BIRTH']]),
           missing('OH','radiative attachment','CAPTURE',[['OH-','BIRTH']])]
 assert all(g['target'] not in ('OH','H2') or g['projectile']!='e' for g in spec['groups'])
 spec['groups'] += [{'projectile':'e','target':'OH','processes':oprocess}]
 write(ROOT/'data/candidate-spec.json',spec)
 builder.build(str(ROOT/'data/candidate-spec.json'),str(ROOT/'data/candidate'))
 policy=module('candidate_policy',BASE/'tools/make_coverage_policy.py');policy.main(str(ROOT/'data/candidate-manifest.json'),str(ROOT/'data/candidate-coverage-policy.json'))
 out=json.loads((ROOT/'data/candidate-manifest.json').read_text());pol=json.loads((ROOT/'data/candidate-coverage-policy.json').read_text())
 assert out['species']==registry['species'] and pol['admitted_unbounded_regions']==[]
 assert out['groups'][:-1]==registry['groups'] and out['processes'][:len(registry['processes'])]==registry['processes']
 new=out['processes'][len(registry['processes']):]
 assert all(not p['source_values_below_threshold_set_to_zero'] for p in new)
 assert all(not r['bounded'] for r in out['region_info'] if r['group'] in ('e+H2','e+OH'))
 contract={'schema':'primary-fragment-table-candidate-v1.2','base_blob_sha256':registry['blob_sha256'],'candidate_blob_sha256':out['blob_sha256'],'candidate_manifest_sha256':sha(ROOT/'data/candidate-manifest.json'),'candidate_policy_sha256':sha(ROOT/'data/candidate-coverage-policy.json'),'private_rate_bank_sha256':sha(ROOT/'data/private-primary-rates.json'),
   'abi':'unchanged PFCHEMV3 version 1; SpeciesRec64/ProcessRec512; 22 original species; no native header edits',
   'OH':{'source_window_eV':[float(orows[0][0]),10.0],'source_knots':len(orows),'source_rate':'calculated','event_model':'diagnostic fixed vertical loss / isotropic electron and fragments','event_implementation_qualified_component_only':True,'loss_eV':7.36,'reference_cost_eV':qo,'target_rest_threshold_eV':str(ot),'ground_fragment_KER_eV':7.36-qo},
   'H2':{'source_window_eV':[4.5,150.0],'source_knots':len(hrows),'reference_cost_eV':qh,'target_rest_threshold_eV':str(ht),'source_rate':'calculated integral over continuum pseudostates','rate_bank_only':True,'native_group_absent':True,'native_sampling_row':'none; positive rate bank is not native event support','event_implementation_qualified_component_only':False,'reason':'no continuum-resolved fragment energy or declared differential closure; loss=Q must not silently force zero relative breakup energy'},
   'new_missing_processes':[p['id'] for p in new if 'no event-qualified data' in p['id']],
   'new_group_regions':sum(r['group'] in ('e+H2','e+OH') for r in out['region_info']),
   'strict_populated_groups_expected':'H2 refused by authentic absent-group scope guard; OH refused by missing-channel reachable support including its known source window',
   'production_qualified':False,'physical_model_validated':False,'network_complete':False,'native_adopted':False}
 write(ROOT/'data/contract.json',contract)
if __name__=='__main__':main()
