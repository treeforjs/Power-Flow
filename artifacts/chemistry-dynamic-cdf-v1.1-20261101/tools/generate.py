"""Version2 blob extension with an explicit dynamic-loss/CDF process descriptor."""
import copy,hashlib,importlib.util,json,math,struct,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE=Path('/home/treefor/projects/runtimes/pic-water-chemistry-v3.1.1-20261031')
CAND=Path('/home/treefor/projects/runtimes/water-fragment-table-candidate-v1.2-20261101')
MODEL=Path('/home/treefor/projects/runtimes/h2-b3sigma-cdf-diagnostic-v1-20261101')
REF=Path('/home/treefor/projects/runtimes/water-fragment-source-inventory-20261101/references')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def module(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m
def write(p,x):Path(p).write_text(json.dumps(x,sort_keys=True,indent=2)+'\n')
def read_sections(raw):
    assert raw[:8]==b'PFCHEMV3';version,count=struct.unpack_from('<II',raw,8);assert version==1
    sections={}
    for i in range(count):
        name,off,n=struct.unpack_from('<24sQQ',raw,16+40*i);name=name.rstrip(b'\0').decode();assert name not in sections
        sections[name]=raw[off:off+n]
    return sections
def encode(sections):
    names=list(sections);pos=(16+40*len(names)+63)//64*64;toc=[];pieces=[]
    for name in names:
        b=sections[name];toc.append(struct.pack('<24sQQ',name.encode(),pos,len(b)));pad=(-len(b))%64;pieces.append(b+b'\0'*pad);pos+=len(b)+pad
    header=b'PFCHEMV3'+struct.pack('<II',2,len(names))+b''.join(toc);return header+b'\0'*((-len(header))%64)+b''.join(pieces)
def main():
    pins=json.loads((ROOT/'INPUTS.json').read_text())
    for p,h in pins['dependencies'].items():assert sha(p)==h,p
    builder=module('base_builder',BASE/'tools/build_tables.py');policy=module('strict_policy',BASE/'tools/make_coverage_policy.py')
    cdf=module('primary_cdf',MODEL/'src/primary_cdf.py').load(REF/'sartori2024-table-I.txt')
    model=json.loads((MODEL/'build/model.json').read_text());assert model['rng_layout_version']==2 and model['omitted_probability']=='0.002'
    spec=json.loads((CAND/'data/candidate-spec.json').read_text());spec['id']='chemistry-dynamic-cdf-v1'
    Q=model['reference_cost_eV'];thr=float(model['minimum_event_target_rest_threshold_eV'])
    source={'type':'raikou','file':str(CAND/'data/private-primary-rates.json'),'channel':'H2:b3Su:integral-rate-only'}
    # The unchanged builder creates grids/majorants/coverage. Its temporary fixed-
    # loss record is never published or executed: only the final v2 blob is valid.
    h={'id':'e + H2(v=0) -> e + H + H via b3Sigma source-CDF diagnostic','kind':'DISSOC','loss_eV':Q,'kinetic_fraction':1.,'sink':'NONE','angular':'ISO','tag':'model_derived','products':[['e','SURVIVOR'],['H','RESERVOIR'],['H','RESERVOIR']],
       'sigma':{'type':'piecewise','segments':[{'source':source,'lo':4.5,'hi':thr,'class':'UNRESOLVED','tag':'model_derived'},{'source':source,'lo':thr,'hi':150.,'class':'MODEL','tag':'model_derived'}]},
       'notes':'Positive unchanged MCCC integral; event closure only central source quantiles0.001..0.999 conditioned on finite-target phase space. Missing low-energy closure remains UNRESOLVED, never zero. Dynamic descriptor replaces temporary builder loss.'}
    def missing(label,kind,products,loss=0.):
        p={'id':'e + H2 '+label+' (no event-qualified data)','kind':kind,'tag':'model_derived','angular':'ISO','products':products,'known_zero_below_eV':0.,'placeholder':True,'sigma':{'type':'piecewise','segments':[{'source':{'type':'none','lo':spec['energy_bottom_eV'],'hi':spec['energy_top_eV'],'locator':label},'class':'UNRESOLVED','tag':'model_derived'}]}}
        if kind in ('EXCITE','DISSOC'):p['loss_eV']=loss
        return p
    return_target=[['e','SURVIVOR'],['H2','GAS']]
    ps=[h,missing('elastic/momentum transfer','ELASTIC',return_target),missing('rotational/vibrational/electronic excitation','EXCITE',return_target),missing('other neutral dissociation','DISSOC',[['e','SURVIVOR'],['H','RESERVOIR'],['H','RESERVOIR']],Q),missing('single and dissociative ionization','IONIZE',[['e','SURVIVOR'],['e','BIRTH'],['H2+','BIRTH']]),missing('dissociative attachment','CAPTURE',[['H-','BIRTH'],['H','RESERVOIR']]),missing('other radiative energy-loss processes','EXCITE',return_target)]
    spec['groups'].append({'projectile':'e','target':'H2','processes':ps})
    write(ROOT/'data/internal-builder-spec.json',spec);builder.build(str(ROOT/'data/internal-builder-spec.json'),str(ROOT/'data/internal-builder-only'))
    raw=(ROOT/'data/internal-builder-only.blob').read_bytes();sections=read_sections(raw)
    manifest=json.loads((ROOT/'data/internal-builder-only-manifest.json').read_text());idx=next(i for i,p in enumerate(manifest['processes']) if p['id']==h['id'])
    proc=bytearray(sections['procs']);struct.pack_into('<I',proc,idx*512,32);struct.pack_into('<I',proc,idx*512+24,2);struct.pack_into('<I',proc,idx*512+88,0);struct.pack_into('<d',proc,idx*512+96,0.);sections['procs']=bytes(proc)
    by={s['name']:i for i,s in enumerate(manifest['species'])}
    desc=struct.pack('<12I8d',1,1,by['e'],by['H2'],by['H'],0,2,96,0,len(cdf),0,0,4.5,150.,.001,.999,.002,Q,thr,0.)
    desc+=bytes.fromhex(model['source_CDF_sha256'])+bytes.fromhex(model['source_rate_sha256'])+bytes.fromhex(sha(MODEL/'build/model.json'))+bytes(48);assert len(desc)==256
    sections['cdf_desc']=desc;sections['cdf_quantile']=struct.pack('<100d',*(float(q) for q,e in cdf));sections['cdf_fragment_eV']=struct.pack('<100d',*(float(e) for q,e in cdf))
    final=encode(sections);(ROOT/'data/chemistry-cdf.blob').write_bytes(final)
    manifest['blob_sha256']=hashlib.sha256(final).hexdigest();manifest['blob_version']=2;manifest['process_schema']='pfc-process-v2-cdf32';manifest['CDF_descriptor_schema']='pfc-cdf-h2-central-v1';manifest['source_CDF_model_sha256']=sha(MODEL/'build/model.json')
    mp=manifest['processes'][idx];mp['kind']='DISSOC_CDF';mp['loss_eV']=None;mp['dynamic_descriptor_index']=0;mp['event_support_threshold_eV']=thr;mp['rate_kinematic_threshold_eV']=mp.pop('kinematic_threshold_eV');mp['whole_channel_supported']=False
    assert not mp['source_values_below_threshold_set_to_zero'];assert all(not x['bounded'] for x in manifest['region_info'] if x['group']=='e+H2')
    write(ROOT/'data/chemistry-cdf-manifest.json',manifest);policy.main(str(ROOT/'data/chemistry-cdf-manifest.json'),str(ROOT/'data/chemistry-cdf-coverage-policy.json'))
    pol=json.loads((ROOT/'data/chemistry-cdf-coverage-policy.json').read_text());assert pol['admitted_unbounded_regions']==[]
    write(ROOT/'data/contract.json',{'schema':'chemistry-dynamic-cdf-v1','blob_sha256':manifest['blob_sha256'],'manifest_sha256':sha(ROOT/'data/chemistry-cdf-manifest.json'),'policy_sha256':sha(ROOT/'data/chemistry-cdf-coverage-policy.json'),'model_sha256':sha(MODEL/'build/model.json'),'process_kind':32,'process_index':idx,'process_record_bytes':512,'descriptor_record_bytes':256,'blob_version':2,'descriptor_version':1,'RNG_layout_version':2,'conditional_event_threshold_eV':thr,'quantiles':[.001,.999],'omitted_probability':.002,'target_v':0,'rate_window_eV':[4.5,150.],'unsupported_selected_event_status':8,'ordinary_fixed_loss_semantics_preserved':True,'production_qualified':False,'physical_model_validated':False,'network_complete':False,'native_adopted':False})
if __name__=='__main__':main()
