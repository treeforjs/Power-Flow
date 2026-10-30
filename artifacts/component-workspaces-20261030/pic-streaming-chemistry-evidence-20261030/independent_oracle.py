from pathlib import Path
import os,json,hashlib,struct,math,csv
from fractions import Fraction
P=Path(__file__).parent;O=Path('/home/treefor/sim-workspace/powerflow/outputs/pic-streaming-chemistry-evidence-20261030')
def need(x,s):
 if not x:raise RuntimeError(s)
need(os.sched_getaffinity(0)=={7} and os.getpriority(os.PRIO_PROCESS,0)==15,'CPU7/nice15 required')
checks=0
names=['selected_positive','selected_zero','selected_unsupported','unselected_positive','unselected_zero','unselected_unsupported','unqualified_no_rate']
def exact(j):return Fraction(int(j['nonnegative_units_hex'],16),1<<1074)
def obs(doc):
 out=[]
 for x in doc['conditional_models']:out.append((('conditional_model',x['channel'],x['model'],'m2'),x['selected_subset'],x['source_rate_known'],False,x['source_sigma_m2']))
 for x in doc['unapplied_observables']:out.append((('unapplied_observable',x['channel'],'',x['unit_SI']),False,x['source_value_known'],False,x['source_value_SI']))
 for x in doc['missing_rates']:out.append((('missing_rate_contract',x['channel'],'','m2'),False,x['status']=='known_unmapped_rate',x['status']=='unknown_or_unqualified',x['sigma_m2']))
 return out
def u64(n):return struct.pack('>Q',n)
def word(s):
 b=s.encode();return u64(len(b))+b
for stem in ['small','legacy-batch','duplicate']:
 inputs=json.loads((O/(stem+'-inputs.json')).read_text());receipt=json.loads((O/(stem+'-receipt.json')).read_text());header=receipt['canonical_header_json'];decoded=json.loads(header)
 need(decoded['bindings']==receipt['bindings'] and decoded['policy']==receipt['policy'],'header binding mismatch')
 need(decoded['layout']==[[r['lane'],r['id'],r['model'],r['unit']] for r in receipt['rows']],'header layout mismatch')
 h=hashlib.sha256(b'CHEM-EVIDENCE-V1'+word(header));buckets={};byte_count=0;total=Fraction()
 for q in inputs:
  E=float(q['energy_eV']);W=float(q['weight']);doc=json.loads(q['provider_canonical_json']);total+=Fraction(W)
  wire=b'Q'+u64(q['ordinal'])+word(q['id'])+word(q['context_sha256'])+struct.pack('>d',E)+struct.pack('>d',W)+word(q['provider_canonical_json']);h.update(wire);byte_count+=len(wire)
  for key,selected,known,unqualified,v in obs(doc):
   category='unqualified_no_rate' if unqualified else ('selected_unsupported' if selected else 'unselected_unsupported') if not known else ('selected_positive' if selected else 'unselected_positive') if v>0 else ('selected_zero' if selected else 'unselected_zero')
   bucket=buckets.setdefault((key,category),[]);bucket.append((q['ordinal'],E,W,v if known else None))
 need(h.hexdigest()==receipt['content_sha256'],'independent wire commitment differs');need(byte_count==receipt['encoded_query_bytes'],'encoded byte counter differs');need(exact(receipt['query_weight_exact'])==total,'exact total query weight differs');checks+=3
 for row in receipt['rows']:
  key=(row['lane'],row['id'],row['model'],row['unit']);partition=Fraction();n=0
  for category in names:
   got=row['categories'][category];raw=buckets.get((key,category),[]);n+=len(raw)
   W=sum((Fraction(r[2]) for r in raw),Fraction());partition+=W
   EW=sum((Fraction(float(r[1]*r[2])) for r in raw),Fraction())
   need(got['count']==len(raw),'category count differs');need(exact(got['weight_exact'])==W,'exact category weight differs');need(exact(got['energy_weight_product_sum'])==EW,'per-term rounded energy product sum differs');checks+=3
   known=[r for r in raw if r[3] is not None];positive=[r for r in raw if r[2]>0]
   if known:
    VW=sum((Fraction(float(r[2]*r[3])) for r in known),Fraction());need(exact(got['source_value_weight_product_sum'])==VW,'source value product sum differs');need(got['source_value_min']==min(r[3] for r in known) and got['source_value_max']==max(r[3] for r in known),'source extrema differ')
   else:need(got['source_value_weight_product_sum'] is None and got['source_value_min'] is None and got['source_value_max'] is None,'unknown source became zero')
   need(got['positive_weight_count']==len(positive),'positive weight count differs')
   need(got['energy_eV_min']==(min(r[1] for r in raw) if raw else None) and got['energy_eV_max']==(max(r[1] for r in raw) if raw else None),'finite energy extrema differ')
   need(got['positive_weight_energy_eV_min']==(min(r[1] for r in positive) if positive else None) and got['positive_weight_energy_eV_max']==(max(r[1] for r in positive) if positive else None),'positive-weight extrema differ')
   need(got['first_query_ordinal']==(raw[0][0] if raw else None) and got['last_query_ordinal']==(raw[-1][0] if raw else None),'query witnesses differ');checks+=5
  need(n==len(inputs) and partition==total,'support partition dropped/renormalized queries');checks+=1
capture=json.loads((O/'capture-receipt.json').read_text());manifest=json.loads((O/'capture-input-manifest.json').read_text());total=Fraction();count=0
with Path(manifest['primary_input']).open() as f:
 for row in csv.DictReader(f):total+=Fraction(float(row['weight']));count+=1
need(count==46726 and capture['queries']==count and exact(capture['query_weight_exact'])==total,'capture exact weight/count differs');checks+=1
energy_total=None
for row in capture['rows']:
 n=sum(x['count'] for x in row['categories'].values());weight=sum((exact(x['weight_exact']) for x in row['categories'].values()),Fraction());energy=sum((exact(x['energy_weight_product_sum']) for x in row['categories'].values()),Fraction())
 need(n==count and weight==total,'capture support partition dropped queries');need(energy_total is None or energy==energy_total,'capture energy partitions inconsistent');energy_total=energy;checks+=2
checkpoints=json.loads((O/'memory-checkpoints.json').read_text());need(all(x['retained_keys']==40 and x['retained_category_slots']==280 for x in checkpoints),'retained storage grew');checks+=1
result={'status':'passed','assertions':checks,'small_fixture_oracle':'independent Fraction support/weights/moments/extrema and SHA256 wire framing','capture_oracle':'exact input count/weight and per-key support partitions','capture_queries':count,'exact_total_weight_numerator':str(total.numerator),'exact_total_weight_denominator':str(total.denominator),'count_backend':'unchanged sealed ExactParticleCount','physics_kernel_changed':False,'new_EEDF_or_evolution':False}
(O/'independent-oracle.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
