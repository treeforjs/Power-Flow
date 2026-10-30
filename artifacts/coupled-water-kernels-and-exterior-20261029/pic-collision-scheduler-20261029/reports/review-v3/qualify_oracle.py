from pathlib import Path
from decimal import Decimal,localcontext
import json,math,os
p=Path(__file__).parent;assert sorted(os.sched_getaffinity(0))==[7]
d=json.loads((p/'high-precision-oracle.json').read_text());ordinary=[]
for r in d['results']:
 if r['native_status']!='ok':continue
 for kind,key in [('energy','reference_energy_eV'),('flux','reference_moller_m_s')]:
  ref=float(Decimal(r[key]));r[kind+'_correctly_rounded_binary64']=ref
  r[kind+'_ULP_difference']=abs(r[kind+'_native']-ref)/math.ulp(ref) if math.isfinite(ref) else None
 if 'ulp' not in r['label'] and not r['energy_underflow_accepted'] and not r['flux_underflow_accepted']:ordinary.append(r)
admitted=[r for r in d['results'] if r['native_status']=='ok']
d['summary']['all_admitted_cases']=len(admitted)
d['summary']['all_admitted_max_energy_ULP_difference']=max(r['energy_ULP_difference'] for r in admitted)
d['summary']['all_admitted_max_flux_ULP_difference']=max(r['flux_ULP_difference'] for r in admitted)
d['summary']['ordinary_representable_cases']=len(ordinary)
d['summary']['ordinary_max_energy_ULP_difference']=max(r['energy_ULP_difference'] for r in ordinary)
d['summary']['ordinary_max_flux_ULP_difference']=max(r['flux_ULP_difference'] for r in ordinary)
d['summary']['subnormal_rounding_note']='A large relative error can be unavoidable binary64 output quantization; ULP comparisons separate that from internal precision loss.'
with localcontext() as ctx:
 ctx.prec=900
 row=next(r for r in d['results'] if r['label']=='boost8keV_7')
 rate=Decimal.from_float(1e200)*Decimal.from_float(1e-20)*Decimal(row['reference_moller_m_s'])
 d['summary']['underflow_rate_counterexample']={'weights':1e200,'volume_m3':1,'sigma_m2':1e-20,'reference_rate_s':str(rate),'reference_rate_float':float(rate),'v2_native_candidates':0,'v3_underflow_rejected':row['native_status']=='rejected','positive_rate_representable':float(rate)>0}
(p/'high-precision-oracle.json').write_text(json.dumps(d,indent=2)+'\n')
print({k:v for k,v in d['summary'].items() if 'large_relative' not in k and k!='underflow_rate_counterexample'})
print('reference underflow-amplified rate',float(rate))
