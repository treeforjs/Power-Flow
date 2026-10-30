from pathlib import Path
from decimal import Decimal,localcontext
import json,csv,math,sys
p=Path(__file__).parent
cases=json.loads((p/'cases.json').read_text());actual={r['label']:r for r in csv.DictReader((p/'native-kinematics.csv').open())}
D=Decimal.from_float;c=D(299792458.);e=D(1.602176634e-19);m=D(9.1093837139e-31)
rows=[]
with localcontext() as ctx:
 ctx.prec=900
 for x in cases:
  a=[D(v)/c for v in x['u']];b=[D(v)/c for v in x['v']]
  gp=(1+sum(v*v for v in a)).sqrt();gn=(1+sum(v*v for v in b)).sqrt()
  # Independent high-precision direct four-vector invariant, not the implementation's rationalization.
  dg=gp*gn-sum(v*w for v,w in zip(a,b))-1
  if x['u']==x['v']:dg=Decimal(0) # exact mathematical identity, avoids a finite-precision sqrt residue
  assert dg>=0,(x['label'],str(dg))
  E=m*c*c*dg/e;F=c*(dg*(dg+2)).sqrt()/(gp*gn)
  result=actual[x['label']];expected_overflow=E>D(sys.float_info.max) or F>D(sys.float_info.max)
  r={'label':x['label'],'reference_energy_eV':str(E),'reference_moller_m_s':str(F),'native_status':result['status'],'expected_overflow_rejection':expected_overflow}
  if result['status']=='ok':
   for key,ref,col in [('energy',E,'energy_eV'),('flux',F,'moller_m_s')]:
    v=float(result[col]);r[key+'_native']=v;r[key+'_underflow_accepted']=bool(ref>0 and v==0)
    r[key+'_relative_error']=float(abs(D(v)-ref)/ref) if ref else (0 if v==0 else math.inf)
    r[key+'_absolute_error']=str(abs(D(v)-ref))
  rows.append(r)
summary={'oracle':'900 decimal digits; direct gamma_p gamma_n - u_p dot u_n/c^2 -1 on exact binary64 inputs','cases':len(rows),'ordinary_overflow_rejections':[r['label'] for r in rows if r['expected_overflow_rejection'] and r['native_status']=='rejected'],'unexpected_rejections':[r['label'] for r in rows if not r['expected_overflow_rejection'] and r['native_status']!='ok'],'positive_energy_underflow_accepted':[r['label'] for r in rows if r.get('energy_underflow_accepted')],'positive_flux_underflow_accepted':[r['label'] for r in rows if r.get('flux_underflow_accepted')],'large_relative_errors':[{'label':r['label'],'energy_relative_error':r.get('energy_relative_error'),'flux_relative_error':r.get('flux_relative_error'),'energy_native':r.get('energy_native'),'flux_native':r.get('flux_native')} for r in rows if r.get('energy_relative_error',0)>2e-12 or r.get('flux_relative_error',0)>2e-12]}
(p/'high-precision-oracle.json').write_text(json.dumps({'summary':summary,'results':rows},indent=2)+'\n')
print(json.dumps(summary,indent=2))
