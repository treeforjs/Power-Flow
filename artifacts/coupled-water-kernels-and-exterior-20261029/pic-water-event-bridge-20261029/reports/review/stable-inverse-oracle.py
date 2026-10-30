from decimal import Decimal as D,localcontext
from pathlib import Path
import os,json
p=Path(__file__).parent
assert sorted(os.sched_getaffinity(0))==[7]
def d(x):return D.from_float(x) if isinstance(x,float) else D(x)
cases=[]
with localcontext() as ctx:
 ctx.prec=120;c=d(299792458.)
 def norm(u):return sum(x*x for x in u).sqrt()
 def gamma(u):return (1+sum(x*x for x in u)/(c*c)).sqrt()
 for line in (p/'stable-inverse.jsonl').read_text().splitlines():
  row=json.loads(line);e=list(map(d,row['e_u']));n=list(map(d,row['n_u']));dot=sum(e[i]*n[i] for i in range(3));a=dot/(c*c*(gamma(n)+1))-gamma(e);rest=[e[i]+a*n[i] for i in range(3)];rn=norm(rest);direct=[x/rn for x in rest]
  error=norm([d(row['direction'][i])-direct[i] for i in range(3)])
  cases.append({'input_eV':row['input_eV'],'boost_u_over_c':row['boost_u_over_c'],'shape':row['shape'],'direction_error':float(error)})
summary={'cases':len(cases),'max_error':max(cases,key=lambda r:r['direction_error']),'former_counterexample':[r for r in cases if r['input_eV']==14 and r['boost_u_over_c']==1e5 and r['shape']==0][0],'max_error_u_over_c_le_0p8':max(r['direction_error'] for r in cases if r['boost_u_over_c']<=.8),'status':'identity_probe_only_no_admission_or_source_change'}
(p/'stable-inverse-oracle.json').write_text(json.dumps({'summary':summary,'cases':cases},indent=2)+'\n');print(json.dumps(summary))
