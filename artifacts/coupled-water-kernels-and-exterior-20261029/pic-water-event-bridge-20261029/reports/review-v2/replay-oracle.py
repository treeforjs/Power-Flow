from pathlib import Path
from decimal import Decimal as D, localcontext
import json,math,os
p=Path(__file__).parent
assert sorted(os.sched_getaffinity(0))==[7]
def d(x):return D.from_float(x) if isinstance(x,float) else D(x)
def gamma(u):return (1+sum(x*x for x in u)/(c*c)).sqrt()
def norm(u):return sum(x*x for x in u).sqrt()
def f(x):return float(x)
cases=[];mutations=[]
with localcontext() as ctx:
 ctx.prec=120
 c=d(299792458.);qe=d(1.602176634e-19);me=d(9.1093837139e-31);mn=d(2.9915e-26)
 for line in (p/'replay-probes.jsonl').read_text().splitlines():
  r=json.loads(line)
  if 'mutation' in r:mutations.append(r);continue
  e=list(map(d,r['e_u']));n=list(map(d,r['n_u']));ge,gn=gamma(e),gamma(n)
  dot=sum(x*y for x,y in zip(e,n));gre=ge*gn-dot/(c*c)
  erest=(gre-1)*me*c*c/qe
  # Independent direct four-vector inverse boost using exact binary64 inputs.
  a=dot/(c*c*(gn+1))-ge
  rest=[e[i]+a*n[i] for i in range(3)]
  expected=norm(rest);direction=[x/expected for x in rest]
  z={'index':r['index'],'input_eV':r['input_eV'],'boost_u_over_c':r['boost_u_over_c'],'shape':r['shape'],'accepted':r['accepted'],'oracle_target_rest_eV':f(erest),'exception':r.get('exception')}
  if 'direction' in r:
   direct=list(map(d,r['direction']));diff=norm([direct[i]-direction[i] for i in range(3)])
   z['direction_vector_error']=f(diff)
   z['direction_norm_error']=f(abs(norm(direct)-1))
  if r['accepted']:
   beforeK=me*sum(x*x for x in e)/(ge+1)+mn*sum(x*x for x in n)/(gn+1)
   afterK=d(0);afterP=[d(0)]*3;afterRest=d(0);ref=d(0)
   for v in r['products']:
    mass=d(v['mass']);u=list(map(d,v['u']));afterK+=mass*sum(x*x for x in u)/(gamma(u)+1);afterRest+=mass*c*c;ref+=d(v['reference_J'])
    for i in range(3):afterP[i]+=mass*u[i]
   weight=d(.25)
   energy=(afterK+ref-beforeK)*weight
   dp=[(afterP[i]-me*e[i]-mn*n[i])*weight for i in range(3)]
   beforeRest=(me+mn)*c*c
   z['oracle_energy_defect_J']=f(energy);z['oracle_energy_defect_eV']=f(energy/weight/qe)
   z['oracle_energy_bound_fraction']=f(abs(energy)/d(r['energy_bound_J']))
   z['oracle_momentum_bound_fraction']=f(max(map(abs,dp))/d(r['momentum_bound']))
   z['invariant_mass_reference_rounding_eV']=f((afterRest-beforeRest-ref)/qe)
   z['declared_energy_defect_J']=r['energy_defect_J'];z['energy_bound_eV_per_reaction']=f(d(r['energy_bound_J'])/weight/qe)
  cases.append(z)
accepted=[r for r in cases if r['accepted']]
directions=[r for r in cases if 'direction_vector_error' in r]
summary={'cases':len(cases),'accepted':len(accepted),'rejected':len(cases)-len(accepted),'max_direction_error_accepted':max(accepted,key=lambda r:r['direction_vector_error']),'max_direction_error_returned':max(directions,key=lambda r:r['direction_vector_error']),'oracle_energy_outside_reported_bound':[r for r in accepted if r['oracle_energy_bound_fraction']>1],'oracle_momentum_outside_reported_bound':[r for r in accepted if r['oracle_momentum_bound_fraction']>1],'max_energy_bound_eV_per_reaction':max(accepted,key=lambda r:r['energy_bound_eV_per_reaction']),'mutation_cases':mutations,'oracle':'120-digit Decimal direct Lorentz four-vectors from exact binary64 input values; independent inverse boost and K+reference/product momentum closure'}
(p/'replay-oracle.json').write_text(json.dumps({'summary':summary,'cases':cases},indent=2)+'\n')
print(json.dumps(summary,indent=2))
