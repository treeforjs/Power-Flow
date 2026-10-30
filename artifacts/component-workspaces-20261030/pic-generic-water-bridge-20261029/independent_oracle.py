from pathlib import Path
import json,decimal,os
from fractions import Fraction
assert os.sched_getaffinity(0)=={7} and os.getpriority(os.PRIO_PROCESS,0)==15
o=Path('/home/treefor/sim-workspace/powerflow/outputs/pic-generic-water-bridge-20261029')
decimal.getcontext().prec=140
D=decimal.Decimal.from_float;C=D(299792458.)
def kinetic(m,u):
 m=D(m);v=sum((D(x)*D(x) for x in u),decimal.Decimal(0))
 return m*v/((1+v/(C*C)).sqrt()+1)
rows=json.loads((o/'outcome-records.json').read_text());checks=0;worstE=decimal.Decimal(0);worstP=decimal.Decimal(0);models=set();levels=set()
for r in rows:
 W=Fraction(r['weight']);h=oh=q=Fraction(0);ea=decimal.Decimal(0);pa=[decimal.Decimal(0)]*3;Ne=0
 for p in r['products']:
  wp=Fraction(p['weight']);assert wp==W
  h+=p['H']*wp;oh+=p['O']*wp;q+=p['q']*wp
  ea+=kinetic(p['mass'],p['u'])+D(p['internal'])+D(p['reference'])
  for d in range(3):pa[d]+=D(p['mass'])*D(p['u'][d])
  Ne+=p['state']=='electron'
  checks+=1
 assert h==2*W and oh==W and q==-W
 assert Ne==(0 if r['capture'] else 1 if r['model'] in ['A1B1_fixed_loss_dipole_ground_v1','IAM_NIST64_piecewise_mu_elastic_v1'] else 2)
 eb=kinetic(r['incoming_e_mass'],r['incoming_e_u'])+kinetic(r['incoming_n_mass'],r['incoming_n_u'])
 dw=D(r['weight']);ed=abs(ea-eb)*dw;bound=D(r['energy_bound']);assert ed<=bound,(r['model'],str(ed/bound))
 worstE=max(worstE,ed/bound)
 for d in range(3):
  before=D(r['incoming_e_mass'])*D(r['incoming_e_u'][d])+D(r['incoming_n_mass'])*D(r['incoming_n_u'][d]);pd=abs(pa[d]-before)*dw;pb=D(r['momentum_bound']);assert pd<=pb,(r['model'],str(pd/pb));worstP=max(worstP,pd/pb)
 checks+=6;models.add(r['model'])
 if r['OH_level']>=0:levels.add(r['OH_level'])
assert len(models)==8 and levels==set(range(5))
result={'status':'passed','method':'Fraction exact binary64 weights/nuclei/charge; 140-decimal-digit independent proper-velocity kinetic energy and momentum','records':len(rows),'assertions':checks,'models':len(models),'OH_levels':sorted(levels),'worst_energy_defect_over_unchanged_bound':str(worstE),'worst_momentum_defect_over_unchanged_bound':str(worstP),'physical_evolution':False}
(o/'independent-oracle.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
