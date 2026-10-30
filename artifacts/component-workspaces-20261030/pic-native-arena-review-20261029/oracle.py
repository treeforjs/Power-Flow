from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal as D,localcontext
import math,json,os
assert sorted(os.sched_getaffinity(0))==[7]
p=Path('/home/treefor/sim-workspace/powerflow/outputs/pic-native-arena-review-20261029')
cases={r['id']:r for r in json.loads((p/'count-cases.json').read_text())};rows=json.loads((p/'count-results.json').read_text())
summary={'count_cases':len(rows),'carrier_admissions':0,'exact_debit_admissions':0,'overdraw_admissions':[],'carrier_sum_errors':[],'nonpositive_or_nonfinite_carriers':[],'debit_errors':[],'maximum_carriers':0,'near_cancel_maximum_cases_admitted':[]}
for row in rows:
 c=cases[row['id']];before=F(c['hi'])+F(c['lo'])
 if row['carrier_accepted']:
  summary['carrier_admissions']+=1;weights=row['carriers'];summary['maximum_carriers']=max(summary['maximum_carriers'],len(weights))
  if sum(map(F,weights),F())!=before:summary['carrier_sum_errors'].append(row['id'])
  if not all(math.isfinite(x) and x>0 for x in weights):summary['nonpositive_or_nonfinite_carriers'].append(row['id'])
  if c['hi']>1e308 and c['lo']<-1e300:summary['near_cancel_maximum_cases_admitted'].append(row['id'])
 if row.get('debit_accepted'):
  summary['exact_debit_admissions']+=1
  if F(c['debit'])>before:summary['overdraw_admissions'].append(row['id'])
  if sum(map(F,row['debit_result']),F())!=before-F(c['debit']):summary['debit_errors'].append(row['id'])
assert not any(summary[k] for k in ('overdraw_admissions','carrier_sum_errors','nonpositive_or_nonfinite_carriers','debit_errors'))
probes=json.loads((p/'probes.json').read_text());byname={x['name']:x for x in probes}
control=byname['actual_parent_publish_rollback'];assert control['accepted']
counts={}
for row in control['packets']:
 n=sum(map(F,row['count']),F());assert sum((F(v['weight']) for v in row['carriers']),F())==n
 counts[row['state']]=counts.get(row['state'],F())+n
ew=F(.031334406939830223);nw=F(64.)
charge=-counts['electron']+counts['H2O+'];H=2*counts['H2O']+2*counts['H2O+'];O=counts['H2O']+counts['H2O+']
assert charge==-ew and H==2*nw and O==nw
with localcontext() as ctx:
 ctx.prec=120
 d=lambda x:D.from_float(x)
 c=d(299792458.);qe=d(1.602176634e-19);mass_e=d(9.1093837139e-31)
 initial_u=list(map(d,control['audit']['primary_u']))
 def K(m,u):return m*sum(x*x for x in u)/((1+sum(x*x for x in u)/(c*c)).sqrt()+1)
 before=K(mass_e,initial_u)*d(float(ew));after=D(0);momentum=[D(0)]*3
 for row in control['packets']:
  for v in row['carriers']:
   mass=d(v['mass_kg']);u=list(map(d,v['u']));w=d(v['weight']);after+=K(mass,u)*w
   if row['state']=='H2O+':after+=d(12.621)*qe*w
   for i in range(3):momentum[i]+=mass*u[i]*w
 dp=[momentum[i]-mass_e*initial_u[i]*d(float(ew)) for i in range(3)]
 energy={'weighted_energy_defect_J':float(after-before),'weighted_energy_receipt_bound_J':control['audit']['energy_roundoff_bound_J'],'weighted_momentum_defect_kg_m_s':[float(x) for x in dp],'charge_e_exact':str(charge),'H_count_exact':str(H),'O_count_exact':str(O),'full_neutral_and_electron_inventory_carrier_sums_exact':True,'oracle_decimal_digits':120}
 assert abs(after-before)<=d(control['audit']['energy_roundoff_bound_J'])
(p/'independent-oracle.json').write_text(json.dumps({'count_summary':summary,'actual_parent_control':energy,'probe_acceptance':{r['name']:r.get('accepted') for r in probes}},indent=2)+'\n')
print(json.dumps(summary));print(json.dumps(energy));print(json.dumps({r['name']:r.get('accepted') for r in probes}))
