"""Independent high-precision check of manufactured receipts, not a rate model."""
from pathlib import Path
from decimal import Decimal,localcontext
import json,math,os,fcntl
ROOT=Path('<workspace>/outputs/pic-profile')
OUT=ROOT/'reaction-transactions-20261029'
assert os.sched_getaffinity(0)=={7}and os.getpriority(os.PRIO_PROCESS,0)==15
with(ROOT/'testing-core7.lock').open('a')as h:
 try:fcntl.flock(h,fcntl.LOCK_EX|fcntl.LOCK_NB)
 except BlockingIOError:pass
 else:raise RuntimeError('CPU7 lease required')
catalog=json.loads((OUT/'tests/manufactured-catalog.json').read_text())
checkpoint=json.loads((OUT/'tests/manufactured-continuation.restart.json').read_text())['payload']
species={s['id']:s for s in catalog['species']};states={s['id']:s for s in catalog['states']}
D=lambda x:Decimal(x)if isinstance(x,int)else Decimal.from_float(float(x))
eps=D(float.fromhex('0x1p-52'));factor=D(256)
with localcontext()as context:
 context.prec=80;c2=D(299792458)**2
 def moments(records):
  result={'H':D(0),'O':D(0),'charge_e':D(0),'kinetic_J':D(0),'internal_J':D(0),'reference_J':D(0),'px':D(0),'py':D(0),'pz':D(0)}
  abs_energy=D(0);abs_momentum=[D(0)]*3
  for state,w,u in records:
   st=states[state];sp=species[st['species']];w=D(w);u=[D(x)for x in u];mass=D(sp['mass_kg'])
   u2=sum(x*x for x in u);K=w*mass*u2/((1+u2/c2).sqrt()+1)
   internal=w*D(st['internal_J']);reference=w*D(st['reference_J'])
   result['kinetic_J']+=K;result['internal_J']+=internal;result['reference_J']+=reference;abs_energy+=abs(K)+abs(internal)+abs(reference)
   result['H']+=w*sp['H'];result['O']+=w*sp['O'];result['charge_e']+=w*sp['charge_e']
   for k,key in enumerate(['px','py','pz']):
    value=w*mass*u[k];result[key]+=value;abs_momentum[k]+=abs(value)
  result['total_J']=result['kinetic_J']+result['internal_J']+result['reference_J']
  return result,abs_energy,abs_momentum
 packets={p['id']:p.copy()for p in checkpoint['genesis']['packets']};parents={p['id']:p.copy()for p in checkpoint['genesis']['parents']}
 def global_records():
  return[(p['state'],p['weight'],p['proper_velocity_m_s'])for p in packets.values()]+[(p['state'],p['remaining_weight'],p['proper_velocity_m_s'])for p in parents.values()]
 initial,_,_=moments(global_records());evidence=[]
 for record in checkpoint['journal']:
  event=record['event'];primary=packets.pop(event['primary_id']);parent=parents[event['parent_id']];w=event['weight']
  before,scale0,p0=moments([(primary['state'],w,primary['proper_velocity_m_s']),(parent['state'],w,parent['proper_velocity_m_s'])])
  after,scale1,p1=moments([(p['state'],p['weight'],p['proper_velocity_m_s'])for p in event['products']])
  for key in ['H','O','charge_e']:assert before[key]==after[key]
  for k,key in enumerate(['px','py','pz']):assert abs(after[key]-before[key])<=factor*eps*max(p0[k]+p1[k],D(1e-300))
  residual=after['total_J']-before['total_J'];limit=factor*eps*(scale0+scale1);assert abs(residual)<=limit
  parent['remaining_weight']-=w
  for p in event['products']:
   assert p['weight']==w;packets[p['id']]={**p,'position_m':parent['position_m']}
  evidence.append({'event_id':event['id'],'energy_residual_J_80digit':str(residual),'energy_limit_J':str(limit),'energy_relative_residual':str(residual/(scale0+scale1)),'momentum_residuals_80digit':[str(after[k]-before[k])for k in ['px','py','pz']],'parent_debit_count':1})
 final,scale,_=moments(global_records());global_energy_residual=final['total_J']-initial['total_J']
 assert all(final[k]==initial[k]for k in ['H','O','charge_e'])
 assert abs(global_energy_residual)<=factor*eps*scale
 assert parents['gas-parent-lot/1']['remaining_weight']==.25
 assert set(packets)=={p['id']for p in checkpoint['snapshot']['packets']}
 result={'status':'passed','precision_decimal_digits':80,'independent_of_cpp_receipt_values':True,'events':evidence,'global_H':str(final['H']),'global_O':str(final['O']),'global_charge_e':str(final['charge_e']),'global_energy_residual_J':str(global_energy_residual),'scope':'Manufactured closed reaction events only, without field/material/circuit/global-PIC energy claims.'}
(OUT/'independent-ledger-check.json').write_text(json.dumps(result,indent=2)+'\n')
print('Independent80-digit atom/charge/momentum/energy checks passed for',len(evidence),'events')
