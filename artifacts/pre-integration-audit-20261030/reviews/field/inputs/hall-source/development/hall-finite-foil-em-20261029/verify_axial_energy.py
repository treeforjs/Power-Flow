from pathlib import Path
from decimal import Decimal,localcontext
import json,sys
run=Path(sys.argv[1]);before=json.loads((run/'split/state.json').read_text());after=json.loads((run/'full/state.json').read_text())
mesh=json.loads((run/'full/manifest.json').read_text())['mesh'];nx=mesh['nx'];ny=mesh['ny'];nz=mesh['nz'];dx=mesh['dx_m'];dy=mesh['dy_m'];dz=mesh['dz_m']
eps=8.8541878188e-12;speed=299792458.;mu=1/(eps*speed*speed)
dual=lambda i,n,h:.5*h if i in (0,n) else h
M=[];H=[]
for k in range(nz+1):
 for j in range(ny+1):
  for i in range(nx):M.append(eps*(dx*dual(j,ny,dy)*dual(k,nz,dz)))
for k in range(nz+1):
 for j in range(ny):
  for i in range(nx+1):M.append(eps*(dual(i,nx,dx)*dy*dual(k,nz,dz)))
for k in range(nz):
 for j in range(ny+1):
  for i in range(nx+1):M.append(eps*(dual(i,nx,dx)*dual(j,ny,dy)*dz))
for k in range(nz):
 for j in range(ny):
  for i in range(nx+1):H.append((dual(i,nx,dx)*dy*dz)/mu)
for k in range(nz):
 for j in range(ny+1):
  for i in range(nx):H.append((dx*dual(j,ny,dy)*dz)/mu)
for k in range(nz+1):
 for j in range(ny):
  for i in range(nx):H.append((dx*dy*dual(k,nz,dz))/mu)
D=Decimal.from_float
with localcontext() as ctx:
 ctx.prec=80
 def energy(state):
  total=Decimal(0)
  for key,weight in [('E',M),('B',H)]:
   assert len(state[key])==len(state[key+'_low'])==len(weight)
   for high,low,w in zip(state[key],state[key+'_low'],weight):
    x=D(high)+D(low);total+=D(w)*x*x/2
  return total
 old=energy(before);new=energy(after);change=new-old
 recorded=D(after['work'][0])-D(before['work'][0]);defect=change-recorded
 limit=Decimal('1e-20')+Decimal('2e-13')*abs(recorded)
 result={'status':'passed' if abs(defect)<=limit else 'failed','method':'independent 80-digit sum of physical quadratic energy from high+low checkpoint fields; SI Hodge weights reconstructed from mesh',
  'before_revision':before['revision'],'after_revision':after['revision'],'exact_decimal_energy_change_J':str(change),'recorded_energy_change_J':str(recorded),
  'defect_J':float(defect),'unchanged_absolute_limit_J':1e-20,'allowed_J':float(limit),'corrective_energy_added_J':0}
 (run/'independent-energy.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result));assert result['status']=='passed'
