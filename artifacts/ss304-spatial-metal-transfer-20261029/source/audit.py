from fractions import Fraction as F
from decimal import Decimal as Q,localcontext
from pathlib import Path
import json,math,sys
path=Path(sys.argv[1]);data=json.loads(path.read_text())
checks=[]
def check(ok,name):
 checks.append({'name':name,'passed':bool(ok)})
 if not ok:raise AssertionError(name)
def exact(x):
 n=int(x['hex'],16);e=x['exponent']
 return F(n)*F(2)**e
def f(x):return F.from_float(float(x))
def dec(x):return Q(x.numerator)/Q(x.denominator)
archive=json.loads((path.parent/'native-batches.json').read_text())
batches={x['binding']['gas_trial_sha256']:x for x in archive}
max_ulp=Q(0)
with localcontext()as ctx:
 ctx.prec=160
 for case,p in enumerate(data['proofs']):
  row=p['row'];before=row['before']['cells'];after=row['after']['cells'];N=len(before)
  M=[exact(x['mass_kg'])for x in before];E=[exact(x['energy_J'])for x in before]
  support=[[exact(v)for v in x['support_impulse']]for x in before]
  for h in row['heat']:E[h['cell']]+=exact(h['energy_J'])
  packetM=F(0);packetE=F(0)
  for r in row['mapping']:
   i=r['cell'];dm=exact(r['outward_mass_kg']);de=exact(r['outward_energy_J'])
   M[i]-=dm;E[i]-=de;packetM+=dm;packetE+=de
  check(sum(M)-sum(exact(x['mass_kg'])for x in before)+packetM==0,f'{case}/exact global mass')
  batch=batches.get(row['binding']['gas_trial_sha256'])
  if batch:
   mapped={x['receipt_id']:x for x in row['mapping']}
   for r in batch['packets']:
    sign=1 if r['emit']else -1
    count=sum(f(x['weight'])for x in r['carriers'])
    check(count==exact(r['exact_atoms']),f'{case}/{r["receipt_id"]}/native exact carrier count')
    mass=f(r['mass_per_atom_kg']);R=f(r['reference_per_atom_J']);u=list(map(f,r['proper_velocity_m_s']))
    energy=sum(f(c['weight'])*(f(c['kinetic_per_atom_J'])+R)for c in r['carriers'])
    entry=mapped[r['receipt_id']];i=entry['cell']
    check(exact(entry['outward_mass_kg'])==sign*count*mass,f'{case}/{r["receipt_id"]}/actual packet mass')
    check(exact(entry['outward_energy_J'])==sign*energy,f'{case}/{r["receipt_id"]}/K plus reference exactly once')
    for d in range(3):support[i][d]-=sign*count*mass*u[d]
    expected=dec(mass)*sum(dec(v)**2 for v in u)/((1+sum(dec(v)**2 for v in u)/Q(299792458)**2).sqrt()+1)
    for c in r['carriers']:
     raw=float(c['kinetic_per_atom_J']);ulps=abs(Q.from_float(raw)-expected)/Q.from_float(math.ulp(raw))
     max_ulp=max(max_ulp,ulps)
     check(ulps<=8,f'{case}/{r["receipt_id"]}/160digit native K under8ULP')
  elapsed=F(0)
  for k,sub in enumerate(row['conduction']):
   dt=exact(sub['duration_s']);elapsed+=dt;net=F(0)
   for edge in sub['fluxes']:
    flux=exact(edge['energy_J']);check(flux==dt*f(edge['conductance_W_K'])*f(edge['temperature_difference_K']),f'{case}/{k}/{edge["a"]}-{edge["b"]}/actual Fourier operand product')
    E[edge['a']]-=flux;E[edge['b']]+=flux;net-=flux;net+=flux
   check(net==0,f'{case}/{k}/pairwise exact conduction closure')
  check(elapsed==f(row['binding']['dt_s']),f'{case}/authoritative local conduction dt')
  for i,x in enumerate(after):
   check(exact(x['mass_kg'])==M[i],f'{case}/{i}/one canonical final mass')
   check(exact(x['energy_J'])==E[i],f'{case}/{i}/one canonical final energy')
   check(list(map(exact,x['support_impulse']))==support[i],f'{case}/{i}/fixed support actual momentum')
   diag=p['after_transfer_diagnostics'][i];T=Q.from_float(diag['T_K'])
   # All exported nonuniform conduction fixtures are in the sealed hot liquid.
   if T>1800:
    expected=dec(M[i])*Q.from_float(794.96)
    error=abs(Q.from_float(diag['capacity_J_K'])-expected)
    check(error<=Q(4)*Q.from_float(math.ulp(diag['capacity_J_K'])),f'{case}/{i}/capacity is current mass times sealed liquid cp')
  check(exact(p['exact_mass_residual'])==0 and exact(p['exact_energy_residual'])==0,f'{case}/exact reported residuals')
  initial=sum(exact(x['energy_J'])for x in before)
  heat=sum(exact(x['energy_J'])for x in row['heat'])
  check(sum(E)-initial+packetE-heat==0,f'{case}/independent full energy equation')
result={'status':'passed','precision_decimal_digits':160,'assertions':len(checks),'maximum_native_kinetic_ULP_difference':str(max_ulp),'checks':checks,'scope':'Exact represented-binary currency and independent highprecision arithmetic; no claim of spatial physical convergence or full-time native source integration.'}
Path(sys.argv[2]).write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k]for k in ['status','assertions','maximum_native_kinetic_ULP_difference']}))
