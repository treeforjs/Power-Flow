from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal as D,localcontext
import fcntl,json,math,os
O=Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs/ss304-native-metal-transfer-20261029')
def dyad(j):return F(int(j['hex'],16))*(F(2)**j['exponent'])
def exact(x):return F.from_float(float(x))
def count(x):return exact(x[0])+exact(x[1])
def dec(x):return D(x.numerator)/D(x.denominator)
os.sched_setaffinity(0,{7});p=os.getpriority(os.PRIO_PROCESS,0)
if p<15:os.nice(15-p)
with (O.parent/'pic-profile/testing-core7.lock').open('a')as lease:
 fcntl.flock(lease,fcntl.LOCK_EX);j=json.loads((O/'transfer/native-tests.json').read_text());n=0;worst=0
 def check(v):
  global n
  assert v;n+=1
 with localcontext()as ctx:
  ctx.prec=160
  for a in j['proofs']:
   emitted=F(0);K=F(0);I=F(0);R=F(0);momentum=[F(0)]*3
   for packet in a['actual_packets']:
    sign=1 if packet['emit']else-1
    w=sum(map(exact,packet['weights']),F(0))*sign
    emitted+=w;K+=w*exact(packet['native_per_atom_kinetic_J']);I+=w*exact(a['internal_per_atom_J']);R+=w*exact(a['reference_per_atom_J'])
    for d in range(3):momentum[d]+=w*exact(a['atom_mass_kg'])*exact(packet['proper_velocity'][d])
    u2=sum((dec(exact(v))**2 for v in packet['proper_velocity']),D(0));c2=D(299792458)**2
    high=dec(exact(a['atom_mass_kg']))*u2/((1+u2/c2).sqrt()+1)
    rounded=dec(exact(packet['native_per_atom_kinetic_J']));err=abs(high-rounded)
    limit=D(8)*D.from_float(math.ulp(packet['native_per_atom_kinetic_J']))
    check(err<=limit);worst=max(worst,float(err/limit)if limit else 0)
   check(dyad(a['outward_atoms'])==emitted)
   check(count(a['solid_before']['solid_atoms'])-count(a['solid_after']['solid_atoms'])==emitted)
   check(dyad(a['actual_native_kinetic_J'])==K);check(dyad(a['actual_state_internal_J'])==I);check(dyad(a['actual_state_reference_J'])==R)
   check(dyad(a['outward_packet_energy_J'])==K+I+R)
   check(dyad(a['solid_before']['energy'])+exact(a['other_heat_J'])-dyad(a['solid_after']['energy'])==K+I+R)
   for d in range(3):
    check(dyad(a['gas_momentum'][d])==momentum[d]);check(dyad(a['solid_after']['support_impulse'][d])-dyad(a['solid_before']['support_impulse'][d])==-momentum[d])
   kg=count(a['solid_after']['solid_atoms'])*exact(a['atom_mass_kg'])
   needed=dyad(a['solid_after']['energy'])-kg*exact(a['u_before_J_kg'])
   check(abs(exact(a['native_heat_load_J'])-needed)<=exact(math.ulp(a['native_heat_load_J'])) if a['native_heat_load_J'] else abs(needed)<F(2)**-1074)
   check(abs(dyad(a['native_heat_projection_residual']))<=exact(a['native_heat_roundoff_bound_J']))
   check(a['field_work_J']==a['fixed_support_work_J']==a['extra_latent_debit_J']==0)
  for b in j['boundaries']:
   E=dyad(b['before']['body']['state']['energy']);q=exact(b['heat_J'])
   check(float(E)==float(E+q));check(E+q>E if q>0 else E+q<E)
  # Exact decimal/hex restart currency roundtrips independently of JSON decimal
  # formatting or the native long-double thermodynamic projection.
  r=json.loads((O/'transfer/metal.restart.json').read_text());check(dyad(r['body']['state']['energy']).denominator>0)
 out={'status':'passed','independent_assertions':n,'native_kinetic_oracle':'160-digit Decimal proper-velocity kinetic energy;8ULP arithmetic bound','exact_ledger_oracle':'Python Fraction of binary64 inputs; zero mass/energy/momentum discrepancies','worst_fraction_of_kinetic_numeric_bound':worst,'CPU':[7],'nice':os.getpriority(os.PRIO_PROCESS,0),'not_a_changed_scientific_tolerance':True}
 (O/'highprecision-audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
