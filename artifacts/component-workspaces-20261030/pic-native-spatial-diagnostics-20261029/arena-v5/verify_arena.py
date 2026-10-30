from pathlib import Path
import datetime,fcntl,hashlib,json,os,subprocess,struct,math
from fractions import Fraction
from decimal import Decimal,localcontext
P=Path('/home/treefor/projects/runtimes/pic-reaction-transactions-20261029/native-audit-fix');R=Path('<workspace>/outputs/pic-profile/reaction-transactions-20261029/native-audit-fix')
assert os.sched_getaffinity(0)=={7}and os.getpriority(os.PRIO_PROCESS,0)==15
with(R.parent.parent/'testing-core7.lock').open('a')as h:
 try:fcntl.flock(h,fcntl.LOCK_EX|fcntl.LOCK_NB)
 except BlockingIOError:pass
 else:raise RuntimeError('CPU7 lease required')
C=json.loads((R/'commands.json').read_text());normal=json.loads((R/'tests/arena-tests.json').read_text());frozen=C['run'][0];B=P/'build';commands=[];objects=[str(B/'native.o')]
for old in C['build'][1:3]:
 command=[('-O1'if x=='-O2'else x)for x in old];command[-1]=str(B/('sanitized-'+Path(command[-1]).name));command[command.index('-c'):command.index('-c')]=['-fsanitize=address,undefined','-fno-omit-frame-pointer'];commands.append(command)
 with(R/('sanitizer-'+Path(command[-1]).stem+'.log')).open('w')as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=300)
 objects.append(command[-1])
link=list(C['build'][-1]);i=link.index(str(B/'native.o'));link=link[:i]+objects+link[i+3:];link[-1]=str(B/'test_arena_sanitized');link[2:2]=['-fsanitize=address,undefined'];commands.append(link)
with(R/'sanitizer-link.log').open('w')as log:subprocess.run(link,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=90)
out=R/'sanitized-tests';out.mkdir(exist_ok=True);env=dict(os.environ,ASAN_OPTIONS='detect_leaks=1:abort_on_error=1',UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1')
with(R/'sanitizer-tests.log').open('w')as log:subprocess.run([frozen,link[-1],C['run'][2],str(out)],stdout=log,stderr=subprocess.STDOUT,env=env,check=True,timeout=120)
assert json.loads((out/'arena-tests.json').read_text())==normal
# Foundation suite and actual ABI probes against the modified native observer.
compile=list(C['build'][2]);compile[compile.index(str(P/'tests/test_arena.cpp'))]=str(P/'source/tests/test_mfem_gap_pic.cpp');compile[-1]=str(B/'foundation.o');commands.append(compile)
with(R/'foundation-build.log').open('w')as log:subprocess.run(compile,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=180)
link=list(C['build'][-1]);i=link.index(str(B/'native.o'));link=link[:i]+[str(B/'native.o'),str(B/'foundation.o')]+link[i+3:];link[-1]=str(B/'test_foundation');commands.append(link)
with(R/'foundation-link.log').open('w')as log:subprocess.run(link,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=90)
with(R/'foundation-tests.log').open('w')as log:subprocess.run([frozen,link[-1]],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=90)
assert 'foundation tests passed'in(R/'foundation-tests.log').read_text()
abi=json.loads((R.parent/'native-integration/abi-comparison.json').read_text())['sizes'];assert normal['ABI']==[80,168,200,abi['hall_baseline'][0]]
# Independent exact binary64 fraction checks of actual embedded native packets.
def F(x):return Fraction.from_float(float(x))
def count(n):return F(n[0])+F(n[1])
checks=[];energies=[]
for path in sorted((R/'tests').glob('case-*.json'))+list((R/'tests').glob('neutral.json')):
 j=json.loads(path.read_text())['payload'];raw=bytes.fromhex(j['native_bytes_hex']);assert hashlib.sha256(raw).hexdigest()==j['native_bytes_sha256'];n=struct.unpack_from('<Q',raw,384)[0];assert len(raw)==392+80*n
 native={}
 for k in range(n):
  v=struct.unpack_from('<Q9d',raw,392+80*k);assert v[0]not in native;native[v[0]]=v
 specs={s['id']:s for s in j['registry']['states']};owned=set()
 for record in j['metadata']['packets']:
  rows=[native[i]for i in record['native_ids']];assert not(owned&set(record['native_ids']));owned.update(record['native_ids']);assert sum((F(v[9])for v in rows),Fraction())==count(record['count'])
  spec=specs[record['state']]
  for v in rows:assert v[1:7]==rows[0][1:7]and v[8]==spec['mass_kg']and v[7]==float(spec['q_e'])*1.602176634e-19 and v[9]>0
 assert owned==set(native)
 for receipt in j['metadata']['sinks']:assert sum((F(v['weight'])for v in receipt['native_carrier_receipts']),Fraction())==count(receipt['count'])
 checks.append({'bundle':path.name,'logical_packets':len(j['metadata']['packets']),'positive_native_carriers':n,'exact_fraction_count_checks':'passed'})
#120-digit local event ledgers from actual recorded operands, no huge reservoir masks.
with localcontext()as ctx:
 ctx.prec=120
 D=lambda x:Decimal.from_float(float(x));c=D(299792458.);qe=D(1.602176634e-19)
 def K(m,u):
  u2=sum((D(x)*D(x)for x in u),Decimal());return D(m)*u2/((1+u2/(c*c)).sqrt()+1)
 specs={s['id']:s for s in json.loads((R/'tests/case-1.json').read_text())['payload']['registry']['states']}
 for proof in normal['proofs']:
  if proof['kind']!='reaction':continue
  assert 'primary_u'in proof and 'target_u'in proof
  e=specs['electron'];n=specs['H2O'];before=K(e['mass_kg'],proof['primary_u'])+D(e['reference_J'])+D(e['internal_J'])+K(n['mass_kg'],proof['target_u'])+D(n['reference_J'])+D(n['internal_J']);after=Decimal();dp=[-D(e['mass_kg'])*D(proof['primary_u'][d])-D(n['mass_kg'])*D(proof['target_u'][d])for d in range(3)];scale=sum((abs(v)for v in dp),Decimal());q=-1;h=-2;o=-1
  for product in proof['products']:
   s=specs[product['state']];after+=K(s['mass_kg'],product['u'])+D(s['reference_J'])+D(s['internal_J']);q-=s['q_e'];h+=s['H'];o+=s['O']
   for d in range(3):v=D(s['mass_kg'])*D(product['u'][d]);dp[d]+=v;scale+=abs(v)
  assert(q,h,o)==(0,0,0);defect=(after-before)*D(proof['event_weight']);assert abs(defect)<=D(proof['energy_roundoff_bound_J']);limit=D(256)*D(2**-52)*scale*D(proof['event_weight']);assert all(abs(v*D(proof['event_weight']))<=limit for v in dp)
  assert count(proof['primary_before'])-F(proof['event_weight'])==count(proof['primary_remainder']);assert count(proof['target_before'])-F(proof['event_weight'])==count(proof['target_remainder']);energies.append({'event':proof['event_id'],'physical_stage_s':proof['physical_stage_s'],'energy_defect_J_120_digit':str(defect),'energy_bound_J':proof['energy_roundoff_bound_J'],'H_O_charge':'exact','momentum':'passed'})
(R/'independent-ledgers.json').write_text(json.dumps({'status':'passed','fraction_native_count_checks':checks,'decimal_precision':120,'local_event_ledgers':energies,'limits':'Exact represented count checks; 120-digit ledger replay of actual binary64 operands. Not differential-model validation, physical rate selection, MCC heating or Maxwell-current qualification.'},indent=2)+'\n')
(R/'verification-commands.json').write_text(json.dumps(commands,indent=2)+'\n')
(R/'verification-result.json').write_text(json.dumps({'status':'passed','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'checks':normal['check_count'],'sanitizers':'address/undefined/leak clean on arena/tests; frozen native/MFEM uninstrumented','normal_sanitized_results_identical':True,'unchanged_native_foundation_suite':'passed','ABI':[80,168,200,abi['hall_baseline'][0]],'independent_fraction_and_120digit_ledgers':'passed','GPU':False},indent=2)+'\n')
print('Canonical arena verification passed')
