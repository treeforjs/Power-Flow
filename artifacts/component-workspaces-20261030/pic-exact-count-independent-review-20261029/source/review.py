from pathlib import Path
from fractions import Fraction as F
import math,json,hashlib,os,fcntl,subprocess,sys,time,tarfile
P=Path(__file__).resolve().parent.parent;S=P/'source';I=P/'sealed-input'
O=Path('<workspace>/outputs/pic-profile/exact-count-independent-review-20261029')
B=P/'build';B.mkdir(exist_ok=True)
FROZEN='/home/treefor/projects/runtimes/arch-base-20261028/run-frozen'
TINY=math.ulp(0.);MAX=float.fromhex('0x1.fffffffffffffp+1023');LIMIT=1<<2098
def val(x):return {'max':MAX,'tiny':TINY,'inf':math.inf,'-inf':-math.inf,'nan':math.nan,'-zero':-0.}.get(x)if isinstance(x,str)else float(x)
def ff(x):return F.from_float(val(x))
def units(x):return int(ff(x)/F(2)**-1074)
def greedy(n):
 out=[]
 while n:
  shift=max(0,n.bit_length()-53);top=n>>shift
  out.append(math.ldexp(float(top),shift-1074));n-=top<<shift
 return out
def op(kind,amount):return {'kind':kind,'amount':amount}
def inputs():
 rows=[]
 def leaves(name,w,ops=[]):rows.append({'id':name,'kind':'leaves','weights':w,'operations':ops})
 def two(name,hi,lo=0,ops=[]):rows.append({'id':name,'kind':'two','hi':hi,'lo':lo,'operations':ops})
 two('zero',0,0,[op('credit','tiny'),op('debit','tiny'),op('debit','tiny'),op('credit',0),op('debit',0)])
 two('minsubnormal','tiny',0,[op('debit','tiny'),op('credit','tiny'),op('credit','tiny'),op('debit',3*TINY)])
 two('negative-zero','-zero','-zero',[op('credit','-zero'),op('debit','-zero')])
 for x in [float.fromhex('0x1.0p-1022'),math.ldexp(1.,-500),1.,math.ldexp(1.,500)]:
  for lo in [0.,math.ulp(x)/2,-math.ulp(x)/2]:two('Count2/'+x.hex()+'/'+lo.hex(),x,lo,[op('debit','tiny'),op('credit','tiny')])
 two('Count2-max','max',0,[op('credit',math.ldexp(1.,970)),op('debit','tiny'),op('credit','tiny')])
 two('Count2-max-positive-small','max',math.ldexp(1.,969))
 two('Count2-max-negative-cancellation','max',-MAX)
 two('Count2-third-residual',1,-1e-20,[op('debit',1e-100),op('credit',1e-100),op('debit',1.)])
 two('Count2-source-real-scale',1e20,0,[op('debit',9.7827104392643903e-59),op('debit',math.nextafter(9.7827104392643903e-59,math.inf)),op('credit',9.7827104392643903e-59)])
 for name,hi,lo in [('negative',-1,0),('zero-low',0,TINY),('too-large-low',1,math.ulp(1.)),('nan','nan',0),('infinite','inf',0),('infinite-low',1,'inf'),('negative-infinite','-inf',0),('max-sum-overflow','max','max')]:two('invalid/'+name,hi,lo)
 for name,w in [('empty',[]),('odd-subnormals',[TINY,2*TINY,4*TINY]),('normal-transition',[math.ldexp(1.,-1022),-0.] ),('two-minnormal-plus-tiny',[math.ldexp(1.,-1021),TINY]),('RNE-even-low',[1.,2**-53]),('RNE-even-high',[math.nextafter(1.,math.inf),2**-53]),('near-one-below',[1-2**-53,2**-54]),('separated',[2**600,1.,2**-600,TINY]),('separated-reverse',[TINY,2**-600,1.,2**600]),('dense-max2098',greedy(LIMIT-1)),('dense-one-minus-min',greedy((1<<1074)-1)),('count42',[TINY]*42),('count43',[TINY]*43),('badzero',[0.]),('badnegative',[-TINY]),('badnan',['nan']),('badinf',['inf']),('credit-range-max',['max','max'])]:
  ops=[op('credit','tiny'),op('debit','tiny')]if name.startswith('dense')else[]
  leaves(name,w,ops)
 leaves('debit-credit-refusals',[1.],[op('debit','nan'),op('debit','inf'),op('debit',-1),op('debit',0),op('credit','nan'),op('credit','inf'),op('credit',-1),op('debit',math.nextafter(1.,math.inf))])
 packet=json.loads(Path('<workspace>/outputs/ss304-spatial-metal-transfer-20261029/native/native-batches.json').read_text())[0]['packets'][0]
 return {'cases':rows,'atom_mass':packet['mass_per_atom_kg'],'reference':packet['reference_per_atom_J'],'manufactured_K':1.23456789e-20,'scope':'count/material exact-currency interface; K is manufactured, no native birth or source-rate admission'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def hold():
 p=Path('/proc/278902');s=(p/'stat').read_text().split(') ',1)[1].split();assert s[0]=='T'and s[19]=='5288339'and'frozen 1'in Path('/sys/fs/cgroup/powerflow-pic-gpu-reserved/cgroup.events').read_text()
 return {'pid':278902,'state':s[0],'start_ticks':s[19],'frozen':True,'exe':os.readlink(p/'exe')}
def oracle(mode):
 data=json.loads((O/('probe-'+mode+'.json')).read_text());checks=[];totals={'cases':len(data['cases']),'admitted_imports':0,'admitted_operations':0,'refused_operations':0,'lossy_or_overflow_Count2_export_refusals':0,'largest_positive_decomposition':0}
 def check(ok,name):
  checks.append({'name':name,'passed':bool(ok)})
  if not ok:raise AssertionError(name)
 def snapshot(j,n,name):
  check(int(j['units_hex'],16)==n,name+'/exact integer')
  carriers=j['carriers'];check(carriers==greedy(n)and sum(units(x)for x in carriers)==n,name+'/canonical positive leaves')
  totals['largest_positive_decomposition']=max(totals['largest_positive_decomposition'],len(carriers))
  check(len(carriers)<=42 and all(math.isfinite(x)and x>0 for x in carriers),name+'/bounds')
  check(j['zero']==(n==0)and j['copy_units_cannot_mutate_value'],name+'/immutable copy and zero')
  quantity=F(n)*F(2)**-1074
  try:rounded=float(quantity)
  except OverflowError:rounded=math.inf
  if math.isfinite(rounded):check(j.get('rounded')==rounded,name+'/independent nearest-even rounding')
  else:check('rounded_refusal'in j,name+'/rounded overflow refusal')
  exportable=False
  if math.isfinite(rounded):
   remainder=quantity-F.from_float(rounded)
   low=float(remainder)
   exportable=F.from_float(low)==remainder
  if exportable:check('two_word'in j and sum((F.from_float(x)for x in j['two_word']),F())==quantity and j['two_word_reimport_equal'],name+'/exact Count2 export')
  else:
   check('two_word_refusal'in j,name+'/lossy Count2 export refused');totals['lossy_or_overflow_Count2_export_refusals']+=1
  def exact(x):return F(int(x['hex'],16))*F(2)**x['exponent']
  check(exact(j['material_exact_atoms'])==quantity,name+'/material port consumes all exact gas leaves')
  check(exact(j['material_mass_kg'])==quantity*F.from_float(data['atom_mass']),name+'/one material mass transfer')
  check(exact(j['material_packet_energy_J'])==quantity*(F.from_float(data['reference'])+F.from_float(data['manufactured_K'])),name+'/actual binary energy operands once')
 for row in data['cases']:
  test=row['input'];name=test['id']
  if test['kind']=='two':
   hi,lo=val(test['hi']),val(test['lo']);valid=math.isfinite(hi)and math.isfinite(lo)and math.isfinite(hi+lo)and hi>=0
   if valid:valid=(lo==0 if hi==0 else abs(lo)<=.5*(math.nextafter(hi,math.inf)-hi))
   n=units(test['hi'])+units(test['lo'])if valid else 0
  else:
   weights=list(map(val,test['weights']));valid=len(weights)<=42 and all(math.isfinite(x)and x>0 for x in weights);n=sum(units(x)for x in weights)if valid else 0
  valid=valid and 0<=n<LIMIT
  check(row['admitted']==valid,name+'/import admission')
  if not valid:continue
  totals['admitted_imports']+=1;snapshot(row['initial'],n,name+'/initial')
  for i,r in enumerate(row['operations']):
   op=r['op'];v=val(op['amount']);kind=op['kind'];allowed=math.isfinite(v)and(v>0 if kind=='debit'else v>=0)
   nextn=n
   if allowed:nextn=n+(-units(v)if kind=='debit'else units(v));allowed=0<=nextn<LIMIT
   check(r['admitted']==allowed,f'{name}/{i}/operation admission');snapshot(r['before'],n,f'{name}/{i}/before')
   if allowed:totals['admitted_operations']+=1;n=nextn;snapshot(r['after'],n,f'{name}/{i}/after')
   else:totals['refused_operations']+=1;check(r['unchanged'],f'{name}/{i}/refusal is unchanged')
 notes=[{'kind':'range-contract-clarification','text':'2098-bit exact domain includes values above DBL_MAX. Positive carriers remain exact; rounded diagnostics/Count2 export may explicitly refuse overflow. This agrees with implementation and existing range contract; callers must not require a scalar to conserve count.'},{'kind':'inherited-Count2-normalization-edge','text':'At hi=DBL_MAX, nextafter(hi,+inf)-hi is infinity, so the shared normalization predicate accepts large negative low limbs, including complete cancellation. Exact value is still conserved and export renormalizes it; this is inherited Count2 admission behavior, not silent count corruption or a new physical source admission. Tightening it needs a shared Count2 contract decision.'},{'kind':'source-material-seam','text':'Gas count uses exact positive leaves; source/adsorbed Count2 remains its own owner. Full-time source preflight must prove exact source inventory delta equals actual emitted native weights. Do not use successful gas debit/export as proof that a separate heat Count2 debit was exact. The production spatial material port can form exact atoms/mass/E from all positive native leaves without Count2 export.'}]
 result={'status':'passed_no_count_arithmetic_or_lossy_export_defect_found_in_bounded_review','mode':mode,'assertions':len(checks),'totals':totals,'checks':checks,'contract_notes':notes,'not_covered':['parent-subflow integration suite owned by root','new full-time native source/motion scope','physical packet weight validity, source-rate selection or C19 physics','all allocation-failure paths and all floating environment modes']}
 (O/('oracle-'+mode+'.json')).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'mode':mode,'assertions':len(checks),'totals':totals}))
os.sched_setaffinity(0,{7});os.nice(max(0,15-os.getpriority(os.PRIO_PROCESS,0)))
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
os.environ['CUDA_VISIBLE_DEVICES']='';mode=sys.argv[1]
with (O.parent/'testing-core7.lock').open('a')as lease:
 fcntl.flock(lease,fcntl.LOCK_EX);before=hold();start=time.time()
 H=json.loads((O/'reviewed-handoff.json').read_text());sourcefiles={str(p.relative_to(I)):sha(p)for folder in ['source','inputs']for p in (I/folder).rglob('*')if p.is_file()};identity=hashlib.sha256(json.dumps(sourcefiles,sort_keys=True).encode()).hexdigest();assert identity==H['source_identity']=='b8a9a27e6a0fcd8c54e49dd341503fed165da460180b608aa0ae024b0436fb0e'
 (O/'inputs.json').write_text(json.dumps(inputs(),indent=2)+'\n')
 flags=['-O1','-fsanitize=undefined','-fno-sanitize-recover=all']if mode=='ubsan'else ['-O2']
 cmd=[FROZEN,'/usr/bin/g++','-std=c++20','-Wall','-Wextra','-Wpedantic','-Werror',*flags,'-I'+str(I/'source/include'),'-I'+str(S),'-isystem',str(S/'vendor'),str(S/'probe.cpp'),'-lcrypto','-o',str(B/('probe-'+mode))]
 with (O/('compile-'+mode+'.log')).open('w')as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=50)
 run=[FROZEN,str(B/('probe-'+mode)),str(O/'inputs.json'),str(O/('probe-'+mode+'.json'))]
 with (O/('run-'+mode+'.log')).open('w')as log:subprocess.run(run,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=10)
 oracle(mode);after=hold();receipt={'mode':mode,'compile':cmd,'run':run,'CPU':[7],'nice':15,'threads':1,'GPU':False,'lease_s':time.time()-start,'production_before':before,'production_after':after,'source_identity':identity,'review_source_files':{str(p.relative_to(S)):sha(p)for p in S.rglob('*')if p.is_file()},'binary_sha256':sha(B/('probe-'+mode))};(O/('receipt-'+mode+'.json')).write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'mode':mode,'lease_s':receipt['lease_s']}))
