from pathlib import Path
import os,json,hashlib,subprocess,math,sys
p=Path(__file__).parent;s=p/'snapshot';o=Path('/home/treefor/sim-workspace/powerflow/outputs/pic-native-arena-review-20261029')
assert sorted(os.sched_getaffinity(0))==[7] and os.getpriority(os.PRIO_PROCESS,0)==15
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
identity=json.loads((o/'source-identity.json').read_text())
for name,want in identity['files'].items():assert sha(s/name)==want,name
for path,want in json.loads((o/'native-library-hashes.json').read_text()).items():assert sha(Path(path))==want,path
cases=[]
for exponent in (-1074,-1022,-1000,-500,-100,-1,0,1,100,500,1000,1023):
 hi=math.ldexp(1.,exponent);ulp=math.ulp(hi)
 for lo in (0.,ulp/4,-ulp/4,math.ulp(0.),-math.ulp(0.)):
  for weight in (hi,hi/2,math.ulp(0.),ulp,1e-300,0.031334406939830223):
   if weight>0 and math.isfinite(hi+lo):cases.append({'id':len(cases),'hi':hi,'lo':lo,'debit':weight})
maximum=sys.float_info.max
for lo in (0.,-math.ulp(maximum)/4,-maximum,-math.nextafter(maximum,0.),-maximum/2):
 for amount in (math.ulp(0.),math.ulp(maximum),maximum,math.nextafter(maximum,0.)):
  cases.append({'id':len(cases),'hi':maximum,'lo':lo,'debit':amount})
cases.extend([{'id':len(cases),'hi':1.,'lo':-math.nextafter(1.,0.),'debit':math.ulp(0.)},{'id':len(cases)+1,'hi':1.,'lo':-1.,'debit':math.ulp(0.)}])
(o/'count-cases.json').write_text(json.dumps(cases)+'\n')
(p/'build').mkdir(exist_ok=True)
base='/home/treefor/projects/runtimes/mfem-4.10-pic-20261029'
includes=['-I'+str(s/'source/include'),'-isystem',str(s/'inputs'),'-I'+str(s/'vendor')]
for name in ('install-mpi/include','tpl/hypre/include','tpl/metis/include'):includes+=['-isystem',base+'/'+name]
compiler=['/usr/bin/mpic++','-std=c++20','-O2','-g','-ffunction-sections','-fdata-sections']+includes+['-DPIC_REACTION_NATIVE_SOURCE_SHA256="'+identity['sha256']+'"']
commands=[];objects=[]
for source,name,strict in ((s/'source/src/mfem_gap_pic.cpp','native',False),(s/'source/src/mfem_pic_arena.cpp','arena',True),(p/'probe.cpp','probes',True)):
 command=compiler+(['-Wall','-Wextra','-Wpedantic','-Werror'] if strict else [])+['-c',str(source),'-o',str(p/'build'/f'{name}.o')];commands.append(command)
 with (o/f'build-{name}.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
 objects.append(str(p/'build'/f'{name}.o'))
libs=[base+'/install-mpi/lib/libmfem.a',base+'/tpl/hypre/lib/libHYPRE.a',base+'/tpl/metis/lib/libmetis.a',base+'/tpl/GKlib/lib/libGKlib.a','-lcrypto']
command=['/usr/bin/mpic++','-Wl,--gc-sections']+objects+libs+['-o',str(p/'build/probes')];commands.append(command)
subprocess.run(command,check=True,stdout=(o/'link.log').open('w'),stderr=subprocess.STDOUT)
command=['/usr/bin/time','-v','-o',str(o/'resources.txt'),str(p/'build/probes'),str(s/'inputs/h2o-gas-20261029.json'),str(o/'owner-tests'),str(o),str(o/'count-cases.json')];commands.append(command)
subprocess.run(command,check=True,stdout=(o/'probes.txt').open('w'),stderr=subprocess.STDOUT)
(o/'commands-independent.json').write_text(json.dumps(commands,indent=2)+'\n')
print((o/'probes.txt').read_text())
