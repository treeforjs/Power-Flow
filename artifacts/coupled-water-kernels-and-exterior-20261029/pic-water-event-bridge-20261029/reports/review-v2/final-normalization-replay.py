from pathlib import Path
import os,json,hashlib,shutil,subprocess
p=Path(__file__).parent;r=Path('/home/treefor/projects/runtimes/pic-water-event-bridge-20261029')
assert sorted(os.sched_getaffinity(0))==[7]
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
assert sha(r/'source/water_parent_bridge.hpp')=='5a3ef505ab8bce21b5f6c0a05e1b9718b2e180535c0b646b2a5ba1a1fe577473'
old=p/'source-before-normalization.hpp'
shutil.copy2(p/'source/water_parent_bridge.hpp',old)
shutil.copy2(r/'source/water_parent_bridge.hpp',p/'source/water_parent_bridge.hpp')
shutil.copy2(r/'README.txt',p/'owner-README.txt')
shutil.copy2(p.parent/'checks.json',p/'owner-checks.json')
bindings=json.loads((p/'input-bindings.json').read_text())
for rec in bindings:
 f=Path(rec['snapshot'])
 if f==p/'source/water_parent_bridge.hpp':rec['sha256']=sha(f);rec['bytes']=f.stat().st_size
 if f==p/'owner-checks.json':rec['original']=str(p.parent/'checks.json');rec['sha256']=sha(f);rec['bytes']=f.stat().st_size
for f,original in ((old,'historical reviewed v2 before normalized-direction bound'),(p/'owner-README.txt',str(r/'README.txt'))):
 bindings.append({'snapshot':str(f),'original':original,'sha256':sha(f),'bytes':f.stat().st_size})
(p/'input-bindings.json').write_text(json.dumps(bindings,indent=2)+'\n')
text=(p/'review.txt').read_text().replace('b014aa73d1b43c894da6757b35c4395a770b3bb87186ea39dafca6c6deadc687',sha(r/'source/water_parent_bridge.hpp'))
oldtext='''must be <=1e-12. The original magnitude check remains. The1e-12 criterion is an
unnormalized-vector conditioning admission. If advertising a certified normalized
direction error from a generic vector bound b, use the conservative normalization
bound2b/(|rest|-b), including final normalization/cast error. The actual admitted
oracle errors here are far below either threshold.'''
newtext='''must satisfy b<|rest| and2b/(|rest|-b)<=1e-12 in the final reviewed source,
accounting conservatively for vector normalization. The original magnitude check
remains. The final binary64 direction cast has an additional ordinary rounding
error; the actual admitted oracle errors here are far below the declared bound.
The preceding b014 header is preserved in source-before-normalization.hpp.'''
assert oldtext in text
text=text.replace(oldtext,newtext)
(p/'review.txt').write_text(text)
subprocess.run(['python3',str(p/'prepare-and-run.py')],check=True)
subprocess.run(['python3',str(p/'finalize.py')],check=True)
