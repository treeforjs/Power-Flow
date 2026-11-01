#!/usr/bin/env python3
import hashlib
import json
import pathlib
import shutil

R=pathlib.Path(__file__).resolve().parent
REPO=R.parent.parent/'Power-Flow'
OUT=R.parent.parent/'.workspace/outputs'/R.name
BASE=REPO/'artifacts/gas-candidate-binding-v1.4-20261101'
ART=REPO/'artifacts/gas-attempt-rollback-v1.5-20261101'
BACK=pathlib.Path('/home/treefor/Documents/Codex/reconstructed-source-20261029')/ART.name
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def encoded(x): return (json.dumps(x,sort_keys=True,indent=2)+'\n').encode()
def save(p,x): p.write_bytes(encoded(x))

assert not ART.exists() and not BACK.exists()
receipts={}
for name,job in [('c2',426),('a2',427),('g2',428)]:
    rec=json.loads((OUT/name/'receipt.json').read_text())
    assert int(rec['job'])==job and rec['passed'] and rec['inputs_unchanged']
    assert rec['frozen_root_verified'] and rec['binary_unchanged']
    for p,h in rec['input_hashes'].items(): assert sha(pathlib.Path(p))==h,p
    receipts[name]=rec
assert receipts['g2']['state_parity']
for rel in ['include/pf_coupled/gas_step_scope.hpp','tests/gas_step_scope.inc']:
    assert (R/'integration'/rel).read_bytes()==(REPO/'native/coupled_pic'/rel).read_bytes()

prefixes=['state ','motion_state ','publication_state ','binding_state ']
def states(p):
    lines=p.read_text().splitlines()
    return {prefix:next(line for line in lines if line.startswith(prefix)) for prefix in prefixes}
parity=[]
for case,mode,oldcase in [('c2','host','h2'),('a2','host','s2'),('g2','host','g1'),('g2','device','g1'),('g2','verified','g1')]:
    old=BASE/'evidence'/oldcase/(mode+'.log');new=OUT/case/(mode+'.log')
    assert states(old)==states(new),(case,mode)
    parity.append({'baseline_log':str(old),'baseline_sha256':sha(old),'new_log':str(new),'new_sha256':sha(new),'states':states(new)})

ART.mkdir()
for name in ['include','src','device','tests','tools','baseline','reference-v1.1','integration']:
    shutil.copytree(R/name,ART/name)
for name in ['BASELINE-v1.4.json','ROLLBACK-CONTRACT.txt','preserve-rollback.py']:
    shutil.copy2(R/name,ART/name)
for case in ['c1','c2','a2','g2']:
    d=ART/'evidence'/case;d.mkdir(parents=True)
    for p in (OUT/case).iterdir():
        if p.suffix in ('.json','.log'):shutil.copy2(p,d/p.name)
for job in [424,426,427,428]:
    p=pathlib.Path('/home/treefor/.local/state/simqueue/logs')/f'{job:06}.log'
    shutil.copy2(p,ART/'evidence'/p.name)
save(ART/'BASELINE-PARITY.json',parity)
sources={str(p.relative_to(ART)):sha(p) for name in ['include','src','device','tests','tools','baseline','reference-v1.1','integration'] for p in sorted((ART/name).rglob('*')) if p.is_file()}
save(ART/'SOURCE.json',{'identity':hashlib.sha256(encoded(sources)).hexdigest(),'files':sources})
checks={case:{p.stem:p.read_text().splitlines()[-1] for p in (OUT/case).glob('*.log') if p.stem in ['host','device','verified']} for case in receipts}
save(ART/'RESULTS.json',{'jobs':{'normal':426,'sanitized':427,'gpu':428},'checks':checks,'baseline_state_parity':True,'joint_publication_qualified':False,'physical_breakdown_qualified':False,'speedup_qualified':False,'maintained_source':'native/coupled_pic','failed_initial_job':424})
files={str(p.relative_to(ART)):sha(p) for p in sorted(ART.rglob('*')) if p.is_file()}
save(ART/'SEAL.json',{'identity':hashlib.sha256(encoded(files)).hexdigest(),'files':files})
shutil.copytree(ART,BACK)
for p in ART.rglob('*'):
    if p.is_file():assert p.read_bytes()==(BACK/p.relative_to(ART)).read_bytes()
print(json.dumps({'artifact':str(ART),'backup':str(BACK),'files':len(files)+1,'checks':checks,'seal_sha256':sha(ART/'SEAL.json')},indent=2))
