#!/usr/bin/python -I
"""Metadata-only gate/transition regressions; no build or numerical work."""
import hashlib,importlib.util,json,os
from pathlib import Path
p=Path(__file__).with_name('owner_timing.py')
s=importlib.util.spec_from_file_location('owner_timing_check',p)
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
os.environ.update(SIMQ_JOB_ID='999999',SIMQ_UNIT='simq-check.service',SIMQ_MEM=str(4*1024**3))
base={'id':999999,'state':'running','partition':'powerflow','borrow':0,'exclusive_timing':1,
      'unit_name':'simq-check.service','cpu_list':'6','mem_bytes':4*1024**3,'request_sha256':'example-request'}
passed=[]
def refuses(name,fn,text):
    try:fn()
    except RuntimeError as e:
        assert text in str(e),(name,str(e));passed.append(name)
    else:raise AssertionError('accepted invalid '+name)
m.check_job(base,6);passed.append('exact_running_metadata')
for key,bad in [('id',999998),('state','queued'),('partition','raikou-mfem'),('borrow',1),
                ('exclusive_timing',0),('unit_name','other.service'),('cpu_list','7'),('mem_bytes',0)]:
    row={**base,key:bad}
    refuses('wrong_'+key,lambda row=row:m.check_job(row,6),'live exclusive')
starting={**base,'state':'starting'}
refuses('starting_never_admitted_to_frozen_phase',lambda:m.check_job(starting,6),'live exclusive')
actual=({},6,'0::/simq-check.service','/sys/fs/cgroup/simq-check.service/memory.max',str(4*1024**3))
class Clock:
    def __init__(self):self.t=0
    def now(self):return self.t
    def sleep(self,dt):self.t+=dt
def exercise(sequence,allocations=None,timeout=1):
    it=iter(sequence);ai=iter(allocations or [actual]*len(sequence));history=[];clock=Clock()
    result=m.await_running(lambda:next(it),lambda:next(ai),history,timeout,clock.now,clock.sleep)
    return result,history
(job,allocation),history=exercise([starting,starting,base])
assert job==base and len(history)==3 and allocation==actual
passed.append('starting_transition_waits_for_running')
refuses('bounded_starting_timeout',lambda:exercise([starting],timeout=0),'bounded wait')
refuses('changed_request_during_transition',lambda:exercise([starting,{**base,'request_sha256':'other'}]),'allocation changed')
refuses('changed_cgroup_during_transition',lambda:exercise([starting,base],[actual,({},6,'0::/other.service',actual[3],actual[4])]),'allocation changed')
refuses('cancellation_during_transition',lambda:exercise([starting,{**base,'state':'cancelling'}]),'live exclusive')
print(json.dumps({'schema':'owner-timing-metadata-checks-v2','launcher_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),
                  'metadata_only':True,'passed':passed,'count':len(passed)},indent=2))
