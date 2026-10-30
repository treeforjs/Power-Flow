#!/usr/bin/env python3
"""Offline fixture provenance and binary parameter conversion; no CUDA calls."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import struct

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('fixture',type=Path)
p.add_argument('--allow-synthetic',action='store_true')
a=p.parse_args();m=json.loads((a.fixture/'manifest.json').read_text())
real=m.get('real_midpoint_capture') is True and m.get('synthetic_self_fields') is False
if not real and not a.allow_synthetic:p.error('Actual midpoint capture is required')
if real:
    assert m['field_callbacks_added']==0 and m['capture_stage']=='post_half_drift_post_first_guarded_solve_pre_kick'
    assert m['native_provider_last_time_s']==m['absolute_native_time_s']==m['voltage_provider_last_time_s']
    assert m['native_lower_time_s']<=m['absolute_native_time_s']<=m['native_upper_time_s']
    before=(a.fixture/'particles.bin').read_bytes();after=(a.fixture/'particles_after_kick.bin').read_bytes()
    reference=(a.fixture/'reference.bin').read_bytes();inputs=(a.fixture/'inputs.bin').read_bytes();n=m['particles']
    assert len(before)==len(after)==len(reference)==80*n and len(inputs)==56*n
    for i in range(n):
        b,c,r,v=before[80*i:80*(i+1)],after[80*i:80*(i+1)],reference[80*i:80*(i+1)],inputs[56*i:56*(i+1)]
        assert b[:32]==c[:32] and b[56:]==c[56:]
        assert c[32:56]==r[48:72] and b[32:56]==v[24:48] and b[8:32]==v[:24]
        assert struct.unpack('<Q',r[72:])[0]==0 and all(map(math.isfinite,struct.unpack('<9d',r[:72])))
params=struct.pack('<QQ4i7d',0x5049434241544348,1 if real else 0,
    *[m[x] for x in ('self_nx','self_ny','native_nx','native_ny')],
    *[m[x] for x in ('width','height','depth','short_z','fraction','voltage','dt')])
target=a.fixture/'parameters.bin'
if target.exists() and target.read_bytes()!=params:raise RuntimeError('Existing binary parameters disagree with manifest')
target.write_bytes(params)
hashes={x.name:hashlib.sha256(x.read_bytes()).hexdigest() for x in a.fixture.iterdir() if x.is_file() and x.name!='sha256.json'}
(a.fixture/'sha256.json').write_text(json.dumps(hashes,indent=2)+'\n')
print(json.dumps({'passed':True,'real_midpoint':real,'gpu_executed':False,'parameters_bytes':len(params),'fixture':str(a.fixture)},indent=2))
