#!/usr/bin/env python3
"""Read-only atlas/actual-like census preparation; no gas reconstruction."""
from fractions import Fraction
import argparse,hashlib,json,os,pathlib,resource,signal,struct,time
from microcase import Grid
FIX=pathlib.Path('/home/treefor/sim-workspace/powerflow/outputs/pic-profile/cuda-gather-investigation/actual-midpoint-fixture')
PT=struct.Struct('<Q9d');PK=struct.Struct('<QQ5d');CELL=struct.Struct('<Qd')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--out',type=pathlib.Path,required=True);a=p.parse_args();assert not a.out.exists() and os.sched_getaffinity(0)=={7} and os.getpriority(os.PRIO_PROCESS,0)>=15;resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));signal.alarm(90);start=time.monotonic();a.out.mkdir()
 pins={name:h for name,h in json.loads((FIX/'sha256.json').read_text()).items()};assert all(sha(FIX/name)==h for name,h in pins.items());raw=(FIX/'particles.bin').read_bytes();assert len(raw)==46726*PT.size;short=struct.unpack_from('<d',(FIX/'options.bin').read_bytes(),152)[0];g=Grid();occupied={};rows=[]
 for item in PT.iter_unpack(raw):
  ident,x,y,z,ux,uy,uz,q,m,w=item;assert w>0 and q==-1.602176634e-19
  # Exact binary-double operands are preserved; no spatial rounding/quantization.
  point=(Fraction(x)*10**9,Fraction(y)*10**9,(Fraction(z)-Fraction(short))*10**9);cid=g.locate(point);occupied[cid]=occupied.get(cid,0)+1;rows.append(PK.pack(ident,cid,w,m,ux,uy,uz))
 (a.out/'projectiles.bin').write_bytes(b''.join(rows))
 with (a.out/'atlas-cells.bin').open('wb') as stream:
  for cid in range(g.total):stream.write(CELL.pack(cid,float(g.volume(cid))))
 with (a.out/'occupied-cells.bin').open('wb') as stream:
  for cid in sorted(occupied):stream.write(CELL.pack(cid,float(g.volume(cid))))
 assert sum(occupied.values())==46726 and all(sha(FIX/name)==h for name,h in pins.items())
 receipt=dict(schema='collision_atlas_actual_like_input/v1',geometry_identity=g.identity,geometry_bins=g.total,occupied_cells=len(occupied),projectile_carriers=len(rows),all_original_positive_carriers_retained=True,fixture_hashes=pins,position_policy='Original x/y binary operands; exact rational z-short_z frame translation for this read-only manufactured atlas test; no native position is mutated/quantized',target_policy='C++ test adds one explicitly manufactured stationary weight8 target per occupied bin; no real gas/ion distribution is inferred',no_physical_rate_or_evolution_claim=True,files={p.name:sha(p) for p in sorted(a.out.iterdir()) if p.is_file()},execution=dict(wall_s=time.monotonic()-start,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,affinity=sorted(os.sched_getaffinity(0)),nice=os.getpriority(os.PRIO_PROCESS,0),time_cap_s=90,address_space_cap_bytes=4*1024**3))
 (a.out/'input-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(dict(projectiles=len(rows),occupied_cells=len(occupied),geometry_bins=g.total,execution=receipt['execution']),indent=2))
if __name__=='__main__':main()
