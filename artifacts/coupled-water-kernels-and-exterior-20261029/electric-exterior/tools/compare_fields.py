#!/usr/bin/env python3
"""Streaming fixed-state electric-field comparisons; no solve or trajectory."""
import argparse,math,pathlib
from static_contract import O,F,FIELDS,PARTICLES,read_json,write_json,record,sha

class Sum:
 def __init__(self):self.s=0.;self.c=0.
 def add(self,v):
  t=self.s+v
  self.c+=(self.s-t)+v if abs(self.s)>=abs(v) else (v-t)+self.s
  self.s=t
 def value(self):return self.s+self.c

def masks(x,y,z):
 inside=lambda l,r,b,t:l<=x<=r and b<=y<=t
 vacuum=not (inside(.002,.007,.0018,.002) or inside(.002,.007,.007,.0072))
 corner_distance=min(math.hypot(x-a,y-b) for a in [.002,.007] for b in [.0018,.002,.007,.0072])
 return dict(all=True,finite_width_gap=.002<=x<=.007 and .002<=y<=.007,above_upper_back=.002<=x<=.007 and y>.0072,foil_corner_vacuum_250um=vacuum and corner_distance<=.00025,midplane_slab=abs(z)<.0005)

class Stats:
 def __init__(self):
  self.n=0;self.parallel_n=0;self.s={k:Sum() for k in ['w','E2a','E2b','dE2','wp','para_a','para_b','abspara_a','abspara_b','para2a','para2b','perp2a','perp2b','dpara2','dperp2']};self.maxdelta=0.;self.maxcomponents=[0.]*3
 def add(self,w,a,b,B):
  self.n+=1;d=[bb-aa for aa,bb in zip(a,b)];d2=math.fsum(v*v for v in d)
  for key,value in [('w',w),('E2a',w*math.fsum(v*v for v in a)),('E2b',w*math.fsum(v*v for v in b)),('dE2',w*d2)]:self.s[key].add(value)
  self.maxdelta=max(self.maxdelta,math.sqrt(d2));self.maxcomponents=[max(old,abs(v)) for old,v in zip(self.maxcomponents,d)]
  norm=math.sqrt(math.fsum(v*v for v in B))
  if norm==0:return
  self.parallel_n+=1;bh=[v/norm for v in B];pa=math.fsum(v*q for v,q in zip(a,bh));pb=math.fsum(v*q for v,q in zip(b,bh));dp=pb-pa
  aperp=[v-pa*q for v,q in zip(a,bh)];bperp=[v-pb*q for v,q in zip(b,bh)];dperp=[v-u for u,v in zip(aperp,bperp)]
  for key,value in [('wp',w),('para_a',w*pa),('para_b',w*pb),('abspara_a',w*abs(pa)),('abspara_b',w*abs(pb)),('para2a',w*pa*pa),('para2b',w*pb*pb),('perp2a',w*math.fsum(v*v for v in aperp)),('perp2b',w*math.fsum(v*v for v in bperp)),('dpara2',w*dp*dp),('dperp2',w*math.fsum(v*v for v in dperp))]:self.s[key].add(value)
 def result(self):
  values={k:v.value() for k,v in self.s.items()};w=values['w'];wp=values['wp']
  def rms(key,weight):return math.sqrt(max(0.,values[key]/weight)) if weight else None
  return dict(macroparticles=self.n,physical_particle_weight=w,weighted_RMS_E_baseline_V_m=rms('E2a',w),weighted_RMS_E_candidate_V_m=rms('E2b',w),weighted_RMS_delta_E_V_m=rms('dE2',w),weighted_relative_L2_delta_E=math.sqrt(values['dE2']/values['E2a']) if values['E2a'] else None,max_delta_E_V_m=self.maxdelta,max_abs_delta_components_V_m=self.maxcomponents,nonzero_B_macroparticles=self.parallel_n,zero_B_macroparticles=self.n-self.parallel_n,parallel_physical_weight=wp,weighted_mean_Eparallel_baseline_V_m=values['para_a']/wp if wp else None,weighted_mean_Eparallel_candidate_V_m=values['para_b']/wp if wp else None,weighted_mean_abs_Eparallel_baseline_V_m=values['abspara_a']/wp if wp else None,weighted_mean_abs_Eparallel_candidate_V_m=values['abspara_b']/wp if wp else None,weighted_RMS_Eparallel_baseline_V_m=rms('para2a',wp),weighted_RMS_Eparallel_candidate_V_m=rms('para2b',wp),weighted_RMS_Eperpendicular_baseline_V_m=rms('perp2a',wp),weighted_RMS_Eperpendicular_candidate_V_m=rms('perp2b',wp),weighted_RMS_delta_Eparallel_V_m=rms('dpara2',wp),weighted_RMS_delta_Eperpendicular_V_m=rms('dperp2',wp))

def compare(particles,baseline,candidate,count):
 particles,baseline,candidate=map(pathlib.Path,[particles,baseline,candidate]);assert count<=46726 and count>0
 assert particles.stat().st_size==count*PARTICLES.size,'Missing/truncated/extra particles'
 assert baseline.stat().st_size==candidate.stat().st_size==count*FIELDS.size,'Missing/truncated/extra field records'
 pins=[record(p) for p in [particles,baseline,candidate]];stats={k:Stats() for k in masks(0,0,0)}
 with particles.open('rb') as ps,baseline.open('rb') as aa,candidate.open('rb') as bb:
  for _ in range(count):
   p=PARTICLES.unpack(ps.read(PARTICLES.size));ab=aa.read(FIELDS.size);bbuf=bb.read(FIELDS.size);a=FIELDS.unpack(ab);b=FIELDS.unpack(bbuf)
   assert all(math.isfinite(v) for v in p[1:]+a+b),'Nonfinite particle/field/potential'
   assert p[8]>0 and p[9]>0,'Invalid captured mass/weight'
   assert ab[24:48]==bbuf[24:48],'B is not bitwise fixed (signed zero included)'
   for key,included in masks(*p[1:4]).items():
    if included:stats[key].add(p[9],a[:3],b[:3],a[3:6])
 for row in pins:assert sha(row['path'])==row['sha256'],'Input changed while streaming'
 return dict(B_bitwise_fixed=True,input_receipts=pins,regions={k:s.result() for k,s in stats.items()})

def default_plan():
 old=O.parent/'cathode-emission-audit-20261029/static-domain-extended';new=O/'parent-static-runs';pairs=[]
 def add(name,a,b,group):pairs.append(dict(name=name,baseline=str(a),candidate=str(b),group=group))
 for width,n in [(13,32),(17,40),(25,48)]:add(f'uniform{width}_vs_graded{width}',old/f'{width}mm',new/f'graded-{width}-{n}','old_uniform_vs_graded_same_extent')
 counts=[40,60,80,96,110,120]
 for a,b in zip(counts,counts[1:]):add(f'nested_{a}_to_{b}',new/f'nested-105-{a}',new/f'nested-105-{b}','nested_extent_fixed_shared_interior_mesh')
 add('nested_first_to_last',new/'nested-105-40',new/'nested-105-120','nested_extent_fixed_shared_interior_mesh')
 for a,b in [(80,120),(120,160),(80,160)]:add(f'width80_margin_{a}_to_{b}',new/f'width80-{a}',new/f'width80-{b}','fixed_extent_outer_mesh_refinement')
 return dict(fixture=str(F),particles=46726,pairs=pairs)

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=pathlib.Path);parser.add_argument('--write-default-plan',type=pathlib.Path);parser.add_argument('--pair',action='append',default=[]);parser.add_argument('--out',type=pathlib.Path);args=parser.parse_args()
 if args.write_default_plan:
  assert not args.write_default_plan.exists();write_json(args.write_default_plan,default_plan());return
 assert args.plan and args.out and not args.out.exists()
 plan=read_json(args.plan);fixture=pathlib.Path(plan['fixture']);assert fixture==F and plan['particles']==46726
 selected=plan['pairs'] if not args.pair else [row for row in plan['pairs'] if row['name'] in args.pair]
 assert selected and len(selected)<=24
 if args.pair:assert set(args.pair)=={row['name'] for row in selected}
 result=dict(status='fixed_state_field_comparison_only',absolute_time_ns=read_json(F/'manifest.json')['absolute_native_time_s']*1e9,B_bitwise_fixed=True,source_particles=record(F/'particles.bin'),plan=record(args.plan),mask_definitions=dict(all='All captured particles',finite_width_gap='x2..7mm,y2..7mm',above_upper_back='x2..7mm,y>7.2mm',foil_corner_vacuum_250um='Outside both finite foil rectangles; Euclidean distance<=250um from any of their eight2Dcorners',midplane_slab='abs(z)<0.5mm'),parallel_definition='E dot B/|B|; Eperp=E-Eparallel B/|B|; zero-B samples counted/excluded from parallel/perpendicular statistics',comparisons=[],limits=['No Poisson or trajectory evolution','Masks are diagnostics and cannot establish emission origin','Finite successive extent differences do not bound the infinite-domain error','Extent changes and fixed-extent exterior refinement reported separately','No density output; fields V/m,B Tesla,physical particle weights dimensionless'])
 for pair in selected:
  a,b=pathlib.Path(pair['baseline']),pathlib.Path(pair['candidate'])
  metadata=[]
  for folder in [a,b]:
   r=read_json(folder/'result.json');assert r['mode']=='gpu' and r['passed'] and r['particles']==46726 and r['particle_state_unchanged'] and r['no_particle_advance'] and r['factorizations']==1 and r['field_solves']==2
   layout=read_json(folder/'work-and-layout.json') if (folder/'work-and-layout.json').exists() else None
   metadata.append(dict(folder=str(folder),result=record(folder/'result.json'),layout=layout))
  compared=compare(F/'particles.bin',a/'particle_fields.bin',b/'particle_fields.bin',46726)
  result['comparisons'].append(dict(**pair,metadata=metadata,**compared))
 write_json(args.out,result);print(f"Compared{len(selected)} fixed-state field pairs; B bitwise fixed; no evolution")

if __name__=='__main__':main()
