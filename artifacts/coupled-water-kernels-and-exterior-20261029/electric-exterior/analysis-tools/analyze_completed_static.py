#!/usr/bin/env python3
"""Bounded independent norm/axis checks and extra pilot extent comparison."""
import pathlib,sys,array,math,struct
O=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(O/'tools'))
from static_contract import F,FIELDS,PARTICLES,read_json,write_json,record,sha,EPS
from compare_fields import compare,masks

def independent_pair(pair):
 a=pathlib.Path(pair['baseline'])/'particle_fields.bin';b=pathlib.Path(pair['candidate'])/'particle_fields.bin'
 terms={key:{'w':[],'total':[],'base':[],'parallel':[],'perp':[]} for key in pair['regions']}
 with (F/'particles.bin').open('rb') as ps,a.open('rb') as aa,b.open('rb') as bb:
  for _ in range(46726):
   p=PARTICLES.unpack(ps.read(PARTICLES.size));x=FIELDS.unpack(aa.read(FIELDS.size));y=FIELDS.unpack(bb.read(FIELDS.size));de=[y[k]-x[k] for k in range(3)];B=x[3:6];bn=math.sqrt(math.fsum(v*v for v in B));d2=math.fsum(v*v for v in de);e2=math.fsum(v*v for v in x[:3]);w=p[9]
   for key,include in masks(*p[1:4]).items():
    if not include:continue
    t=terms[key];t['w'].append(w);t['total'].append(w*d2);t['base'].append(w*e2)
    if bn:
     dp=math.fsum(de[k]*B[k] for k in range(3))/bn;t['parallel'].append(w*dp*dp);t['perp'].append(w*max(0.,d2-dp*dp))
 result={}
 for key,t in terms.items():
  r=pair['regions'][key];w=math.fsum(t['w']);assert len(t['w'])==r['macroparticles'] and w==r['physical_particle_weight']
  total=math.sqrt(math.fsum(t['total'])/w);para=math.sqrt(math.fsum(t['parallel'])/r['parallel_physical_weight']);perp=math.sqrt(math.fsum(t['perp'])/r['parallel_physical_weight']);base=math.sqrt(math.fsum(t['base'])/w)
  # A numerical comparison bound only:64eps times the baseline evaluation
  # scale covers the different dot/subtraction/norm operation ordering.
  bound=64*EPS*max(base,1.)
  errors=[abs(v-r[k]) for v,k in [(total,'weighted_RMS_delta_E_V_m'),(para,'weighted_RMS_delta_Eparallel_V_m'),(perp,'weighted_RMS_delta_Eperpendicular_V_m')]]
  assert max(errors)<=bound
  result[key]=dict(independent_RMS_delta_E_V_m=total,independent_RMS_delta_Eparallel_V_m=para,independent_RMS_delta_Eperp_V_m=perp,max_abs_difference_from_streaming_metrics_V_m=max(errors),operation_scale_check_bound_V_m=bound,passed=True)
 return result

def main():
 data=read_json(O/'field-comparison.json');assert len(data['comparisons'])==12
 tail=next(row for row in data['comparisons'] if row['name']=='nested_110_to_120');independent=independent_pair(tail)
 axes={}
 for axis in ['x','y']:
  prev=None;pn=None;bounds=[]
  for n in [40,60,80,96,110,120]:
   path=O/f'nested-105-{n}-sealed/axes_{axis}.bin';a=array.array('d');a.frombytes(path.read_bytes())
   if prev is not None:assert a[n-pn:n-pn+len(prev)]==prev
   bounds.append(dict(margin=n,low_m=a[0],high_m=a[-1],full_width_mm=(a[-1]-a[0])*1000,record=record(path)));prev=a;pn=n
  axes[axis]=dict(shared_interior_bitwise_equal=True,cases=bounds)
 extra=[]
 for n in [110,120]:
  a=O/'parent-static-runs/width80-160';b=O/f'parent-static-runs/nested-105-{n}'
  result=compare(F/'particles.bin',a/'particle_fields.bin',b/'particle_fields.bin',46726)
  extra.append(dict(name=f'width80_margin160_vs_nested{n}',baseline=str(a),candidate=str(b),**result))
 write_json(O/'pilot-extra-static-comparisons.json',dict(status='extra_frozen_state_comparisons_no_evolution',comparisons=extra))
 summary=[]
 for row in data['comparisons']:
  entry=dict(name=row['name'],group=row['group'],regions={})
  for key in ['all','finite_width_gap','above_upper_back','foil_corner_vacuum_250um']:
   r=row['regions'][key];entry['regions'][key]={k:r[k] for k in ['macroparticles','physical_particle_weight','weighted_RMS_delta_E_V_m','weighted_relative_L2_delta_E','max_delta_E_V_m','weighted_RMS_delta_Eparallel_V_m','weighted_RMS_delta_Eperpendicular_V_m']}
  summary.append(entry)
 write_json(O/'scientific-comparison-audit.json',dict(status='CPU_scientific_masks_norms_and_nested_metadata_reviewed',comparison=record(O/'field-comparison.json'),independent_tail_math_fsum_checks=independent,nested_axes=axes,summary=summary,extra_comparison=record(O/'pilot-extra-static-comparisons.json'),masks_include_all_z_except_midplane=True,masks_overlap=True,corner_mask_is_2D_corners_not_all_foil_faces_or_3D_endfaces=True,zero_B_samples_in_tail=tail['regions']['all']['zero_B_macroparticles'],relative_L2_denominator='Weighted baseline total E magnitude; can be dominated by driven bias. Absolute parallel/perpendicular differences and source guards remain separate.',no_new_residual_tolerance=True,no_Poisson_or_GPU_or_trajectory_execution=True,no_infinite_domain_or_evolving_plasma_convergence_claim=True))
 print('Exact full widths(mm):',[row['full_width_mm'] for row in axes['x']['cases']])
 for row in extra:print(row['name'],{key:{'RMS_V_m':row['regions'][key]['weighted_RMS_delta_E_V_m'],'relative':row['regions'][key]['weighted_relative_L2_delta_E']} for key in ['all','finite_width_gap','above_upper_back','foil_corner_vacuum_250um']})

if __name__=='__main__':main()
