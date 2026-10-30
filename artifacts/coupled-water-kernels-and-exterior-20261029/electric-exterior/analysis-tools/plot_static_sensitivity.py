#!/usr/bin/env python3
"""White scientific figures from completed fixed-state comparison receipts."""
import argparse,hashlib,json,math,pathlib

REGIONS=[('all','All particles','#007c83'),('finite_width_gap','Finite-width gap','#0072b2'),('above_upper_back','Above upper back','#d55e00'),('foil_corner_vacuum_250um','250 µm corner vacuum','#7b3294')]
METRICS=[('weighted_RMS_delta_E_V_m','Weighted RMS |ΔE|'),('weighted_RMS_delta_Eparallel_V_m','Weighted RMS ΔE∥'),('weighted_RMS_delta_Eperpendicular_V_m','Weighted RMS |ΔE⊥|')]
GROUPS=['nested_extent_fixed_shared_interior_mesh','fixed_extent_outer_mesh_refinement','old_uniform_vs_graded_same_extent']

def sha(path):return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()
def finite(v):
 if isinstance(v,float):assert math.isfinite(v),'Nonfinite comparison metric'
 elif isinstance(v,dict):
  for value in v.values():finite(value)
 elif isinstance(v,list):
  for value in v:finite(value)
def load(path):
 result=json.loads(pathlib.Path(path).read_text());finite(result);return result
def geometry(row,index):
 m=row['metadata'][index];p=pathlib.Path(m['result']['path']);assert sha(p)==m['result']['sha256'],'Result receipt changed'
 r=load(p);layout=m['layout']
 if layout:
  bounds=layout['electric_bounds_m'];width=1000*(bounds[1]-bounds[0]);margin=layout['margin_cells_each_side']
 else:width=float(r['case_full_width_mm']);margin=None
 return dict(width_mm=width,margin=margin,nodes=r['nodes'])
def selected_rows(data):
 rows={group:[] for group in GROUPS}
 for row in data['comparisons']:
  assert row['B_bitwise_fixed']
  if row['group'] not in rows:continue
  if row['name'] in ['nested_first_to_last','width80_margin_80_to_160']:continue
  a,b=geometry(row,0),geometry(row,1)
  for key,_,_ in REGIONS:
   region=row['regions'][key]
   for metric,_ in METRICS:
    value=region[metric];assert value is None or value>=0
  if row['group']=='nested_extent_fixed_shared_interior_mesh':
   assert b['width_mm']>a['width_mm'] and a['margin']<b['margin']
  elif row['group']=='fixed_extent_outer_mesh_refinement':
   assert abs(a['width_mm']-80)<1e-10 and abs(b['width_mm']-80)<1e-10 and a['margin']<b['margin']
  else:assert abs(a['width_mm']-b['width_mm'])<1e-10
  rows[row['group']].append(dict(data=row,baseline=a,candidate=b))
 for group,values in rows.items():
  assert values,f'Missing comparison group: {group}'
  values.sort(key=lambda v:v['candidate']['margin'] if group=='fixed_extent_outer_mesh_refinement' else v['candidate']['width_mm'])
 return rows

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--comparison',type=pathlib.Path,required=True);parser.add_argument('--out-dir',type=pathlib.Path,required=True);args=parser.parse_args()
 assert not args.out_dir.exists(),'Fresh plotting directory required'
 data=load(args.comparison);assert data['status']=='fixed_state_field_comparison_only' and data['B_bitwise_fixed']
 rows=selected_rows(data);input_hash=sha(args.comparison)
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 plt.rcParams.update({'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white','font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.labelcolor':'#171717','text.color':'#171717','xtick.color':'#171717','ytick.color':'#171717','axes.grid':True,'grid.color':'#d9d9d9','grid.linewidth':.6,'grid.alpha':.8,'axes.axisbelow':True,'svg.fonttype':'none'})
 args.out_dir.mkdir();time_ns=data['absolute_time_ns']
 footer=f'Frozen charge and particle positions at {time_ns:.9f} ns; no trajectories or infinite-domain convergence claim.\nGap/back/corner masks include all z positions; regions overlap.'
 outputs=[]
 def finish(fig,name,title):
  fig.suptitle(title,fontsize=13,y=.99)
  fig.text(.5,.015,footer,ha='center',fontsize=8,color='#333333')
  fig.tight_layout(rect=(0,.14,1,.88))
  for extension in ['png','svg']:
   p=args.out_dir/(name+'.'+extension);fig.savefig(p,dpi=180);outputs.append(dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)))
  plt.close(fig)
 def delta_figure(group,name,title):
  values=rows[group];fig,axes=plt.subplots(1,3,figsize=(13.2,4.8))
  if group=='fixed_extent_outer_mesh_refinement':x=list(range(len(values)));labels=[f"{v['baseline']['margin']}→{v['candidate']['margin']}" for v in values];xlabel='Margin cells per side (baseline → new)'
  else:x=[v['candidate']['width_mm'] for v in values];labels=None;xlabel='New full electric box width (mm)'
  for ax,(metric,label) in zip(axes,METRICS):
   all_y=[]
   for region,legend,color in REGIONS:
    y=[v['data']['regions'][region][metric] for v in values]
    assert all(value is not None for value in y),'No nonzero-B directional samples in a required region'
    y=[v*1e-6 for v in y];all_y.extend(y);ax.plot(x,y,'o-',color=color,label=legend,lw=1.7,ms=4)
   if all(v>0 for v in all_y) and max(all_y)/min(all_y)>20:ax.set_yscale('log')
   ax.set_ylabel(label+' (MV/m)');ax.set_xlabel(xlabel)
   if labels:ax.set_xticks(x,labels)
   if group=='nested_extent_fixed_shared_interior_mesh':ax.set_xscale('log');ax.set_xticks(x,[f'{v:.1f}' for v in x]);ax.minorticks_off()
  axes[0].legend(frameon=False,fontsize=8)
  
  if group=='nested_extent_fixed_shared_interior_mesh':
   chain=[values[0]['baseline']['width_mm']]+[v['candidate']['width_mm'] for v in values]
   fig.text(.5,.075,'Each point is the change from the preceding extent: '+ ' → '.join(f'{v:.3f}' for v in chain)+' mm. Actual numeric x positions use a log scale.',ha='center',fontsize=8,color='#333333')
  finish(fig,name,title)
 delta_figure(GROUPS[0],'extent-summary','Static frozen-charge sensitivity: successive boundary extensions\nShared interior cells; 5% outward cell growth')
 delta_figure(GROUPS[1],'outer-refinement-summary','Static frozen-charge sensitivity: refine exterior cells at fixed 80 mm width')
 delta_figure(GROUPS[2],'uniform-vs-graded-summary','Static frozen-charge control: uniform versus graded mesh at the same extent')
 # Absolute gap/back levels distinguish a small change from a small field.
 extent=rows[GROUPS[0]];fig,axes=plt.subplots(1,3,figsize=(13.2,4.3))
 absolute=[('weighted_RMS_E_baseline_V_m','weighted_RMS_E_candidate_V_m','Weighted RMS |E|'),('weighted_mean_abs_Eparallel_baseline_V_m','weighted_mean_abs_Eparallel_candidate_V_m','Weighted mean |E∥|'),('weighted_RMS_Eperpendicular_baseline_V_m','weighted_RMS_Eperpendicular_candidate_V_m','Weighted RMS |E⊥|')]
 for region,legend,color in REGIONS[1:3]:
  x=[extent[0]['baseline']['width_mm']]+[v['candidate']['width_mm'] for v in extent]
  for ax,(first,following,label) in zip(axes,absolute):
   y=[extent[0]['data']['regions'][region][first]]+[v['data']['regions'][region][following] for v in extent]
   assert all(v is not None for v in y);ax.plot(x,[v*1e-6 for v in y],'o-',color=color,label=legend,lw=1.8,ms=4);ax.set_xscale('log');ax.set_xticks(x,[f'{v:.1f}' for v in x]);ax.minorticks_off();ax.set_xlabel('Full electric box width (mm)');ax.set_ylabel(label+' (MV/m)')
 axes[0].legend(frameon=False);finish(fig,'gap-back-field-levels','Gap and back: absolute field levels in the same frozen charge state')
 assert sha(args.comparison)==input_hash,'Comparison changed while plotting'
 receipt=dict(status='white_scientific_fixed_state_plots',comparison=dict(path=str(args.comparison),sha256=input_hash),plot_script=dict(path=str(pathlib.Path(__file__).resolve()),sha256=sha(__file__)),outputs=outputs,groups_separate=GROUPS,geometry=[dict(group=g,baseline=v['baseline'],candidate=v['candidate']) for g,values in rows.items() for v in values],mask_scope='Gap/back/corner masks include all captured z positions; midplane is a separate diagnostic in the comparison JSON. Masks overlap and do not establish emission origin.',field_units='MV/m on plots; exact conversion1MV/m=1e6V/m; no magnitude normalization',density_generated=False,no_Poisson_or_trajectory_execution=True,convergence_claim=False)
 (args.out_dir/'plot-receipt.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
 print(f'Wrote {len(outputs)} white scientific plot files; no convergence claim')

if __name__=='__main__':main()
