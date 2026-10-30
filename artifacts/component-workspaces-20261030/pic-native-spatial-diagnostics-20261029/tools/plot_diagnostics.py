#!/usr/bin/env python3
import argparse,pathlib,json,hashlib,math
def main():
 p=argparse.ArgumentParser();p.add_argument('--diagnostics',type=pathlib.Path,required=True);p.add_argument('--out-dir',type=pathlib.Path,required=True);a=p.parse_args();assert not a.out_dir.exists();data=json.loads(a.diagnostics.read_text());assert data['schema']=='native_spatial_diagnostics/v1'
 import numpy as np
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 from matplotlib.colors import LogNorm
 plt.rcParams.update({'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white','font.size':10,'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
 a.out_dir.mkdir();files=[];omitted=[]
 def save(fig,name):
  for ext in ['png','svg']:
   path=a.out_dir/(name+'.'+ext);fig.savefig(path,dpi=170,facecolor='white');files.append({'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
  plt.close(fig)
 edges=np.asarray(data['energy_edges_eV']);fig,axes=plt.subplots(1,2,figsize=(12,4.4));tails=[]
 for state,h in data['energy_distributions'].items():
  counts=[v['rounded'] for v in h['bin_counts']];pdf=h['ensemble_EEDF_per_eV'];outside=[h['underflow_count']['rounded'],h['overflow_count']['rounded']];unresolved=h['unresolved_energy_count']['rounded']
  if any(v is None for v in counts+pdf+outside+[unresolved]):omitted.append(state+' energy plot: nonrepresentable double; exact receipt retained');continue
  axes[0].stairs(counts,edges,label=state,lw=1.7);axes[1].stairs(pdf,edges,label=state,lw=1.7);tails.append(f"{state}: outside bins={sum(outside):g}; unresolved energy={unresolved:g}")
 for ax in axes:ax.set_xscale('symlog',linthresh=1e-6);ax.set_yscale('log');ax.set_xlabel('Actual kinetic energy (eV)');ax.grid(alpha=.22);ax.legend(frameon=False)
 axes[0].set_ylabel('Physical particle count per bin');axes[1].set_ylabel('Ensemble energy PDF (eV⁻¹)');fig.suptitle(f"Native midpoint at {data['time_s']*1e9:.9f} ns: weighted energy histogram and EEDF")
 fig.text(.5,.025,'All species weight remains in the PDF denominator; empty bins are absent on the log display.\n'+'; '.join(tails)+'; no Maxwellian fit or physical breakdown qualification.',ha='center',fontsize=8);fig.tight_layout(rect=(0,.12,1,.92));save(fig,'energy-histogram-EEDF')
 g=data['diagnostic_geometry'];nx=len(g['xs_m'])-1;ny=len(g['ys_m'])-1;nz=len(g['zs_m'])-1
 if nz!=1:omitted.append('2D density plots require one declared z cell; all 3D cell densities remain in JSON')
 if nz==1:
  for state in data['species_counts']:
   if state in data['outside_diagnostic_grid_counts'] or state in data['zero_volume_cell_counts']:omitted.append(state+' density plot: positive weight has no physical diagnostic volume');continue
   rows=[r for r in data['species_cell_densities'] if r['state']==state];has_slab=data['observation_slab_m'] is not None
   if any(r['density_cm3'] is None or (has_slab and r['observation']['density_cm3'] is None) for r in rows):omitted.append(state+' density plot: nonrepresentable density or zero observation volume');continue
   full=np.zeros((ny,nx));slab=np.zeros((ny,nx))
   for row in rows:
    j,i=divmod(row['cell_id'],nx);full[j,i]=row['density_cm3']
    if has_slab:slab[j,i]=row['observation']['density_cm3']
   positive=np.concatenate([full[full>0],slab[slab>0]]) if has_slab else full[full>0]
   if not len(positive):continue
   lo=float(positive.min());hi=float(positive.max());hi=max(hi,lo*1.01);fig,axes=plt.subplots(1,2 if has_slab else 1,figsize=(10 if has_slab else 5,5.0),squeeze=False)
   titles=[f"Physical full-height volume ({(g['zs_m'][-1]-g['zs_m'][0])*1000:g} mm)",f"Observation slab volume ({(data['observation_slab_m'][1]-data['observation_slab_m'][0])*1000:g} mm)" if has_slab else '']
   for ax,values,title in zip(axes[0],[full,slab],titles):
    mesh=ax.pcolormesh(np.asarray(g['xs_m'])*1000,np.asarray(g['ys_m'])*1000,np.ma.masked_where(values==0,values),norm=LogNorm(lo,hi),cmap='viridis',shading='flat');ax.set_aspect('equal');ax.set_xlabel('PIC x (mm)');ax.set_ylabel('PIC y (mm)');ax.set_title(title);fig.colorbar(mesh,ax=ax,label=state+' number density (cm⁻³)')
    for box in g.get('solid_boxes_m',[]):ax.add_patch(plt.Rectangle((box[0]*1000,box[2]*1000),(box[1]-box[0])*1000,(box[3]-box[2])*1000,facecolor='#dddddd',edgecolor='#333333',lw=.6))
   fig.suptitle(f"{state} density at {data['time_s']*1e9:.9f} ns: actual weights, distinct sampling volumes");fig.text(.5,.035,'Common absolute cm⁻³ color scale. White vacuum pixels have zero represented weight; solids are gray.\nInstantaneous diagnostic only; face proximity does not identify source origin or return.',ha='center',va='bottom',fontsize=8);fig.tight_layout(rect=(0,.17,1,.92));save(fig,state+'-density-cm3')
 (a.out_dir/'plot-receipt.json').write_text(json.dumps({'input':str(a.diagnostics),'input_sha256':hashlib.sha256(a.diagnostics.read_bytes()).hexdigest(),'plot_source_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),'files':files,'omitted_plots':omitted,'density_units':'cm^-3','energy_bins_unrenormalized':True,'new_evolution_or_GPU':False,'full_mechanism_qualification':False},indent=2)+'\n')
 print('White diagnostic plots written; no physical validation claim')
if __name__=='__main__':main()
