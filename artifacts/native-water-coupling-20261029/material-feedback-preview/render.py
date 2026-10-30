from pathlib import Path
import csv,json,hashlib,os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
assert os.sched_getaffinity(0)=={7}
S=Path('/home/treefor/sim-workspace/powerflow/outputs/ss304-material-feedback-pilot-20261029')
O=Path(__file__).parent
results=json.loads((S/'results.json').read_text());source_hashes={}
plt.rcParams.update({'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(2,2,figsize=(12,8.2),layout='constrained')
colors=['#7c7d80','#df8c1b','#087f8c'];styles=[':','--','-']
for k,c in enumerate(results['cases']):
 f=S/c['case']/'trajectory.csv';digest=hashlib.sha256(f.read_bytes()).hexdigest()
 assert digest==c['trajectory_sha256'];source_hashes[str(f)]=digest
 with f.open() as h: rows=[{n:float(v) for n,v in r.items()} for r in csv.DictReader(h)]
 assert len(rows)==c['accepted_steps']+1
 t=[r['time_s']*1e9 for r in rows]
 label=['32 × 4; max Δt 250 ps','64 × 8; max Δt 250 ps','64 × 8; max Δt 125 ps'][k]
 for ax,name,scale in [(axes[0,0],'max_T_K',1),(axes[0,1],'sigma_min_S_m',1e-6),(axes[1,0],'thermal_energy_J',1),(axes[1,1],'max_liquid_fraction',1)]:
  ax.plot(t,[r[name]*scale for r in rows],styles[k],color=colors[k],label=label,lw=1.8)
 for ax in axes.flat:ax.set_xlim(0,50);ax.set_xlabel('Time (ns)');ax.grid(alpha=.18)
axes[0,0].set(title='Hottest accepted cell',ylabel='Temperature (K)',ylim=(250,3150))
axes[0,0].axhline(3000,color='#a23b3b',ls='--',lw=1)
axes[0,0].text(1,3040,'3000 K model limit',color='#a23b3b',fontsize=9)
axes[0,0].legend(loc='upper left',bbox_to_anchor=(0.0,.93),frameon=False,fontsize=9)
axes[0,1].set(title='Minimum conductivity used by the field solve',ylabel='Conductivity (MS/m)')
axes[1,0].set(title='Solid thermal energy gain',ylabel='Energy (J)')
axes[1,1].set(title='Maximum local liquid fraction',ylabel='Liquid fraction',ylim=(-.03,1.08))
fig.suptitle('SS304 heating with conductivity feedback — actual C19 drive',fontsize=15)
fig.supxlabel('Solid-only reduced MQS; fixed shape and mass; conditional liquid properties.\nTraces stop at their last accepted state. No plasma, vapor, mechanics, or spatial convergence claim.',fontsize=10)
fig.savefig(O/'material_feedback.png',dpi=170)
fig.savefig(O/'material_feedback.pdf')
receipt={'source_results_sha256':hashlib.sha256((S/'results.json').read_bytes()).hexdigest(),'verified_trajectories':source_hashes,'accepted_trace_only':True,'rejected_states_plotted':False,'new_simulation':False,'white_background':True,'files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [O/'material_feedback.png',O/'material_feedback.pdf']}}
(O/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':'rendered','figure':str(O/'material_feedback.png'),'trajectories':len(source_hashes)}))
