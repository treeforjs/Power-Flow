from pathlib import Path
import os,sys,json,csv,hashlib,time,resource,collections
OUT=Path(__file__).resolve().parent
os.environ['MPLCONFIGDIR']=str(OUT/'matplotlib-cache')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
W=Path('/home/treefor/Documents/Codex/2026-10-28/we')
F=W/'outputs/pic-profile/cuda-gather-investigation/actual-midpoint-fixture'
A=W/'outputs/pic-profile/em-field-consistency-audit/alignment'
E=W/'outputs/pic-profile/native-event-log-prototype/real-seed-1ps/events'
V=W/'outputs/pic-profile/native-event-analysis-v2'
R=W/'outputs/pic-profile/event-analysis-v2-independent-review'
inputs={};started=time.monotonic()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def bind(p,expected=None):
 p=Path(p);h=sha(p)
 if expected is not None:assert h==expected,str(p)
 inputs[str(p)]={'sha256':h,'bytes':p.stat().st_size};return p

def main():
 # Own preview outputs may be regenerated after visual layout inspection.
 fixture=json.loads(bind(F/'manifest.json').read_text());fs=json.loads(bind(F/'sha256.json').read_text())
 bind(F/'particles.bin',fs['particles.bin']);assert sha(F/'manifest.json')==fs['manifest.json']
 decoded=json.loads(bind(A/'decode-validation.json').read_text());bind(A/'decode_fixture.cpp');bind(A/'read-only-receipt.json')
 header=bind('/home/treefor/projects/runtimes/mfem-gpu-pic-20261029/source/include/raikou/mfem_gap_pic.h');bind('/home/treefor/projects/runtimes/pic-cuda-gather-investigation-20261029/source/snapshot_io.hpp')
 bind('/home/treefor/projects/runtimes/mfem-gpu-pic-20261029/source/examples/mfem_heating_history_pic.cpp')
 bind(W/'outputs/pic-profile/native-event-log-prototype/source/src/mfem_gap_pic.cpp')
 bind(W/'outputs/pic-profile/native-event-log-prototype/source/src/pic_native_events.hpp')
 assert decoded['GapParticle_bytes']==80 and decoded['GapParticle_weight_offset']==72 and decoded['native_types_used'] and decoded['all_reference_and_particle_identities_match']
 # Layout is the frozen header plus prior same-header/compiler C++ receipt,
 # then checked field-for-field against that typed decoder's exact output.
 dtype=np.dtype({'names':['id','position','u','q','mass','weight'],'formats':['<u8',('<f8',3),('<f8',3),'<f8','<f8','<f8'],'offsets':[0,8,32,56,64,72],'itemsize':80})
 p=np.fromfile(F/'particles.bin',dtype=dtype);assert len(p)==fixture['particles']==46726
 dc=np.genfromtxt(bind(A/'particles-with-fields.csv'),delimiter=',',names=True,dtype=None,encoding='utf-8')
 assert np.array_equal(p['id'],dc['id']) and np.array_equal(p['weight'],dc['weight'])
 for i,k in enumerate(('x_m','y_m','z_m')):assert np.array_equal(p['position'][:,i],dc[k])
 for i,k in enumerate(('ux_m_s','uy_m_s','uz_m_s')):assert np.array_equal(p['u'][:,i],dc[k])
 assert np.array_equal(p['q'],dc['charge_C']) and np.array_equal(p['mass'],dc['mass_kg'])
 assert np.all(p['weight']>0) and np.all(np.isfinite(p['position']))
 supplied=json.loads(bind(W/'outputs/cathode-emission-audit-20261029/snapshot-distances.json').read_text())
 v2=json.loads(bind(V/'handoff.json').read_text());analysis=json.loads(bind(V/'real-seed-event-analysis.json').read_text());review=json.loads(bind(R/'handoff.json').read_text())
 assert v2['real_seed_checks']==800 and v2['corruption_cases_rejected']==35 and analysis['status']=='passed'
 assert review['status']=='passed_no_actionable_findings'
 cursor=json.loads(bind(E/'cursor.json').read_text());summary=json.loads(bind(E/'summary.json').read_text());seal=analysis['checkpoint_seal']
 data=bind(E/'events.csv',seal['csv_prefix_sha256']).read_bytes()
 assert len(data)==cursor['csv_prefix_bytes']==seal['csv_prefix_bytes'] and summary['complete'] and summary['dropped_rows']==0
 with (E/'events.csv').open(newline='') as f:rows=list(csv.DictReader(f))
 births=[r for r in rows if r['kind']=='thermal_birth'];assert len(births)==192 and len({r['particle_id'] for r in births})==192
 pos=np.array([[float(r[k]) for k in ('x_m','y_m','z_m')] for r in births]);weights=np.array([float(r['weight']) for r in births]);be=np.array([int(r['source_electrode']) for r in births]);bf=np.array([int(r['source_surface']) for r in births]);bt=np.array([float(r['time_global_s']) for r in births])
 geom={0:(.002,.007,.0018,.002),1:(.002,.007,.007,.0072)}
 for r,(x,y,z),electrode,face in zip(births,pos,be,bf):
  lo,hi,bottom,top=geom[electrode];assert lo<=x<=hi and bottom<=y<=top and -.0125<=z<=.0125
  assert face in (0,1,2,3) and r['source_known']=='1' and r['normal_convention']=='source_outward'
  assert [y==(top if electrode==0 else bottom),y==(bottom if electrode==0 else top),x==lo,x==hi][face]
  assert r['birth_event_seq']==r['event_seq'] and r['half']=='0' and float(r['flight_s'])==0
 cloud=p['position'];cw=p['weight'];slab=np.abs(cloud[:,2])<.0005
 assert slab.sum()==supplied['sampling']['numerical_1mm_slab']['macros']==1981
 def dist(rect):
  lo,hi,bottom,top=rect;dx=np.maximum(np.maximum(lo-cloud[:,0],cloud[:,0]-hi),0);dy=np.maximum(np.maximum(bottom-cloud[:,1],cloud[:,1]-top),0);return np.hypot(dx,dy)
 distance=np.minimum(dist(geom[0]),dist(geom[1]))*1e6
 def quantile(a,w,q):
  idx=np.argsort(a,kind='stable');return float(a[idx[np.searchsorted(np.cumsum(w[idx]),q*sum(w),side='left')]])
 for key,mask in (('full_height',np.ones(len(p),dtype=bool)),('numerical_1mm_slab',slab)):
  target=supplied['sampling'][key]
  assert np.isclose(sum(cw[mask]),target['physical_electrons'],rtol=2e-13,atol=1e-11)
  for q in (.05,.5,.95):assert abs(quantile(distance[mask],cw[mask],q)-target['distance_to_nearest_foil_um'][str(q)])<1e-9
 window=[(summary['global_origin_s']+summary[k]['ledger']['time_s'])*1e9 for k in ('initial','final')]
 assert np.all(bt*1e9>=window[0]) and np.all(bt*1e9<=window[1])
 phys=(cloud[:,:2]-.0045)*1e6;birth_xy=(pos[:,:2]-.0045)*1e6
 crop=[-3300,-1700,1600,3500];near=slab&(phys[:,0]>=crop[0])&(phys[:,0]<=crop[1])&(phys[:,1]>=crop[2])&(phys[:,1]<=crop[3])
 plt.rcParams.update({'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white','font.size':9,'axes.spines.top':False,'axes.spines.right':False})
 fig,axes=plt.subplots(2,2,figsize=(12,9),dpi=150)
 fig.subplots_adjust(left=.08,right=.96,bottom=.10,top=.89,wspace=.29,hspace=.38)
 colors={0:'#008d83',1:'#777777',2:'#6956a7',3:'#db7722'};names={0:'Gap face',1:'Back face',2:'Left face',3:'Right face'}
 def foils(ax):
  for e,(lo,hi,bo,to) in geom.items():ax.add_patch(Rectangle(((lo-.0045)*1e6,(bo-.0045)*1e6),(hi-lo)*1e6,(to-bo)*1e6,facecolor='#e6e6e6',edgecolor='#343434',lw=.9,zorder=1))
  ax.set_aspect('equal',adjustable='box');ax.set_xlabel('Physical x (µm)');ax.set_ylabel('Physical y (µm)');ax.grid(alpha=.12,zorder=0)
 ax=axes[0,0];foils(ax)
 for e,marker in ((1,'o'),(0,'^')):
  for f in range(4):
   m=(be==e)&(bf==f)
   if m.any():ax.scatter(birth_xy[m,0],birth_xy[m,1],s=13,c=colors[f],marker=marker,lw=.35,edgecolors='white',zorder=3)
 ax.set(xlim=(-2900,2900),ylim=(-2950,2950),title=f'A  Exact birth positions | {window[0]:.6f}–{window[1]:.6f} ns')
 ax.text(.02,.78,'192 births: upper cathode 96; lower anode 96\nAll z across the full 25 mm height\nBirth markers overlap in this x–y projection',transform=ax.transAxes,va='top',fontsize=8)
 ax.text(0,2440,'Upper cathode',ha='center',va='top',fontsize=8);ax.text(0,-2440,'Lower anode',ha='center',va='bottom',fontsize=8)
 handles=[Line2D([0],[0],marker='o',color='none',markerfacecolor=colors[f],label=names[f]+f' ({sum(bf==f)})',markersize=5) for f in (0,2,3)]
 ax.legend(handles=handles,loc='center',ncol=1,frameon=False,fontsize=8,title='Face colors; ○ upper / △ lower',title_fontsize=8)
 ax=axes[0,1];foils(ax);ax.scatter(phys[slab,0],phys[slab,1],s=3,c='#2c6e8f',alpha=.7,lw=0,zorder=3)
 ax.set(xlim=(-4500,4500),ylim=(-4500,4500),title='B  Saved cloud | t = 26.4784435 ns')
 ax.text(.02,.97,'|z| < 500 µm: all 1,981 saved macros\nUniform color; source origins are not assigned',transform=ax.transAxes,va='top',fontsize=8)
 ax.add_patch(Rectangle((crop[0],crop[2]),crop[1]-crop[0],crop[3]-crop[2],fill=False,edgecolor='#b55b27',lw=1.2,zorder=4))
 ax=axes[1,0];foils(ax);ax.scatter(phys[near,0],phys[near,1],s=9,c='#2c6e8f',alpha=.85,lw=0,zorder=3)
 ax.set(xlim=crop[:2],ylim=crop[2:],title='C  Same saved slab | upper-left cathode edge')
 ax.text(.02,.97,f'{sum(near)} macros inside this fixed crop\nFoil thickness = 200 µm',transform=ax.transAxes,va='top',fontsize=8)
 ax.text(-2110,2600,'Cathode solid',ha='center',va='center',fontsize=8)
 ax=axes[1,1]
 for mask,label,color,ls in ((np.ones(len(p),dtype=bool),'Full height','#848b94','--'),(slab,'|z| < 500 µm','#2c6e8f','-')):
  dd=distance[mask];ww=cw[mask];idx=np.argsort(dd,kind='stable');ax.step(dd[idx],100*np.cumsum(ww[idx])/sum(ww),where='post',label=label,color=color,ls=ls,lw=1.5)
 med=supplied['sampling']['numerical_1mm_slab']['distance_to_nearest_foil_um']['0.5'];rho=supplied['sampling']['numerical_1mm_slab']['local_instantaneous_gyro_radius_um']['0.5']
 ax.axvline(med,color='#2c6e8f',lw=.8,alpha=.7)
 ax.set(xlim=(0,1000),ylim=(0,100),xlabel='Distance to nearest foil rectangle (µm)',ylabel='Cumulative physical electron weight (%)',title='D  Saved cloud distances | t = 26.4784435 ns')
 ax.grid(alpha=.18);ax.legend(loc='lower right',frameon=False,fontsize=8)
 ax.text(.47,.49,f'Slab median distance: {med:.1f} µm\nLocal instantaneous ρ median: {rho:.2f} µm\nρ is not a trajectory extent.',transform=ax.transAxes,fontsize=8)
 fig.suptitle('Logged emission faces and a later particle cloud',fontsize=15,y=.971)
 fig.text(.5,.928,'Different time windows and samples; no trajectory or birth-ID connection is implied.',ha='center',fontsize=10)
 fig.text(.5,.035,'Raw x–y projections and weighted distance steps; no smoothing. Physical mapping: x,y = PIC x,y − 4500 µm.',ha='center',fontsize=8)
 png=OUT/'birth_positions_and_saved_cloud.png';fig.savefig(png);plt.close(fig)
 info={'status':'read_only_source_vs_cloud_preview','birth_record_window_ns':window,'actual_birth_time_range_ns':[float(bt.min()*1e9),float(bt.max()*1e9)],'births':192,'birth_physical_electrons':float(sum(weights)),'birth_source_counts':{f'{e}:{f}':int(sum((be==e)&(bf==f))) for e in (0,1) for f in range(4)},'source_enums':{'electrode_0':'lower anode','electrode_1':'upper cathode','surface_0':'gap','surface_1':'back','surface_2':'left','surface_3':'right'},'birth_faces_exactly_verified':True,'snapshot_time_ns':fixture['absolute_native_time_s']*1e9,'slab_macro_count':int(sum(slab)),'slab_physical_electrons':float(sum(cw[slab])),'near_cathode_crop_bounds_um':crop,'near_cathode_crop_macro_count':int(sum(near)),'slab_distance_median_um':med,'slab_instantaneous_rho_median_um':rho,'physical_coordinate_map':'x_um=(PIC x_m−0.0045)*1e6; y_um=(PIC y_m−0.0045)*1e6; z unchanged','raw_binary_decode_matches_existing_typed_cpp_decoder_exactly':True,'distance_quantiles_reproduce_input_report':True,'current_event_analysis_status':analysis['status'],'current_event_analysis_checks':v2['real_seed_checks'],'independent_event_review_status':review['status'],'origins_inferred':False,'trajectories_inferred':False,'density_estimated':False,'histogram_or_scatter_smoothing':False,'no_new_simulation':True,'native_foil_geometry_m':geom,'cdf_display_xlim_um':[0,1000],'cdf_fraction_above_display_range':{'full_height':float(sum(cw[distance>1000])/sum(cw)),'slab':float(sum(cw[slab&(distance>1000)])/sum(cw[slab]))},'wall_s':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'cpu_affinity':sorted(os.sched_getaffinity(0)),'nice':os.getpriority(os.PRIO_PROCESS,0)}
 assert info['cpu_affinity']==[7] and info['nice']>=15
 for path,r in inputs.items():assert sha(path)==r['sha256'],path
 bind(__file__)
 (OUT/'source-cloud-summary.json').write_text(json.dumps(info,indent=2)+'\n');(OUT/'input-bindings.json').write_text(json.dumps(inputs,indent=2)+'\n')
 (OUT/'handoff.txt').write_text(f"All 192 exact thermal birth records verified on source faces during {window[0]:.9f}–{window[1]:.9f} ns, 96 per electrode. Full-height birth projection is distinct from the saved |z|<0.5 mm cloud at {info['snapshot_time_ns']:.9f} ns, 1981 macros. No IDs joined and no inherited origins inferred. Actual 200 µm foil thickness shown. Slab physical-weighted median distance {med:.9f} µm; reported local instantaneous gyro-radius median {rho:.9f} µm is not an orbit/history extent. Binary particle positions/weights match the prior native typed decoder exactly. Current v2 event analysis and independent review passed. All consumed inputs unchanged. CPU7/nice15; no PIC/GPU/field solve.\n")
 (OUT/'output-manifest.json').write_text(json.dumps({x.name:{'sha256':sha(x),'bytes':x.stat().st_size} for x in OUT.iterdir() if x.is_file() and x.name!='output-manifest.json'},indent=2)+'\n')
 print(json.dumps(info),flush=True)
if __name__=='__main__':main()
