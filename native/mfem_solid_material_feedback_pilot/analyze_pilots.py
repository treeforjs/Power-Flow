from pathlib import Path
import json,csv,hashlib,math,sys,os
ROOT=Path(sys.argv[1]);os.environ['MPLCONFIGDIR']=str(ROOT/'.matplotlib');names=['mesh32-dt250ps','mesh64-dt250ps','mesh64-dt125ps'];report={'schema':'solid-material-feedback-pilot-results-v1','full_plasma_run':False,'GPU_work':False,'cases':[]};allrows={}
for name in names:
 p=ROOT/name;summary=json.loads((p/'summary.json').read_text());rows=[{k:float(v) for k,v in r.items()}for r in csv.DictReader((p/'trajectory.csv').open())];allrows[name]=rows
 def first(T):
  return next(({'time_s':r['time_s'],'Tmax_K':r['max_T_K'],'current_A':r['current_A']}for r in rows if r['max_T_K']>=T),None)
 max_scaled_em=0;max_current=0
 for prev,r in zip(rows,rows[1:]):
  scale=abs(r['terminal_work_increment_J'])+abs(r['magnetic_energy_J']-prev['magnetic_energy_J'])+r['joule_increment_J']+abs(r['BE_magnetic_damping_J']);max_scaled_em=max(max_scaled_em,abs(r['EM_residual_J'])/scale if scale else 0)
 for r in rows:max_current=max(max_current,abs(r['upper_current_A']-r['current_A']),abs(r['lower_current_A']+r['current_A']))
 initial,last=rows[0],rows[-1];closure=last['cumulative_work_J']-(last['magnetic_energy_J']-initial['magnetic_energy_J'])-(last['thermal_energy_J']-initial['thermal_energy_J'])-last['cumulative_damping_J']
 item={'case':name,'mesh':summary['mesh'],'max_dt_s':summary['max_dt_s'],'accepted_steps':summary['accepted_steps'],'accepted_time_support_s':[initial['time_s'],last['time_s']],'last_Tmax_K':last['max_T_K'],'last_Tmin_K':last['min_T_K'],'last_current_A':last['current_A'],'last_max_phase':int(last['max_phase']),'last_max_liquid_fraction':last['max_liquid_fraction'],'remaining_total_hot_headroom_J':last['hot_headroom_J'],'minimum_local_hot_headroom_J':last['min_node_hot_headroom_J'],'field_sigma_min_last_interval_S_m':last['sigma_min_S_m'],'initial_sigma_S_m':initial['sigma_min_S_m'],'first_773_15K_mechanics_data_bound_crossing':first(773.15),'first_solidus_crossing':first(1672),'first_liquidus_crossing':first(1727),'last_energy':{'cumulative_terminal_work_J':last['cumulative_work_J'],'magnetic_change_J':last['magnetic_energy_J']-initial['magnetic_energy_J'],'thermal_change_J':last['thermal_energy_J']-initial['thermal_energy_J'],'cumulative_BE_magnetic_damping_J':last['cumulative_damping_J'],'cumulative_balance_residual_J':closure},'checks':{'max_absolute_current_reconstruction_error_A':max_current,'max_operand_scaled_EM_residual':max_scaled_em,'max_absolute_thermal_residual_J':max(abs(r['thermal_residual_J'])for r in rows),'max_absolute_Joule_mapping_error_J':max(abs(r['joule_mapping_error_J'])for r in rows),'max_native_linear_relative_residual':max(r['linear_relative_residual']for r in rows)},'stop':summary['stop'],'trajectory_sha256':hashlib.sha256((p/'trajectory.csv').read_bytes()).hexdigest()}
 report['cases'].append(item)
# Compare only actual common accepted schedule times; no interpolated field states.
common={}
for r in allrows[names[-1]]:common[r['time_s']]=r
comparisons=[]
for name in names[:-1]:
 differences=[{'time_s':r['time_s'],'Tmax_difference_K':r['max_T_K']-common[r['time_s']]['max_T_K']}for r in allrows[name]if r['time_s']in common]
 comparisons.append({'case':name,'reference':'mesh64-dt125ps','common_actual_times':len(differences),'last_common':differences[-1]if differences else None,'maximum_absolute_peak_temperature_difference_K':max((abs(x['Tmax_difference_K'])for x in differences),default=None),'convergence_claim':False})
report['comparison']=comparisons;report['restart_proof']=json.loads((ROOT/'restart-proof.json').read_text());report['limits']=['Reduced z-invariant A_z MQS; not Hall axial total-EM or full Maxwell/plasma','Fixed geometry and mass, insulating thermal boundaries, no vapor/flow/mechanics/radiation/sources','Nominal conditional q2 material/rho8000 Cp-as-Cv; not measured foil calibration','Coarse sharp-corner peaks and guard times are mesh and interval dependent; no spatial convergence claim','C19 plot digitization and unshifted clock registration; initial DC state at first supported knot does not reconstruct prehistory','Rejected bound diagnostic is an energy residual certificate, never an accepted clamped/extrapolated temperature trace'];(ROOT/'results.json').write_text(json.dumps(report,indent=2)+'\n')
try:
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 fig,ax=plt.subplots(2,2,figsize=(10,7),layout='constrained')
 for name in names:
  rows=allrows[name];t=[r['time_s']*1e9 for r in rows];label=name.replace('mesh','').replace('-dt',', ')
  ax[0,0].plot(t,[r['max_T_K']for r in rows],label=label,lw=1.5)
 ax[0,0].axhspan(1672,1727,color='.5',alpha=.15);ax[0,0].text(.02,1740,'mushy interval 1672–1727 K',transform=ax[0,0].get_yaxis_transform(),fontsize=7,va='bottom')
 ax[0,0].axhline(3000,color='.5',ls=':',lw=.8);ax[0,0].text(.02,3000,'selected material limit',transform=ax[0,0].get_yaxis_transform(),fontsize=7,va='bottom')
 ax[0,0].set(ylabel='Maximum accepted temperature (K)',xlabel='Published C19 clock (ns)');ax[0,0].legend(fontsize=8)
 for name in names:
  rows=allrows[name];ax[0,1].plot([r['time_s']*1e9 for r in rows],[r['sigma_min_S_m']/1e6 for r in rows],label=name,lw=1.3)
 ax[0,1].set(ylabel='Minimum sigma used in field solve (MS/m)',xlabel='Published C19 clock (ns)')
 rows=allrows[names[-1]];t=[r['time_s']*1e9 for r in rows]
 for key,label in [('cumulative_work_J','Terminal work'),('cumulative_joule_J','Joule to heat'),('cumulative_damping_J','BE magnetic damping')]:ax[1,0].plot(t,[r[key]for r in rows],label=label)
 ax[1,0].plot(t,[r['magnetic_energy_J']-rows[0]['magnetic_energy_J']for r in rows],label='Magnetic energy change');ax[1,0].set(ylabel='Energy since supported start (J)',xlabel='Published C19 clock (ns)');ax[1,0].legend(fontsize=8)
 for key,label in [('EM_residual_J','EM'),('thermal_residual_J','Thermal'),('joule_mapping_error_J','Joule mapping')]:ax[1,1].plot(t,[r[key]for r in rows],label=label,lw=.8)
 ax[1,1].set(ylabel='Per-interval balance residual (J)',xlabel='Published C19 clock (ns)');ax[1,1].legend(fontsize=8);ax[1,1].ticklabel_format(axis='y',style='sci',scilimits=(0,0))
 for a in ax.flat:a.grid(alpha=.2)
 fig.suptitle('Solid-only C19 material feedback: fixed geometry, no particle sources',fontsize=12)
 fig.savefig(ROOT/'material-feedback-trajectories.png',dpi=170);plt.close(fig)
 report['figure']='material-feedback-trajectories.png'
except ImportError:report['figure']=None;report['figure_note']='matplotlib unavailable in frozen runtime; numeric CSV retained; no dependency installed'
(ROOT/'results.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
