from pathlib import Path
import json,sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
run=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True)
rows=[json.loads(line) for line in (run/'benchmark.log').read_text().splitlines() if line.startswith('{')]
pilot=next(r for r in rows if r.get('family')=='axial_native_field_enthalpy')
h=pilot['history'];t=np.array([(r['time_s']-7.49788055492e-10)*1e12 for r in h])
a=lambda key:np.array([r[key] for r in h])
plt.rcParams.update({'font.size':9,'axes.facecolor':'white','figure.facecolor':'white','axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(2,2,figsize=(10,6.5),layout='constrained')
ax=axes[0,0];i0=6012.07136047
for key,label,color in [('current_A','Feed − initial current','#222222'),('region_conduction_A','Conduction − initial current','#7446a7'),('region_plasma_A','Plasma','#1982a5'),('region_displacement_A','Displacement','#d87527')]:
 ax.plot(t,a(key)-(i0 if key in ('current_A','region_conduction_A') else 0),label=label,color=color,lw=1.6)
ax.set(title='Total current partitions at the declared port',ylabel='Current components / changes (A)');ax.legend(fontsize=8)
ax=axes[0,1];ax.plot(t,a('terminal_V'),color='#225ea8',lw=2);ax.ticklabel_format(axis='y',style='plain',useOffset=False)
ax.set(title='Terminal voltage is predicted',ylabel='Lower − upper terminal voltage (V)')
ax=axes[1,0];ax.plot(t,(a('max_T_K')-300)*1e6,color='#b44e25',lw=2)
ax.set(title='Actual Joule work heats the native enthalpy owner',ylabel='Maximum temperature rise (µK)')
ax=axes[1,1];ax.plot(t,np.cumsum(a('joule_J'))*1e6,label='Field Joule debit',color='#147f6b',lw=2)
ax.plot(t,np.cumsum(a('thermal_delta_J'))*1e6,label='Thermal energy gain',color='#333333',lw=1.5,ls='--')
ax.set(title='The same Joule transfer is booked once',ylabel='Cumulative energy (µJ)');ax.legend(fontsize=8)
for ax in axes.flat:ax.set_xlabel('Time since supported start (ps)');ax.grid(alpha=.18)
fig.suptitle('Axial field / native particle / thermal pilot — 32 steps at 10 fs',fontsize=13)
fig.savefig(out/'axial-coupling.png',dpi=180,facecolor='white');fig.savefig(out/'axial-coupling.svg',facecolor='white');plt.close(fig)
svg=out/'axial-coupling.svg';svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
ref=json.loads((run/'refinement.json').read_text())
fig,ax=plt.subplots(figsize=(6.6,4.2),layout='constrained')
for key,label,color in [('E_difference_over_actual_response','Electric field','#225ea8'),('B_difference_over_actual_response','Magnetic field','#b44e25')]:
 ax.loglog([10,5],ref['time'][key],marker='o',label=label,color=color)
ax.set(xlabel='Coarser step in successive comparison (fs)',ylabel='Difference / actual field response',title='Temporal refinement of the coupled discrete system')
ax.grid(which='both',alpha=.2);ax.legend();ax.set_xticks([5,10],labels=['5','10']);ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
fig.savefig(out/'axial-time-refinement.png',dpi=180,facecolor='white');plt.close(fig)
print(out/'axial-coupling.png')
