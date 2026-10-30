from pathlib import Path
import json,sys,os
os.environ['MPLCONFIGDIR']='/tmp/hall-matplotlib'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
run=Path(sys.argv[1]);rows=json.loads((run/'results.json').read_text())
plt.rcParams.update({'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white','axes.spines.top':False,'axes.spines.right':False,'font.size':10,'axes.grid':True,'grid.alpha':.2})
native=[r for r in rows if r['family']=='native_wave']
fine=next(r for r in native if r['ny']==64 and r['steps']==1600)
hist=json.loads((run/(fine['name']+'.json')).read_text());t=np.array([r['t_s'] for r in hist])*1e12
fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained')
for s,label in enumerate(['electron','ion × 25']):
 scale=1e9*(1 if s==0 else 25)
 ax[0,0].plot(t,[r['xi_m'][s]*scale for r in hist],label=label)
 ax[0,0].plot(t,[r['exact_xi_m'][s]*scale for r in hist],color='black',ls='--',lw=.7)
ax[0,0].set(xlabel='Time (ps)',ylabel='Displacement amplitude (nm)',title='Self-consistent two-species response')
ax[0,0].legend();ax[0,0].text(.03,.03,'Dashed: linear cold-fluid solution',transform=ax[0,0].transAxes,fontsize=9)
q=1.602176634e-19;n=fine['density_cm3']*1e6
for s,label in enumerate(['electron','ion']):
 ax[0,1].plot(t,[(1 if s else -1)*q*n*r['vz_m_s'][s] for r in hist],label=label)
ax[0,1].set(xlabel='Time (ps)',ylabel='Signed Jz amplitude (A/m²)',title='Transverse current from the native Lorentz push');ax[0,1].legend()
mesh=sorted([r for r in native if r['mass_ratio']==25 and r['ppc_y']==4 and r['steps']==400],key=lambda x:x['ny'])
for key,label in [('initial_phi_relative_L2_error','Poisson potential'),('max_displacement_error_scaled','ion response')]:
 vals=[r[key][1] if isinstance(r[key],list) else r[key] for r in mesh]
 ax[1,0].loglog([r['ny'] for r in mesh],vals,'o-',label=label)
ax[1,0].set(xlabel='Cells across the wave',ylabel='Normalized error',title='Mesh refinement');ax[1,0].legend()
time=sorted([r for r in native if r['ny']==64 and r['ppc_y']==4],key=lambda x:x['steps'])
ax[1,1].loglog([r['dt_s']*1e15 for r in time],[r['max_energy_relative_defect'] for r in time],'o-',label='max |Δ(K+U)| / U₀')
ax[1,1].set(xlabel='Timestep (fs)',ylabel='Relative energy defect',title='Energy error decreases as Δt²');ax[1,1].legend()
fig.suptitle(f'Manufactured Poisson control • nₑ = nᵢ = {fine["density_cm3"]:.3e} cm⁻³ • Bx = {fine["B_T"]:.3f} T • mᵢ/mₑ = 25',fontsize=12)
fig.savefig(run/'species-poisson.png',dpi=170);plt.close(fig)
em=[r for r in rows if r['family']=='em_wave'];fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained')
for current,faraday,label in [(True,True,'coupled'),(False,True,'current feedback disabled'),(True,False,'Faraday update disabled')]:
 r=next(r for r in em if r['N']==64 and r['dt_omega_pe']==.01 and r['helicity']==1 and r['current_feedback']==current and r['faraday_feedback']==faraday)
 data=json.loads((run/(r['name']+'.json')).read_text());ts=[a['t_omega_pe'] for a in data]
 ax[0,0].plot(ts,[a['E_real']/1e-6 for a in data],label=label,lw=1)
 ax[0,1].plot(ts,[a['total_energy_normalized']-1 for a in data],label=label,lw=1)
ax[0,0].set(xlabel='Time (ωpe t)',ylabel='Re(E) / E₀',title='Electromagnetic feedback controls the mode');ax[0,0].legend(fontsize=8)
ax[0,1].set(xlabel='Time (ωpe t)',ylabel='(U − U₀) / U₀',title='One owner for field and species energy')
mesh=sorted([r for r in em if r['current_feedback'] and r['faraday_feedback'] and r['dt_omega_pe']==.01 and r['helicity']==1],key=lambda x:x['N'])
ax[1,0].loglog([r['N'] for r in mesh],[r['dispersion_relative_error'] for r in mesh],'o-')
ax[1,0].set(xlabel='Samples per periodic wavelength',ylabel='Frequency relative error',title='Staggered-curl Fourier mode converges as h²')
limits=next(r for r in rows if r['family']=='hall_limit')['rows']
ax[1,1].loglog([r['k_de'] for r in limits],[r['relative_difference'] for r in limits],'o-')
ax[1,1].set(xlabel='k dₑ',ylabel='|ωHall − ωMaxwell| / ωMaxwell',title='Approach to the Hall-MHD limit')
fig.suptitle('Manufactured cold two-fluid Maxwell mode • physical proton/electron mass ratio • no collisions or chemistry',fontsize=12)
fig.savefig(run/'electromagnetic-feedback.png',dpi=170);plt.close(fig)
print('Wrote two white-background figures')
