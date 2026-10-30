#!/usr/bin/env python3
import pathlib,json,csv,argparse,hashlib,math
import numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from assess import DATA,DIAG,CHANNELS,evaluate,QE,EPS0
def main():
 p=argparse.ArgumentParser();p.add_argument('--assessment-dir',type=pathlib.Path,required=True);a=p.parse_args();out=a.assessment_dir;assert not (out/'collision-and-resolution-scales.png').exists()
 data=json.loads((out/'assessment.json').read_text());gas=json.loads(DATA.read_text());channels={c['id']:c for c in gas['channels']};diag=json.loads(DIAG.read_text())
 plt.rcParams.update({'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white','font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
 fig,axes=plt.subplots(2,2,figsize=(13.8,10.4));energy=np.geomspace(.01,1e4,1400)
 selected=[('H2O:munoz2007:recommended_total','Total scattering (recommended)'),('H2O:ion_yield:H2O+','Intact-parent ion yield'),('H2O:nist_BEB_valence','Valence BEB model'),('H2O:elastic_momentum','Momentum transfer')]
 for cid,label in selected:
  known,sigma=evaluate(channels[cid],energy);axes[0,0].plot(energy[known],sigma[known],label=label,lw=1.6)
 axes[0,0].set(xscale='log',yscale='log',xlabel='Incident electron energy (eV)',ylabel='Specific observable σ (m²)',title='Supported observables; overlapping curves are not summed');axes[0,0].legend(frameon=False,fontsize=8)
 h=diag['energy_distributions']['electron'];edges=np.array(diag['energy_edges_eV']);counts=np.array([r['rounded'] for r in h['bin_counts']]);cdf=np.cumsum(counts)/h['physical_count']['rounded'];axes[0,1].step(edges[1:],cdf,where='post',lw=1.7,color='#0072b2')
 for x,label,color in [(13.5,'Parent table onset','#009e73'),(1000,'Parent table ceiling','#d55e00')]:axes[0,1].axvline(x,color=color,ls='--',lw=1.1,label=label)
 axes[0,1].set(xscale='log',xlim=(.01,1e4),ylim=(0,1.02),xlabel='Actual kinetic energy (eV)',ylabel='Cumulative represented electron fraction',title=f"Frozen {data['time_s']*1e9:.9f} ns energy distribution");axes[0,1].text(.03,.83,f"57.697% above 1 keV\n0.628% at/below 12.621 eV",transform=axes[0,1].transAxes,fontsize=9);axes[0,1].legend(frameon=False,fontsize=8,loc='lower right')
 rows=list(csv.DictReader((out/'finite-column-optical-depth.csv').open()))
 for cid,label in selected[:3]:
  points={float(r['column_m2']):float(r['partial_known_support_formal_single_normal_crossing_attenuation']) for r in rows if r['channel']==cid};xs=sorted(points);missing=100*data['collision_observables'][cid]['unsupported_population_fraction'];axes[1,0].plot(xs,[points[x] for x in xs],marker='.',lw=1.3,label=f'{label} ({missing:.2f}% unsupported weight)')
 background=5e-5*133.32236842105263/(1.380649e-23*300)*.005;axes[1,0].axvline(background,color='#666666',ls=':',lw=1.1,label='Assumed 5 mm background column');axes[1,0].axvline(4.29e18,color='#333333',ls='--',lw=1.1,label='θ₀=.33, fully released column')
 axes[1,0].set(xscale='log',yscale='log',xlabel='Hypothetical neutral column (m⁻²)',ylabel='Known-support 1 − exp(−column × σ)',title='Straight normal crossing at fixed recorded energy');axes[1,0].legend(frameon=False,fontsize=7.2,loc='upper left');axes[1,0].text(.98,.03,'Yield/model curves are formal optical scales;\nonly total scattering is a collision attenuation.',ha='right',transform=axes[1,0].transAxes,fontsize=8)
 ne=np.geomspace(3e5,1e20,400)
 for T in [1.,10.,100.]:axes[1,1].plot(ne,1e6*np.sqrt(EPS0*T/(ne*1e6*QE)),label=f'Assumed Maxwellian {T:g} eV')
 axes[1,1].axhline(25,color='#333333',ls='--',lw=1.2,label='25 μm field cell');axes[1,1].axvline(data['resolution_summary']['max_ne_cm3'],color='#888888',ls=':',lw=1.2,label='Old electron-only peak density');axes[1,1].set(xscale='log',yscale='log',xlabel='Hypothetical electron density (cm⁻³)',ylabel='Classical Maxwellian Debye length (μm)',title='Resolution scale only; no shielding imposed');axes[1,1].legend(frameon=False,fontsize=8)
 for ax in axes.flat:ax.grid(alpha=.2)
 fig.suptitle('Conditional water collisionality and resolution — static recorded electrons and fields',fontsize=14)
 fig.text(.5,.025,'Stationary ground-state H₂O test targets; unsupported weight stays in every ensemble denominator.\nColumns and Maxwellian plasma are explicit scenarios, not reconstructed gas, equilibrium, or an explanation of the halo.',ha='center',va='bottom',fontsize=9);fig.tight_layout(rect=(0,.095,1,.95));files=[]
 for ext in ['png','svg']:
  path=out/('collision-and-resolution-scales.'+ext);fig.savefig(path,dpi=165,facecolor='white');files.append(dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
 plt.close(fig);(out/'plot-receipt.json').write_text(json.dumps(dict(input_sha256=hashlib.sha256((out/'assessment.json').read_bytes()).hexdigest(),source_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),files=files,white_background=True,density_units='cm^-3',physical_evolution=False),indent=2)+'\n')
 print('White conditional scale plot complete')
if __name__=='__main__':main()
