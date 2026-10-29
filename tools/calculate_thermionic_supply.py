"""Evaluate unrestricted thermionic supply from explicit parameters and saved surface temperatures."""
from pathlib import Path
import argparse,csv,hashlib,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from mhdlab.constants import KB,QE
from mhdlab.thermionic_supply import IDEAL_RICHARDSON_A_M2_K2 as A0,thermionic_supply
from mhdlab.mfem_surface import load_mfem_surface_samples


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-run',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args();out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=False)
    phis=[4.,4.3,4.5,5.];temperatures=np.unique(np.r_[np.arange(300.,15001.,25.),[1000.,1500.,1700.,2000.,2500.,3000.]])
    rows=[]
    for phi in phis:
        for factor in [.1,1.]:
            current=thermionic_supply(temperatures,work_function_ev=phi,richardson_a_m2_k2=A0*factor)
            for T,j in zip(temperatures,current):
                rows.append(dict(surface_temperature_K=float(T),work_function_eV=phi,Richardson_factor=factor,Richardson_A_m2_K2=A0*factor,supply_current_density_A_m2=float(j),supply_current_density_A_cm2=float(j/1e4),outgoing_number_flux_m2_s=float(j/QE),mean_flux_weighted_initial_kinetic_energy_eV=float(2*KB*T/QE),space_charge_restriction_applied=False))
    with (out/'temperature_supply_scan.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    arrays={};face_summaries=[];inputs={}
    for face in ['top','bottom','left','right']:
        path=args.source_run/f'{face}_surface.json';surface,extraction=load_mfem_surface_samples(path);T=surface.temperature_k
        current=thermionic_supply(T,work_function_ev=4.3,richardson_a_m2_k2=A0)
        area=np.diff(surface.z_edges_m)[:,None]*np.diff(surface.x_edges_m)[None,:]
        arrays[face+'_temperature_K']=T;arrays[face+'_supply_J_A_m2']=current;arrays[face+'_number_flux_m2_s']=current/QE;arrays[face+'_time_s']=surface.time_s
        for i,t in enumerate(surface.time_s):
            face_summaries.append(dict(face=face,time_ns=float(t*1e9),area_m2=float(area.sum()),temperature_min_K=float(T[i].min()),temperature_max_K=float(T[i].max()),temperature_area_mean_K=float(np.sum(T[i]*area)/area.sum()),formal_supply_current_A=float(np.sum(current[i]*area)),formal_supply_J_area_mean_A_m2=float(np.sum(current[i]*area)/area.sum()),formal_supply_J_max_A_m2=float(current[i].max()),work_function_eV=4.3,effective_Richardson_A_m2_K2=A0,physical_collected_current_established=False))
        inputs[str(path.resolve())]=hashlib.sha256(path.read_bytes()).hexdigest()
        inputs[str((args.source_run/f'{face}.csv').resolve())]=hashlib.sha256((args.source_run/f'{face}.csv').read_bytes()).hexdigest()
    with (out/'surface_supply_history.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=face_summaries[0]);w.writeheader();w.writerows(face_summaries)
    np.savez_compressed(out/'surface_supply_arrays.npz',**arrays)
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(10,6));plot_T=np.linspace(1000,4000,500)
    for phi,color in zip(phis,['#2c7bb6','#b2182b','#6a51a3','#238b45']):
        j=thermionic_supply(plot_T,work_function_ev=phi,richardson_a_m2_k2=A0)/1e4
        ax.semilogy(plot_T,j,label=f'φ = {phi:g} eV',color=color,lw=2.5 if phi==4.3 else 1.6)
        if phi==4.3:ax.fill_between(plot_T,.1*j,j,color=color,alpha=.12,label='4.3 eV: A* = 0.1–1 × ideal A₀')
    ax.set(xlabel='Surface temperature (K)',ylabel='Unrestricted thermionic current density (A/cm²)',title='Electron emission supply versus temperature | before space-charge restriction')
    ax.grid(alpha=.23);ax.legend(loc='lower right');fig.text(.11,.037,'J = A* T² exp[−φ/(kBT)]. Curves use prescribed constant work functions; no field enhancement or space-charge cap.',fontsize=9);fig.text(.11,.012,'4.3 eV: well-outgassed SS304 reference. Effective A*, contaminated-surface properties and hot-phase behavior remain uncalibrated.',fontsize=9);fig.tight_layout(rect=(0,.075,1,1));fig.savefig(out/'temperature_supply_curve.png',dpi=160)
    fig,ax=plt.subplots(figsize=(10,5.7))
    for face in ['top','bottom','left','right']:
        s=[r for r in face_summaries if r['face']==face];ax.semilogy([r['time_ns'] for r in s],[r['formal_supply_J_area_mean_A_m2']/1e4 for r in s],label=face)
    ax.set(xlabel='Simulation time (ns)',ylabel='Formal face-average thermionic supply (A/cm²)',title='Saved MFEM boundary temperatures → formal unrestricted electron supply');ax.grid(alpha=.2);ax.legend(title='Faces of the same upper foil');fig.text(.11,.035,'Benchmark only: φ = 4.3 eV, A* = ideal A₀. Evaluate J per exported temperature patch, then area-average.',fontsize=9);fig.text(.11,.012,'Hot-corner extrapolations are not physical collected currents: no phase change, field enhancement or space-charge transport is solved.',fontsize=9);fig.tight_layout(rect=(0,.075,1,1));fig.savefig(out/'formal_surface_supply_history.png',dpi=160)
    summary=dict(status='unrestricted_thermionic_supply_parameter_study',formula='J=A_eff*T^2*exp(-phi_eV*e/(k_B*T))',ideal_Richardson_A_m2_K2=A0,work_function_scan_eV=phis,Richardson_factors=[.1,1.],benchmark_work_function_eV=4.3,space_charge_cap_applied=False,field_enhancement_applied=False,actual_foil_work_function_or_effective_Richardson_calibrated=False,physical_collected_current_established=False,source_temperature='Exported material-boundary patch-average absolute temperature, not interior-cell temperature or temperature increment',nonlinear_averaging_scope='J evaluated at each exported patch-average T. Not quadrature-resolved integration of J(T) inside a patch.',formula_reference='https://goldbook.iupac.org/terms/view/09297',work_function_reference=dict(doi='10.1063/1.1708797',value_eV=4.3,conditions='well-outgassed polycrystalline304 stainless steel, clean1e-9Torr vacuum',access='publisher-deposited abstract through Crossref; full article not retrieved',applicability='Benchmark, not a measured work function of the present contaminated foil'),incoming_distribution_scope='A half-range Maxwellian thermal-emission boundary would have flux-weighted mean kinetic energy2kBT; this is not a bulk electron temperature or solved gap EEDF.',limitations=['Temperature-dependent work-function/oxide/adsorbate evolution unmodeled','No field emission or Schottky barrier lowering','No returning particles, surface charge, sheath, circuit or net collected current','Heating source uses provisional constant properties and no melting/vaporization; hot-corner emission is formal extrapolation','No photon or ion-induced secondary emission included'],source_sha256=inputs,tool_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    for r in rows:
        if r['work_function_eV']==4.3 and r['Richardson_factor']==1 and r['surface_temperature_K'] in [1000,1500,1700,2000,2500,3000]:print(json.dumps(r))
if __name__=='__main__':main()
