"""Export the sourced H2 radiative graph and probe native transient branching.

This supplies a radiative subgraph, not a complete neutral CR model. The native
control uses actual curated A values and an explicitly authored initial level
population; no electron excitation or measured-spectrum fit is inferred.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.meta_path[:] = [f for f in sys.meta_path if type(f).__module__ != '_editable_skbc_mhdlab']
from mhdlab.resolved_emission_observation import ResolvedTransition, project_resolved_emission

HC = 6.62607015e-34 * 299792458.0
EV = 1.602176634e-19


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def state_id(state, vibrational, rotational):
    return f'H2:{state}:v{int(vibrational)}:J{int(rotational)}'


def radiative_graph(path):
    rows = list(csv.DictReader(path.open(newline='')))
    catalog, reactions, transitions = {}, [], []
    for index, row in enumerate(rows):
        identifiers = []
        for side in ('upper', 'lower'):
            key = state_id(row[side+'_state'], row[side+'_v'], row[side+'_j'])
            identifiers.append(key)
            catalog[key] = {'state_id':key, 'species':'H2', 'electronic_state':row[side+'_state'],
                            'v':int(row[side+'_v']), 'J':int(row[side+'_j']),
                            'quantum_basis':'MOLAT labels retained; parity/nuclear-spin mapping to collision data remains required'}
        upper, lower = identifiers
        a, wave, wavenumber = (float(row[name]) for name in
            ('einstein_a_per_s', 'vacuum_wavelength_nm', 'transition_wavenumber_per_cm'))
        if not (math.isfinite(a) and a >= 0 and math.isfinite(wave) and wave > 0 and wavenumber > 0):
            raise ValueError('invalid radiative source value')
        if not math.isclose(wave*wavenumber, 1e7, rel_tol=2e-14):
            raise ValueError('wavelength does not match retained wavenumber')
        if a == 0:
            continue
        identity = f'MOLAT:{row["source_file"]}:row{index+2}'
        energy_ev = HC/(wave*1e-9)/EV
        reactions.append({'id':identity, 'reactants':{upper:1}, 'products':{lower:1},
            'rate':{'type':'constant','A':a}, 'source':'data/cr_model/molat_h2/raw/'+row['source_file'],
            'status':'curated_radiative_subset',
            'energy_eV_per_event':{'gas/internal':-energy_ev,'gas/radiation':energy_ev},
            'external_energy_eV_per_event':0.0})
        transitions.append({'transition_id':identity,'upper_state_id':upper,'lower_state_id':lower,
                            'vacuum_wavelength_nm':wave,'einstein_a_per_s':a})
    mechanism = {'schema_version':1, 'id':'H2_MOLAT_radiative_bound_bound_subset', 'status':'research',
        'species':[{'id':key,'elements':{'H':2},'charge':0,'phase':'gas'} for key in sorted(catalog)],
        'reactions':reactions,
        'provenance':[{'source':'MOLAT Abgrall et al. 1994 retained radiative tables',
                       'csv_sha256':sha(path),'metadata_sha256':sha(path.with_name('metadata.json'))}],
        'coverage':{'included_families':['H2 B/C+- to X bound-bound spontaneous branches'],
            'missing_physics':['electron/heavy excitation and de-excitation','omitted electronic-state cascades',
                'state-dependent quenching and predissociation','ionization/recombination/dissociation channels',
                'nascent source populations and recoil','electrons/EEDF and coupled energy closure',
                'radiative transfer and opacity','physical level energies/weights and collision-basis mapping'],
            'physical_full_cr_ready':False}}
    return mechanism, list(catalog.values()), transitions


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raikou-root',type=Path,required=True)
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    out=args.output_dir.resolve()
    if out.exists(): parser.error('use a new output directory')
    out.mkdir(parents=True)
    executable=args.executable.resolve()
    revision=subprocess.check_output(['git','-C',str(args.raikou_root),'rev-parse','HEAD'],text=True).strip()
    external_status=subprocess.check_output(['git','-C',str(args.raikou_root),'status','--porcelain'],text=True)
    write(out/'source_status.json',{'source_revision_reference':revision,'git_status_porcelain':external_status,
        'external_source_modified_by_probe':False,'binary_source_revision_match_established':False,
        'executable':str(executable),'executable_sha256':sha(executable)})

    def execute(operation, request, name):
        write(out/(name+'_request.json'),request)
        proc=subprocess.run([str(executable),operation,'--input','-'],input=json.dumps(request),
                            text=True,capture_output=True,timeout=60,check=False)
        (out/(name+'_stdout.json')).write_text(proc.stdout,encoding='utf-8')
        (out/(name+'_stderr.txt')).write_text(proc.stderr,encoding='utf-8')
        if proc.returncode: raise RuntimeError(f'{operation} failed: {proc.stderr[:300]}')
        return json.loads(proc.stdout)

    capabilities=execute('capabilities',{},'capabilities')
    if 'chemistry.advance' not in capabilities.get('operations',[]):
        raise ValueError('runtime lacks chemistry.advance')
    graph,catalog,transitions=radiative_graph(ROOT/'data/cr_model/molat_h2/transitions.csv')
    write(out/'h2_radiative_mechanism.json',graph)
    write(out/'h2_state_catalog.json',{'status':'radiative_basis_only','states':catalog,'transitions':transitions})
    graph_validation=execute('mechanism.validate',graph,'h2_radiative_validate')
    if graph_validation.get('valid') is not True: raise ValueError('native radiative graph rejected')

    water=json.loads((args.raikou_root/'mechanisms/h2o_research.json').read_bytes())
    water_validation=execute('mechanism.validate',water,'h2o_research_validate')
    upper=state_id('B',0,0)
    selected=[r for r in graph['reactions'] if upper in r['reactants']]
    selected_lines=[t for t in transitions if t['upper_state_id']==upper]
    ids=sorted({key for r in selected for terms in ('reactants','products') for key in r[terms]})
    control=dict(graph,id='H2_MOLAT_single_upper_decay_control',status='manufactured',
                 species=[s for s in graph['species'] if s['id'] in ids],reactions=selected)
    control['coverage']=dict(graph['coverage'],initial_populations='Authored upper density1e12/m3; all other states empty')
    write(out/'control_mechanism.json',control)
    initial={key:1e12 if key==upper else 0. for key in ids}
    times=np.array([0,.01,.025,.05,.1,.2,.4,.6,.8,1.,1.5,2.,3.,4.,6.,10.,20.,40.,60.])*1e-9
    history=[initial]
    native_energy=0.; accepted=0
    options={'relative_tolerance':1e-7,'absolute_tolerance_m3':1e-6,'maximum_steps':100000}
    for i,dt in enumerate(np.diff(times)):
        result=execute('chemistry.advance',{'mechanism':control,'state':{'densities_m3':history[-1]},
            'environment':{'gas_temperature_K':300.,'electron_temperature_eV':0.},
            'dt_s':float(dt),'options':options,'execution':{'device':'cpu','mpi':False,'threads':1}},f'advance_{i:03}')
        if result.get('success') is not True or not math.isclose(result['advanced_time_s'],float(dt),rel_tol=1e-12):
            raise ValueError('native chemistry did not complete control interval')
        state=result['state']['densities_m3']
        if set(state)!=set(initial) or any(not math.isfinite(v) or v<0 for v in state.values()):
            raise ValueError('invalid native state readback')
        history.append(state)
        native_energy += result['integrated_energy_J_m3'].get('gas/radiation',0.)
        accepted += result['accepted_steps']
    densities={key:np.array([row[key] for row in history]) for key in ids}
    total_a=sum(r['rate']['A'] for r in selected)
    exact={upper:1e12*np.exp(-total_a*times)}
    for r in selected:
        lower=next(iter(r['products']))
        exact[lower]=1e12*r['rate']['A']/total_a*(-np.expm1(-total_a*times))
    error=max(float(np.max(np.abs(densities[key]-exact[key]))/1e12) for key in ids)
    conservation=max(abs(sum(row.values())/1e12-1) for row in history)
    expected_energy=sum(1e12*r['rate']['A']/total_a*(-math.expm1(-total_a*times[-1]))
                        *r['energy_eV_per_event']['gas/radiation']*EV for r in selected)
    energy_error=abs(native_energy/expected_energy-1)
    np.savez_compressed(out/'control_populations.npz',time_s=times,state_ids=np.array(ids),
                        number_density_m3=np.stack([densities[key] for key in ids],axis=1))
    lines=[ResolvedTransition(t['transition_id'],upper,t['vacuum_wavelength_nm'],t['einstein_a_per_s'])
           for t in selected_lines]
    spectra=[]
    for start,stop in [(0.,3e-9),(3e-9,6e-9)]:
        spectra.append(project_resolved_emission(time_s=times,upper_state_density_m3={upper:densities[upper][:,None]},
            path_weights_m=np.array([1e-3]),transitions=lines,bin_edges_nm=np.linspace(115,165,1001),
            gate_start_s=start,gate_end_s=stop,fwhm_nm=.1,time_registration='Authored control photon_time=native_time; no shot registration',
            population_provenance='Native chemistry.advance on sourced radiative subset with authored initial upper population',
            geometry_provenance='Authored single-cell1mm path; no collection-optics calibration',
            status='manufactured',optical_model='optically_thin'))
    write(out/'control_gate_spectra.json',spectra)
    plot(out,times,densities,upper,spectra)
    source_files=['include/raikou/chemistry.h','src/chemistry.cpp','include/raikou/collisional_radiative_model.h',
                  'src/collisional_radiative_model.cpp','mechanisms/h2o_research.json']
    report={'schema':'mhdlab-neutral-cr-start-v1','created_utc':datetime.now(timezone.utc).isoformat(),
        'status':'manufactured_population_control_with_sourced_radiative_subset','source_revision_reference':revision,
        'source_files_sha256':{name:sha(args.raikou_root/name) for name in source_files},
        'executable_sha256':sha(executable),'binary_source_revision_match_established':False,
        'external_repository_initially_clean':not external_status.strip(),'external_repository_modified':False,
        'full_cr_physical_ready':False,'neutral_density_or_temperature_inferred':False,
        'radiative_graph':{'species_count':len(graph['species']),'reaction_count':len(graph['reactions']),
            'upper_levels':len({t['upper_state_id'] for t in transitions}),
            'lower_levels':len({t['lower_state_id'] for t in transitions}),
            'lines_in115_165_nm':sum(115<=t['vacuum_wavelength_nm']<=165 for t in transitions),
            'native_validation':graph_validation['valid'],'complete_kinetics':False},
        'water_scaffold':{'species_count':len(water['species']),'reactions':len(water['reactions']),
            'native_validation':water_validation.get('valid'),'coverage':water_validation.get('coverage')},
        'control':{'upper_state_id':upper,'radiative_branches':len(selected),'total_A_s':total_a,
            'lifetime_s':1/total_a,'authored_initial_density_m3':1e12,'authored_path_m':1e-3,
            'maximum_population_error_over_initial':error,'maximum_molecule_balance_relative_error':conservation,
            'radiation_event_energy_relative_error':energy_error,'native_accepted_steps':accepted},
        'passed':error<1e-3 and conservation<1e-10 and energy_error<1e-3,
        'tool_sha256':sha(__file__),'input_csv_sha256':sha(ROOT/'data/cr_model/molat_h2/transitions.csv')}
    report['artifacts_sha256']={p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()}
    write(out/'report.json',report)
    print(json.dumps({'report':str(out/'report.json'),'passed':report['passed'],'control':report['control']}))
    if not report['passed']:
        raise SystemExit(1)


def plot(out,times,densities,upper,spectra):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4.2),layout='constrained')
    axes[0].plot(times*1e9,densities[upper]/1e12,label=upper)
    for name in sorted((key for key in densities if key!=upper),key=lambda key:densities[key][-1],reverse=True)[:4]:
        axes[0].plot(times*1e9,densities[name]/1e12,label=name)
    axes[0].set(xlim=(0,6),xlabel='Control time [ns]',ylabel='Fraction of initial H2 molecules',title='Native transient radiative branching')
    axes[0].legend(fontsize=7)
    for spec,label in zip(spectra,['0-3ns control gate','3-6ns control gate']):
        edges=np.array(spec['bin_edges_nm'])
        axes[1].stairs(spec['bin_mean_radiance_w_m2_sr_nm'],edges,label=label)
    axes[1].set(xlabel='Vacuum wavelength [nm]',ylabel='Control radiance [W m$^{-2}$ sr$^{-1}$ nm$^{-1}$]',title='Sourced H2 branches; authored population/path')
    axes[1].legend(fontsize=8)
    fig.suptitle('Radiative-subset control only: electron excitation and full neutral CR remain missing')
    fig.savefig(out/'neutral_cr_radiative_control.png',dpi=150)
    plt.close(fig)


if __name__=='__main__': main()
