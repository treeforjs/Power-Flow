"""Pinned source extraction; generated numerical arrays stay local outside git."""
import hashlib
import importlib.util
import json
from decimal import Decimal as D, Context, ROUND_HALF_EVEN, localcontext
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
CHEM=Path('/home/treefor/projects/runtimes/pic-water-chemistry-v3.1.1-20261031')
REF=Path('/home/treefor/projects/runtimes/water-fragment-source-inventory-20261101/references')
CTX=Context(prec=70,rounding=ROUND_HALF_EVEN,Emin=-999999,Emax=999999)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    pins=json.loads((ROOT/'INPUTS.json').read_text())
    for p,h in pins['dependencies'].items():assert sha(Path(p))==h,p
    path=Path('/home/treefor/projects/runtimes/oh-primary-supplement-v1.2-20261101/src/oh_primary.py')
    sp=importlib.util.spec_from_file_location('oh_primary',path);module=importlib.util.module_from_spec(sp)
    sys.modules['oh_primary']=module;sp.loader.exec_module(module)
    channel=module.load_channel(REF/'PSSTab264csupp_e-OH-data.xlsx','excitation_C')
    assert channel.label=='X 2\\Pi --> A 4\\Sigma-' and channel.value_unit=='bohr^2' and channel.energy_unit=='eV'
    assert len(channel.rows)==1451
    positive=[(e,s) for e,s in channel.rows if s>0]
    assert positive[0][0]==D('7.3536999999999999') and channel.rows[-1][0]==D(10)
    assert max(s for _,s in channel.rows)==D('0.47727999999999998')
    assert all(s>0 for e,s in channel.rows if e>=D('7.3605'))
    registry=json.loads((CHEM/'data/tables/chemistry-v3.1-manifest.json').read_text())
    policy=json.loads((CHEM/'data/tables/chemistry-v3.1-coverage-policy.json').read_text())
    assert policy['admitted_unbounded_regions']==[]
    byname={s['name']:s for s in registry['species']}
    selected=[byname[n] for n in ('e','OH','H','O')]
    q=byname['H']['ref_eV']+byname['O']['ref_eV']-byname['OH']['ref_eV']
    # Constants match the unchanged registry; never adjust references to source thresholds.
    a0=registry['constants']['bohr_m']
    rows=[(e,s) for e,s in channel.rows if e>=D('7.3605')]
    with localcontext(CTX):
        m,M=[D.from_float(s['mass_eV']) for s in selected[:2]]
        loss=D.from_float(7.36);threshold=loss*(2*(m+M)+loss)/(2*M)
        assert rows[0][0]>threshold
        converted=[(float(e),float(s*D.from_float(a0)**2)) for e,s in rows]
    for s in selected:assert s['charge'] in (-1,0)
    assert selected[0]['charge']+selected[1]['charge']==selected[0]['charge']+selected[2]['charge']+selected[3]['charge']
    for atom in ('H','O'):
        assert selected[1]['atoms'].get(atom,0)==selected[2]['atoms'].get(atom,0)+selected[3]['atoms'].get(atom,0)
    lines=['// GENERATED LOCAL PRIMARY DATA; do not commit.','#pragma once','namespace ohdiag {',
           'inline constexpr double loss_eV=7.36;',f'inline constexpr double reference_cost_eV={q.hex()};',
           f'inline constexpr unsigned point_count={len(rows)};', 'inline const pfc::SpeciesRec species[4]={']
    for s in selected:
        nuclei=','.join(str(s['atoms'].get(n,0)) for n in ('H','O','Fe','Cr','Ni','spare'))
        kind=0 if s['name']=='e' else 3
        lines.append('{'+f"{s['mass_eV'].hex()},{s['ref_eV'].hex()},{s['charge']},{kind},"+'{'+nuclei+'},{0,0,0,0}},')
    lines+=['};','inline const double energies[]={'+','.join(e.hex() for e,_ in converted)+'};',
            'inline const double cross_sections_m2[]={'+','.join(s.hex() for _,s in converted)+'};','}']
    (ROOT/'build/generated_channel.hpp').write_text('\n'.join(lines)+'\n')
    metadata={'schema':'oh-a4sigma-diagnostic-v1','species':selected,'source_point_count':len(channel.rows),
              'supported_point_count':len(rows),'source_window_eV':['7.3605','10'],'reference_cost_eV':q,
              'source_vertical_loss_decimal_eV':'7.36','solver_loss_binary64_eV':7.36,
              'solver_target_rest_threshold_eV':str(threshold),'bohr_radius_m':a0,
              'source_workbook_sha256':sha(REF/'PSSTab264csupp_e-OH-data.xlsx'),
              'generated_header_sha256':sha(ROOT/'build/generated_channel.hpp'),
              'source_rate_class':'calculated','event_kinematics_class':'explicit diagnostic model',
              'electron_angular_model':'isotropic in pair CM','fragment_angular_model':'isotropic two-body breakup in excited-target rest frame',
              'ground_fragment_KER_eV':7.36-q,'production_qualified':False,'physical_model_validated':False,
              'native_adopted':False,'network_complete':False,'whole_encounter_group_supported':False}
    (ROOT/'build/model.json').write_text(json.dumps(metadata,indent=2,sort_keys=True)+'\n')
if __name__=='__main__':main()
