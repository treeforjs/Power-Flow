"""Generate local pinned operands; primary numerical arrays stay outside git."""
from decimal import Decimal as D, Context, ROUND_HALF_EVEN, localcontext
import hashlib, importlib.util, json, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
CHEM=Path('/home/treefor/projects/runtimes/pic-water-chemistry-v3.1.1-20261031')
REF=Path('/home/treefor/projects/runtimes/water-fragment-source-inventory-20261101/references')
CTX=Context(prec=70,rounding=ROUND_HALF_EVEN,Emin=-999999,Emax=999999)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def module(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m
def main():
    pins=json.loads((ROOT/'INPUTS.json').read_text())
    for p,h in pins['dependencies'].items():assert sha(Path(p))==h,p
    cdf=module('primary_cdf',ROOT/'src/primary_cdf.py').load(REF/'sartori2024-table-I.txt')
    h2=module('h2_primary',Path('/home/treefor/projects/runtimes/h2-dissociation-primary-v1.1-20261101/src/h2_primary.py'))
    raw=(REF/'MCCC-el-H2-b3Su_DE.X1Sg_vi=0.txt').read_bytes();assert hashlib.sha256(raw).hexdigest()==h2.SOURCE_SHA256
    rate=h2._parse_source(raw);rows=[(e,s) for e,s in rate.rows if e>=D('4.5')]
    assert len(rows)==53 and rows[-1][0]==150 and all(s>0 for _,s in rows)
    reg=json.loads((CHEM/'data/tables/chemistry-v3.1-manifest.json').read_text());by={s['name']:s for s in reg['species']}
    selected=[by[n] for n in ('e','H2','H')];q=2*by['H']['ref_eV']-by['H2']['ref_eV'];a0=reg['constants']['bohr_m']
    assert selected[1]['atoms']=={'H':2} and selected[2]['atoms']=={'H':1}
    assert selected[0]['charge']==-1 and selected[1]['charge']==selected[2]['charge']==0
    with localcontext(CTX):
        m,M=[D.from_float(s['mass_eV']) for s in selected[:2]];Q=D.from_float(q)
        lossmin=Q+2*cdf[0][1];threshold=lossmin*(2*(m+M)+lossmin)/(2*M)
        base_threshold=Q*(2*(m+M)+Q)/(2*M);converted=[(float(e),float(s*D.from_float(a0)**2)) for e,s in rows]
    lines=['// GENERATED LOCAL PRIMARY DATA; do not commit.','#pragma once','namespace h2diag {',f'inline constexpr double reference_cost_eV={q.hex()};',f'inline constexpr unsigned point_count={len(rows)}, cdf_count={len(cdf)};', 'inline const pfc::SpeciesRec species[3]={']
    for s in selected:
        nu=','.join(str(s['atoms'].get(n,0)) for n in ('H','O','Fe','Cr','Ni','spare'));kind=0 if s['name']=='e' else 3
        lines.append('{'+f"{s['mass_eV'].hex()},{s['ref_eV'].hex()},{s['charge']},{kind},"+'{'+nu+'},{0,0,0,0}},')
    lines+=['};','inline const double energies[]={'+','.join(e.hex() for e,_ in converted)+'};','inline const double cross_sections_m2[]={'+','.join(s.hex() for _,s in converted)+'};','inline const double quantiles[]={'+','.join(float(x).hex() for x,_ in cdf)+'};','inline const double fragment_eV[]={'+','.join(float(e).hex() for _,e in cdf)+'};','}']
    (ROOT/'build/generated_channel.hpp').write_text('\n'.join(lines)+'\n')
    meta={'schema':'h2-b3sigma-source-cdf-diagnostic-v1','species':selected,'reference_cost_eV':q,'bohr_radius_m':a0,'rate_knots':len(rows),'source_rate_window_eV':['4.5','150'],'unshifted_source_onset_eV':str(rate.threshold),'registry_target_rest_threshold_eV':str(base_threshold),'cdf_table':'I','target_state':'H2(X1Sigma_g,vi=0)','product_state':'two ground H(1s) atoms','quantile_interval':['0.001','0.999'],'omitted_probability':'0.002','cdf_points':len(cdf),'fragment_energy_unit':'eV per atom in fragment-pair rest frame','source_fragment_interval_eV':[str(cdf[0][1]),str(cdf[-1][1])],'minimum_event_target_rest_threshold_eV':str(threshold),'CDF_conditioning':'renormalize only published central quantile interval intersected with exact finite-target fragment phase space; no tail extrapolation','interpolation':'linear in published quantile and per-fragment eV; MCCC integral separately linear in incident eV and sigma','angular_model':'isotropic scattered electron in collision CM; isotropic H-pair breakup in pair rest frame','source_KER_class':'primary Franck-Condon/delta approximation; not measured spectrum','rate_class':'primary MCCC integral; unchanged by conditional KER normalization','source_CDF_sha256':sha(REF/'sartori2024-table-I.txt'),'source_rate_sha256':h2.SOURCE_SHA256,'rng_layout_version':2,'partition_word':96,'ordinary_kinematic_range':[32,96],'generated_header_sha256':sha(ROOT/'build/generated_channel.hpp'),'production_qualified':False,'physical_model_validated':False,'network_complete':False,'native_adopted':False,'whole_channel_supported':False}
    (ROOT/'build/model.json').write_text(json.dumps(meta,indent=2,sort_keys=True)+'\n')
if __name__=='__main__':main()
