"""Pinned primary-source thermochemical audit and proposed windows, no admission."""
from decimal import Decimal as D, Context, localcontext, ROUND_HALF_EVEN, InvalidOperation, DivisionByZero, Overflow
from html.parser import HTMLParser
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
CHEM=Path('/home/treefor/projects/runtimes/pic-water-chemistry-v3.2.2-auditprep-20261101')
REF=Path('/home/treefor/projects/runtimes/water-fragment-source-inventory-20261101/references')
CTX=Context(prec=70,rounding=ROUND_HALF_EVEN,Emin=-999999,Emax=999999,capitals=1,clamp=0,flags=[],traps=[InvalidOperation,DivisionByZero,Overflow])

def load_module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m)
    return m

class PrimaryRows(HTMLParser):
    def __init__(self):super().__init__();self.row=None;self.cell=None;self.rows={}
    def handle_starttag(self,t,attrs):
        if t=='tr':self.row=[]
        if t in ('td','th') and self.row is not None:self.cell=[]
    def handle_data(self,text):
        if self.cell is not None:self.cell.append(text)
    def handle_endtag(self,t):
        if t in ('td','th') and self.cell is not None:
            self.row.append(''.join(self.cell).strip());self.cell=None
        if t=='tr' and self.row is not None:
            if len(self.row)==9 and '*' in self.row[-1]:
                assert self.row[-1] not in self.rows,'duplicate primary species identity'
                self.rows[self.row[-1]]=self.row
            self.row=None

def threshold(q,m,M):
    with localcontext(CTX):
        return q*(2*(m+M)+q)/(2*M)

def available_cm(T,m,M):
    with localcontext(CTX):
        return 2*M*T/(((m+M)**2+2*M*T).sqrt()+m+M)

def main():
    pins=json.loads((ROOT/'INPUTS.json').read_text())
    for path,sha in pins['dependencies'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==sha,'input changed: '+path
    assert b'15867.862' in (REF/'nist-oxygen-levels.html').read_bytes()
    manifest=json.loads((CHEM/'data/tables/chemistry-v3.2-manifest.json').read_text())
    policy=json.loads((CHEM/'data/tables/chemistry-v3.2-coverage-policy.json').read_text())
    assert policy['admitted_unbounded_regions']==[]
    spec=json.loads((CHEM/'data/tables/chemistry-v3.2-spec.json').read_text())
    selected_path=CHEM/'data/sources/thermo-atct/atct-selected.json'
    assert hashlib.sha256(selected_path.read_bytes()).hexdigest()==manifest['input_files_sha256']['data/sources/thermo-atct/atct-selected.json']
    selected=json.loads(selected_path.read_text());primary=PrimaryRows()
    primary.feed((CHEM/'data/sources/thermo-atct/atct-main.html').read_text())
    species={s['name']:s for s in manifest['species']}
    native=json.loads(Path('/home/treefor/projects/runtimes/pic-water-chemistry-v3.1.1-20261031/data/tables/chemistry-v3.1-manifest.json').read_text())
    native_species={x['name']:x for x in native['species']}
    for name in ('e','H','H2','OH','OH-','O','O-'):
        for field in ('mass_eV','ref_eV','charge','atoms'):
            assert species[name][field]==native_species[name][field], 'native/prepared registry mismatch'
    factor=1.0/96.48533212331001  # Exact same binary64 operation as pinned builder.
    verified={}
    for s in spec['species']:
        if s['name'] not in ['H','H2','OH','OH-','O','O-']:continue
        row=selected[s['atct']];full=primary.rows[row['atct_id']]
        assert D(full[3])==D(row['DHf0_kJ_mol']) and D(full[4])==D(row['DHf298_kJ_mol'])
        assert species[s['name']]['ref_eV']==float(row['DHf0_kJ_mol'])*factor
        verified[s['name']]={'ATcT_id':row['atct_id'],'DHf0_kJ_mol':full[3],
                              'DHf298_kJ_mol':full[4],'registry_ref_eV':species[s['name']]['ref_eV']}
    singlet=primary.rows['17778-80-2*2'];triplet=primary.rows['17778-80-2*1']
    assert 'singlet' in singlet[1] and 'triplet' in triplet[1]
    assert D(triplet[3])==D(selected['O (g)']['DHf0_kJ_mol'])
    assert not any(s['name'] in ('O(1D)','O_1D') for s in manifest['species'])
    # Use the same per-species binary64 reference convention as the registry.
    hypothetical_1D_ref=float(singlet[3])*factor
    q_h2=2*species['H']['ref_eV']-species['H2']['ref_eV']
    q_oh=species['O']['ref_eV']+species['H']['ref_eV']-species['OH']['ref_eV']
    q_dea=species['O-']['ref_eV']+species['H']['ref_eV']-species['OH']['ref_eV']
    q_oh1d=hypothetical_1D_ref+species['H']['ref_eV']-species['OH']['ref_eV']
    h2m=load_module('h2_primary',Path('/home/treefor/projects/runtimes/h2-dissociation-primary-v1.1-20261101/src/h2_primary.py'))
    ohm=load_module('oh_primary',Path('/home/treefor/projects/runtimes/oh-primary-supplement-v1.2-20261101/src/oh_primary.py'))
    h2=h2m.load_primary(REF/'MCCC-el-H2-b3Su_DE.X1Sg_vi=0.txt',CHEM)
    oh=ohm.inspect(REF/'PSSTab264csupp_e-OH-data.xlsx')
    me=D.from_float(species['e']['mass_eV']);mh2=D.from_float(species['H2']['mass_eV']);moh=D.from_float(species['OH']['mass_eV'])
    windows={}
    for name,q,M,rows in [('H2_b3Su_to_2H',q_h2,mh2,h2._data.rows),
                          ('OH_unlabelled_to_H_O_ground_hypothesis',q_oh,moh,oh['dissociation'].rows),
                          ('OH_unlabelled_to_H_O_1D_hypothesis',q_oh1d,moh,oh['dissociation'].rows),
                          ('OH_DEA_to_H_Ominus',q_dea,moh,oh['DEA_B'].rows)]:
        Q=D.from_float(q);T0=threshold(Q,me,M)
        lower=next(e for e,v in rows if v>0 and e>=T0)
        assert available_cm(lower,me,M)>Q
        assert abs(available_cm(T0,me,M)-Q)<D('1e-55')
        windows[name]={'reference_cost_eV':q,'target_rest_incident_threshold_eV':str(T0),
                       'first_positive_source_knot_above_threshold_eV':str(lower),
                       'proposed_unshifted_source_window_eV':[str(lower),str(rows[-1][0])],
                       'source_knots_below_threshold_with_positive_sigma':sum(e<T0 and v>0 for e,v in rows),
                       'kinematic_window_only':True,'event_admission_qualified':False,
                       'excluded_interval_treatment':'Explicit refusal; no extrapolation or invented zero',
                       'rate_or_energy_axis_changed':False}
    with localcontext(CTX):
        delta_atct=(D(singlet[3])-D(triplet[3]))/D('96.48533212331001')
        delta_nist=D('15867.862')*D('6.62607015e-34')*D('299792458')*100/D('1.602176634e-19')
        assert abs(delta_nist-delta_atct)<D('0.00002'),'independent O(1D) energy mismatch'
        room_q=(D(selected['H (g)']['DHf298_kJ_mol'])+D(selected['O (g)']['DHf298_kJ_mol'])-D(selected['OH (g)']['DHf298_kJ_mol']))/D('96.48533212331001')
        ea_o=D.from_float(species['O']['ref_eV'])-D.from_float(species['O-']['ref_eV'])
        # Column correspondence evidence is only an arithmetic comparison.
        c=dict(oh['excitation_C'].rows);d=dict(oh['excitation_D'].rows)
        residual=[(e,v-c[e]-d[e]) for e,v in oh['dissociation'].rows if e in c and e in d]
        mismatch=[(e,v) for e,v in residual if abs(v)>D('1e-12')]
    with localcontext(CTX):
        vertical={}
        for key,loss in [('excitation_C',D('7.36')),('excitation_D',D('7.89'))]:
            t0=threshold(loss,me,moh)
            rows=oh[key].rows
            lo=next(e for e,v in rows if v>0 and e>=t0)
            entry={'source_column':key,'paper_vertical_loss_eV':str(loss),
                   'target_rest_incident_threshold_eV':str(t0),
                   'unshifted_proposed_window_eV':[str(lo),str(rows[-1][0])],
                   'ground_fragment_relative_KER_eV':str(loss-D.from_float(q_oh)),
                   'kinetic_fraction_for_ground_atomic_products':1,
                   'status':'proposal only; source rate is calculated, event kinematics model-derived',
                   'electron_angular_distribution':'not supplied; any chosen law must remain explicit and unqualified',
                   'production_admitted':False}
            if key=='excitation_D':
                index=1350-3
                bad_e,bad_v=rows[index]
                assert bad_e==D('9.2993000000000006') and bad_v==D(1974)
                entry['source_anomaly']={'cell':'Excitation!D1350','energy_cell':'A1350',
                                         'energy_eV':str(bad_e),'literal_value_bohr2':str(bad_v),
                                         'value_changed':False,'proposed_imputation':None}
                entry['unshifted_proposed_window_eV']=None
                entry['candidate_separate_complete_segments_eV']=[[str(lo),str(rows[index-1][0])],[str(rows[index+1][0]),str(rows[-1][0])]]
                entry['gap_handling']='Refuse entire open interval between the bounding adjacent knots; never interpolate across anomaly'
            vertical[key]=entry
        maxima={key:{'max_value':str(max(v for _,v in channel.rows)),
                      'max_value_abscissa':str(max(channel.rows,key=lambda row:row[1])[0]),
                      'value_unit':channel.value_unit}
                for key,channel in oh.items() if isinstance(channel,ohm.Channel)}
    result={'schema':'water-fragment-reference-and-window-audit-v1.1' ,'verified_registry_references':verified,
            'native_v3_1_1_vs_prepared_v3_2_2_species_fields_identical':['e','H','H2','OH','OH-','O','O-'],
            'proposed_vertical_breakup_models':vertical,'OH_source_column_maxima':maxima,
            'reference_convention':'Pinned ATcT v1.202 0 K formation enthalpies; electron reference zero. Registry uses stable per-species reference differences, not rest-mass subtraction.',
            'O_1D':{'already_in_registry':False,'available_primary_ATcT_id':'17778-80-2*2',
                    'available_DHf0_kJ_mol':singlet[3],'hypothetical_ref_eV':hypothetical_1D_ref,
                    'excitation_from_ATcT_eV':str(delta_atct),'NIST_level_wavenumber_cm_minus1':'15867.862',
                    'excitation_from_NIST_eV':str(delta_nist),'registry_modified':False},
            'OH_reference_comparison':{'ground_dissociation_cost_0K_eV':q_oh,'ground_dissociation_cost_298K_eV':str(room_q),
                                       'source_paper_D0_eV':'3.88','oxygen_atomic_electron_affinity_from_registry_eV':str(ea_o),
                                       'source_DEA_model_OH_molecular_electron_affinity_eV':'1.56',
                                       'note':'OH molecular electron affinity and O atomic electron affinity are distinct operands. A 0 K/298 K swap does not explain the paper D0 difference.'},
            'windows':windows,
            'OH_dissociation_column_check':{'shared_nodes':len(residual),'nodes_differing_from_excitation_C_plus_D_by_more_than_1e_minus12_bohr2':len(mismatch),
                                           'first_different_energy_eV':str(mismatch[0][0]) if mismatch else None,
                                           'max_absolute_difference_bohr2':str(max(abs(v) for _,v in residual)),
                                           'branch_identity_proven':False,'note':'Comparison is not a primary product-state label and cannot authorize a branching assignment.'},
            'controls':['selected registry refs equal retained primary 0 K cells','existing unbounded admissions empty',
                        'actual table source hashes verified by qualified loaders','all window endpoints have positive invariant excess energy',
                        'inverse invariant threshold reproduces reference cost to 1e-55 eV','independent ATcT/NIST O(1D) energy agreement',
                        'native and prepared water references and masses identical','actual D1350 anomaly located without modification'],
            'production_qualified':False,'physical_model_validated':False,'network_complete':False,'native_adopted':False}
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=='__main__':main()
