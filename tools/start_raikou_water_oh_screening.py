"""Sourced single-channel OH-yield screening through Raikou; not a complete water mechanism."""
from pathlib import Path
from decimal import Decimal
import argparse,json,hashlib,csv,math,subprocess
import numpy as np,pymupdf
ROOT=Path('/home/treefor/Projects/Power-Flow');out=Path('outputs/raikou-water-oh-screening');out.mkdir(exist_ok=False)
pdf=Path('work/song2021-water-cross-sections.pdf');digest=hashlib.sha256(pdf.read_bytes()).hexdigest();assert digest=='d94ae0c73ab745a3e5b958d096f2c774b123b8257877e2d71911990e9b83b595'
page=pymupdf.open(pdf)[17];words=page.get_text('words');manual=[('10','.103'),('15','.330'),('20','.481'),('30','.893'),('50','1.30'),('75','1.44'),('100','1.41'),('150','1.36'),('200','1.20'),('250','1.09'),('300','.96')];rows=[]
energy_words=sorted([w for w in words if 43<w[0]<60 and 150<w[1]<270],key=lambda w:w[1]);assert len(energy_words)==11
for (e,v),word in zip(manual,energy_words):
 other=[w for w in words if 225<w[0]<250 and abs(w[1]-word[1])<1];assert len(other)==1
 assert Decimal(word[4])==Decimal(e) and Decimal(other[0][4])==Decimal(v)
 rows.append(dict(energy_eV=float(e),cross_section_m2=float(Decimal(v)*Decimal('1e-20'))))
with (out/'song2021_table10_OH_yield.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
source=dict(doi='10.1063/5.0035315',source_pdf_sha256=digest,pdf_page=18,article_page='023103-17',table=10,verified_numeric_rows=11,cross_section_units='10^-16 cm^2 converted by1e-20 to m^2',table_uncertainty_relative=.4,text_uncertainty_relative=.36,uncertainty_policy='Retain40percent caption and note36percent text conflict; factors0.6/1/1.4 are sensitivity cases, not confidence intervals.',energy_support_eV=[10,300],zero_filled_or_extrapolated=False,physical_scope='Measured neutral OH yield; Eq3 suggests OH+H+e as most probable channel. OH internal states and recoil/energy partitions unresolved. Do not add an unpartitioned OH(A) production channel to this inclusive yield.')
(out/'source-data.json').write_text(json.dumps(source,indent=2)+'\n')
mechanism=dict(schema_version=1,id='H2O_Song2021_inclusive_OH_yield_screening',status='research',species=[dict(id='e_g',elements={},charge=-1,phase='gas'),dict(id='H2O',elements={'H':2,'O':1},charge=0,phase='gas'),dict(id='OH',elements={'H':1,'O':1},charge=0,phase='gas'),dict(id='H',elements={'H':1},charge=0,phase='gas')],reactions=[dict(id='H2O:inclusive_OH_yield',reactants={'e_g':1,'H2O':1},products={'e_g':1,'OH':1,'H':1},rate={'type':'electron_eedf','key':'H2O:inclusive_OH_yield'},source='https://doi.org/10.1063/5.0035315; Table10 and Eq3',status='research')],provenance=[source],coverage={'included_families':['Conditional inclusive neutral OH-production screening'],'missing_physics':['All other reaction channels','Resolved OH/H internal states and fragment recoil','Reaction energy partitions','Self-consistent electron-energy/field evolution','Surface desorption injection and gas transport coupling'],'validity':'Single-channel prescribed monoenergetic electron bath, exact table energies only. Not a complete water fragmentation prediction.'})
(out/'mechanism.json').write_text(json.dumps(mechanism,indent=2)+'\n')
exe=Path('/home/treefor/Projects/raikou/install/bin/raikou-native')
def native(op,req,name):
 (out/(name+'-request.json')).write_text(json.dumps(req,indent=2)+'\n');p=subprocess.run([str(exe),op,'--input','-'],input=json.dumps(req),text=True,capture_output=True,timeout=120);(out/(name+'-stdout.json')).write_text(p.stdout);(out/(name+'-stderr.txt')).write_text(p.stderr)
 if p.returncode:raise RuntimeError(p.stderr)
 return json.loads(p.stdout)
validation=native('mechanism.validate',mechanism,'validate');assert validation['valid']
results=[];worst=0
for entry in rows:
 E=entry['energy_eV'];nominal_k=entry['cross_section_m2']*math.sqrt(2*1.602176634e-19*E/9.1093837015e-31)
 for ne in [1e16,1e18,1e20]:
  for scale in [.6,1.,1.4]:
   k=nominal_k*scale;tag=f'E{E:g}-ne{ne:g}-scale{scale:g}'
   request=dict(mechanism=mechanism,state={'densities_m3':{'e_g':ne,'H2O':1e20,'OH':0.,'H':0.}},environment={'gas_temperature_K':300.,'electron_temperature_eV':2*E/3,'electron_rates':{'H2O:inclusive_OH_yield':k}},dt_s=60e-9,execution={'threads':1,'device':'cpu'},options={'relative_tolerance':1e-8,'absolute_tolerance_m3':1.,'maximum_step_s':1e-9,'maximum_steps':100000})
   if E==20 and ne==1e18 and scale==1:
    ev=native('chemistry.evaluate',request,'bound-evaluate');assert ev['coverage']['rates_ready'];assert ev['coverage']['unchecked_energy_reactions']==['H2O:inclusive_OH_yield']
   r=native('chemistry.advance',request,tag);assert r['success'] and r['advanced_time_s']==60e-9,r
   n=r['state']['densities_m3'];fraction=n['OH']/1e20;exact=-math.expm1(-ne*k*60e-9);error=abs(fraction-exact);worst=max(worst,error)
   hb=abs(2*n['H2O']+n['OH']+n['H']-2e20)/2e20;ob=abs(n['H2O']+n['OH']-1e20)/1e20;charge=abs(n['e_g']-ne)/ne
   assert max(hb,ob,charge)<1e-10 and error<2e-4
   results.append(dict(energy_eV=E,prescribed_electron_density_cm3=ne/1e6,cross_section_factor=scale,rate_coefficient_m3_s=k,accepted_time_s=r['advanced_time_s'],native_fragmented_parent_fraction=fraction,analytic_single_channel_fraction=exact,absolute_fraction_error=error,H_nuclei_relative_defect=hb,O_nuclei_relative_defect=ob,electron_number_relative_defect=charge,accepted_steps=r['accepted_steps']))
 print('finished E',E,flush=True)
with (out/'conditional_dissociation_scan.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=results[0]);w.writeheader();w.writerows(results)
summary=dict(status='Raikou_single_channel_conditional_screening_complete_not_full_fragmentation',native_cases=len(results),accepted_time_per_case_s=60e-9,max_absolute_fraction_error=worst,native_executable=str(exe),native_executable_sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),source_data=source,prescribed_electron_inputs_not_measured=True,initial_parent_density_cm3=1e14,initial_parent_density_scope='Authored normalization for a linear one-channel response; not inferred from RGA or assigned a physical monolayer count.',energy_ready=False,full_water_mechanism_ready=False,remaining=['Bind complete measured electron rates including neutral dissociation/excitation/attachment without double counting','State-resolved branching and fragment-energy/recoil closure','Physical monolayer definition and surface-specific desorption parameters','Electron density and EEDF/field histories','Collisional transport and gas/surface feedback'],parent_source_equations={'areal_H2O':'initial_h2o_monolayers * sigma_ML','total_H2O_for_current_four_faces':'1.04e-4 m^2 * initial_h2o_monolayers * sigma_ML','H_nuclei':'2 * total_initial_H2O','O_nuclei':'total_initial_H2O','sigma_ML_molecules_m2':None,'initial_h2o_monolayers':'adjustable; no best fit fixed'})
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
