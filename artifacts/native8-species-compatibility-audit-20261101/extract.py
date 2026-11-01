#!/usr/bin/python -I
"""Read sealed species operands; no simulation or shared-source mutation."""
from pathlib import Path
import hashlib,json,struct
ROOT=Path('/home/treefor/projects/runtimes')
OUT=Path(__file__).resolve().parent
C=ROOT/'pic-water-common-commit-v8.0.2-cpu-20261031'
N=ROOT/'pic-reaction-transactions-20261029/native-full-time-v8.3-admitted-export-20261101'
S=ROOT/'pic-electrode-heat-source-v7-20261031'
T=ROOT/'ss304-phase-thermal-owner-v1.3-20261031'
CHEM=ROOT/'pic-water-chemistry-v3.1.1-20261031'
files=[CHEM/'data/tables/chemistry-v3.1.blob',CHEM/'data/tables/chemistry-v3.1-manifest.json',
       CHEM/'tools/build_tables.py',CHEM/'tools/make_production_spec.py',
       C/'source/identity_v8.hpp',C/'source/genesis_v8.hpp',C/'source/joint_step_v8.hpp',
       N/'source/src/chem_v3_native.cpp',N/'source/src/mfem_pic_full_time_chemistry.inc',
       N/'source/src/mfem_pic_full_time.inc',N/'source/src/mfem_pic_arena.cpp',
       N/'source/include/raikou/mfem_pic_full_time.h',N/'source/include/raikou/mfem_gap_pic.h',
       S/'native/mfem_electrode_heat_source_v7/source_owner.hpp',
       S/'native/mfem_electrode_heat_source_v7/surface_model_v7.hpp',
       T/'source/phase_law_host.hpp',T/'source/phase_law.h',T/'source/phase_owner.hpp',
       ROOT/'pic-electrode-heat-feedback-20261029/native/mfem_electrode_heat_feedback_v2/water_scenario.hpp',
       ROOT/'pic-electrode-heat-feedback-20261029/native/mfem_electrode_heat_feedback_v2/source-kernel/include/surface_emission.hpp',
       ROOT/'pic-water-event-bridge-20261029/source/water_parent_bridge.hpp',
       ROOT/'pic-arrival-owner-v1.1-20261030/include/pf_arrival_v11/native_v6_births.hpp',
       ROOT/'hall-corrective-v3.3-20261031/source/census_charge.hpp',
       ROOT/'hall-corrective-v3.3-20261031/source/motion_receipt.hpp',
       ROOT/'native-gas-adapter-v1-20261101/include/native_gas_stage.hpp']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
pins={str(p):sha(p) for p in files}
blob=files[0].read_bytes();manifest=json.loads(files[1].read_text())
expected='70fcc2adc2fdaa62abed7ed02d68ffc09b1b18554eda8fa99c5451120ae454b8'
if sha(files[0])!=expected or manifest['blob_sha256']!=expected:raise RuntimeError('wrong sealed blob')
if blob[:8]!=b'PFCHEMV3':raise RuntimeError('wrong blob magic')
version,nsections=struct.unpack_from('<II',blob,8)
if version!=1:raise RuntimeError('unknown blob version')
sections={}
for k in range(nsections):
 name,offset,size=struct.unpack_from('<24sQQ',blob,16+40*k)
 if offset+size>len(blob):raise RuntimeError('section outside blob')
 sections[name.rstrip(b'\0').decode()]=blob[offset:offset+size]
names=sections['strings'].decode().splitlines()
Q=1.602176634e-19;C2=299792458.0*299792458.0
rows=[]
for i,s in enumerate(manifest['species']):
 row=struct.unpack_from('<ddiI6I4I',sections['species'],64*i)
 if (row[0],row[1],row[2],names[i])!=(s['mass_eV'],s['ref_eV'],s['charge'],s['name']):
  raise RuntimeError('manifest differs from actual blob species row')
 rows.append({'index':i,'state':'pfc-v3/'+s['name'],'mass_eV':row[0],'mass_kg':row[0]*Q/C2,
              'reference_eV':row[1],'reference_J':row[1]*Q,'charge_e':row[2],
              'nuclei_H_O_Fe_Cr_Ni_spare':list(row[4:10]),'static_internal_J':0.0})
old={'pfc-v3/e':9.1093837139e-31,'pfc-v3/H2O':2.9915e-26,
     'pfc-v3/H2O+':2.9915e-26-9.1093837139e-31+12.621*Q/C2}
compat=[]
for s in rows:
 if s['state'] in old:
  a=old[s['state']];compat.append({'state':s['state'],'pilot_mass_kg':a,'blob_mass_kg':s['mass_kg'],
                                 'relative_pilot_minus_blob':a/s['mass_kg']-1,
                                 'bit_equal':struct.pack('<d',a)==struct.pack('<d',s['mass_kg'])})
metals=[]
for name,m in [('Fe',55.845),('Cr',51.9961),('Ni',58.6934),('Mn',54.938043),('Si',28.085),('C',12.011)]:
 phase=m*1e-3/6.02214076e23;native=m*1.66053906660e-27
 chem=next((s['mass_kg'] for s in rows if s['state']=='pfc-v3/'+name),None)
 metals.append({'element':name,'phase_mass_kg':phase,'metal_v1_mass_kg':native,'chemistry_mass_kg':chem,
                'relative_phase_minus_native':phase/native-1,'relative_phase_minus_chemistry':phase/chem-1 if chem else None})
result={'schema':'native8-species-compatibility-audit-v1','blob_sha256':expected,
        'verified_binary_species_rows':len(rows),'species':rows,'pilot_compatibility':compat,
        'parent_ionization_Q_eV':12.621,'blob_water_ionization_Q_eV':rows[2]['reference_eV']-rows[1]['reference_eV'],
        'metal_mass_conventions':metals,'model_changes':False,'numerical_tests_run':False}
(OUT/'mapping.json').write_text(json.dumps(result,indent=2)+'\n')
(OUT/'INPUTS-SHA256.json').write_text(json.dumps({'files':pins},indent=2)+'\n')
print(json.dumps({'verified_species_rows':len(rows),'input_files':len(pins),'blob_sha256':expected}))
