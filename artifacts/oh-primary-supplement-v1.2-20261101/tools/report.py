import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from decimal import Decimal

ROOT=Path(__file__).resolve().parents[1]
REF=Path('/home/treefor/projects/runtimes/water-fragment-source-inventory-20261101/references')
MANIFEST=Path('/home/treefor/projects/runtimes/pic-water-chemistry-v3.2.2-auditprep-20261101/data/tables/chemistry-v3.2-manifest.json')
spec=importlib.util.spec_from_file_location('oh_primary',ROOT/'src/oh_primary.py')
module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
channels=module.inspect(REF/'PSSTab264csupp_e-OH-data.xlsx')
raw=MANIFEST.read_bytes()
assert hashlib.sha256(raw).hexdigest()=='decc4f78e3f6ebd8c51eda9e662bb72d36597124c0b7026bea46dbc31c9ad9bd'
manifest=json.loads(raw);species={s['name']:s for s in manifest['species']}
ground_cost=(species['O']['ref_eV']+species['H']['ref_eV'])-species['OH']['ref_eV']
dea_cost=(species['O-']['ref_eV']+species['H']['ref_eV'])-species['OH']['ref_eV']
def below(channel, threshold):
    # Count an incompatibility, not an accepted collision probability.
    t=Decimal.from_float(threshold)
    return sum(1 for e,v in channels[channel].rows if e<t and v>0)
result={
 'schema':'oh-primary-numerical-inventory-v1.2',
 'source_doi':'10.1088/1361-6595/ab364c',
 'workbook_sha256':module.WORKBOOK_SHA256,
 'publisher_data_url':'https://iopscience.iop.org/article/10.1088/1361-6595/ab364c/data',
 'publisher_filename':'PSSTab264csupp_e-OH-data.xlsx',
 'license':'Publisher supplement page explicitly states CC BY 3.0 with authors/title/journal/DOI attribution; raw workbook remains local in this increment.',
 'channels':{k:v.summary() for k,v in channels.items()},
 'compatibility':{
  'registry_ground_O_H_reference_cost_eV':ground_cost,
  'registry_Ominus_H_reference_cost_eV':dea_cost,
  'paper_dissociation_D0_eV':3.88,
  'paper_DEA_model_OH_electron_affinity_eV':1.56,
  'dissociation_positive_nodes_below_registry_ground_product_cost':below('dissociation',ground_cost),
  'DEA_positive_nodes_below_registry_product_cost':below('DEA_B',dea_cost),
  'excited_O_1D_registry_state_present':any(s in species for s in ['O(1D)','O_1D']),
  'dissociation_branch_resolution':'Workbook has one column without a product-state label. Paper Figure 5 discusses ground O(3P) and excited O(1D); whether the column is one branch or a sum is unverified.',
  'dissociation_model':'Published estimate assuming electronic excitation into specified repulsive states produces dissociation; not an experimentally validated complete cross section.',
  'DEA_model':'Published approximate resonance/Morse-potential estimate; paper assumes only O-minus production. No threshold clamp or product law adopted.',
  'elastic_model':'Fixed-equilibrium-geometry integral cross section with Born high-partial-wave correction; not a differential scattering sampler.',
  'rate_sheet_units':'Workbook prints bohr^2, while primary paper Figure 4 labels cm^3/s. Do not silently reassign units; normalization/source-column correspondence remains to be verified.',
  'source_energy_to_native_target_frame_qualified':False,
  'product_state_threshold_and_differential_law_qualified':False},
 'remaining_gaps':['OH ionization not supplied by this workbook','sub-grid and upper-energy tails unresolved','product energy sharing and angles unresolved','excitation/product-state registry mapping unresolved','strict unknown/recombination guards unchanged'],
 'production_qualified':False,'physical_model_validated':False,'network_complete':False,'native_adopted':False,
 'raw_arrays_redistributed':False}
print(json.dumps(result,sort_keys=True,indent=2))
