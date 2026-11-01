import importlib.util
import json
from pathlib import Path
import sys

root=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('h2_primary',root/'src/h2_primary.py')
m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m)
c=m.load_primary(
 '/home/treefor/projects/runtimes/water-fragment-source-inventory-20261101/references/MCCC-el-H2-b3Su_DE.X1Sg_vi=0.txt',
 '/home/treefor/projects/runtimes/pic-water-chemistry-v3.2.2-auditprep-20261101')
result={'identity':c.identity(),'threshold_assessment':c.threshold_assessment(),
        'rate_array_redistributed':False,
        'full_network_admission':False,
        'remaining_missing_H2_channels':['elastic','rotational','vibrational',
            'other electronic excitation/dissociation','attachment','resolved ionization products'],
        'remaining_energy_coverage':'No below-grid or above-150-eV rate; source/target-frame mapping unresolved',
        'other_missing_groups':['e+OH','e+O'],
        'unknown_recombination_and_existing_region_guards':'unchanged',
        'event_model':'No product energy-sharing, scattering-angle or recoil law admitted'}
print(json.dumps(result,indent=2,sort_keys=True))
