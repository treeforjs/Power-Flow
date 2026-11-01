"""Compile-only mapping fixture against frozen v8.8 transaction declarations."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from tooling import *
resource_guard()
N=RT/'native-full-time-v8.8-provisional-target-owner-20261101'
G=Path('/home/treefor/projects/runtimes/gas-target-owner-v1.7.1-census-transport-20261101')
manifest=json.loads((N/'SOURCE-SHA256.json').read_text())
need(manifest['source_identity']=='9538ae6d7f1cffb4c0f446eba7df0449234a34eeb14169a1127d18e52853d6a7','frozen v8.8 identity')
source={name:sha256(N/name) for name in manifest['files']}
need(source==manifest['files'],'frozen v8.8 source changed')
inputs=sealed_inputs()
libraries={str(f):sha256(f) for f in map(Path,LIBS) if f.is_file()}
adoption=tree_hashes(P/'adoption')
gas=tree_hashes(G/'include')
out=OUT/sys.argv[1]
need(not out.exists(),'preserve evidence')
r=Runner(out,'adapter')
inc=['-I'+str(P/'adoption'),'-I'+str(N/'source/include'),'-I'+str(N/'source/src'),'-I'+str(N/'forms/include'),'-I'+str(N/'forms/include/pf_v7'),'-isystem',str(ARRIVAL/'include/pf_arrival'),'-isystem',str(CHEM/'include/pf_chem_v3'),'-isystem',str(G/'include'),'-isystem',str(N/'inputs'),'-I'+str(VENDOR)]+SYSTEM
r.run(['/usr/bin/mpic++','-std=c++20','-O1','-ffunction-sections','-fdata-sections','-Wl,--gc-sections','-ffp-contract=off','-fno-fast-math','-Wall','-Wextra','-Wpedantic','-Werror',*inc,P/'adoption/test_kind32.cpp',*LIBS,'-o',out/'test_kind32'],'compile.log',timeout=60)
r.run([out/'test_kind32'],'test.log',timeout=30)
need((out/'test.log').read_text().strip()=='kind32 mapping controls 10','actual controls')
need(libraries=={f:sha256(f) for f in libraries},'library changed')
need(source=={n:sha256(N/n) for n in source} and inputs==sealed_inputs() and adoption==tree_hashes(P/'adoption') and gas==tree_hashes(G/'include'),'input changed')
(out/'qualification.json').write_text(json.dumps({'passed':True,'checks':10,'transaction_header_identity':manifest['source_identity'],'native_source':source,'chemistry_inputs':inputs,'adoption':adoption,'gas_headers':gas,'libraries':libraries,'harness_sha256':sha256(__file__),'binary_sha256':sha256(out/'test_kind32'),'commands':r.records,'compile_only_adapter':True,'combined_native_gas_publication_qualified':False,'production_qualified':False},indent=2)+'\n')
print('kind32 adapter controls passed against frozen v8.8 headers')
