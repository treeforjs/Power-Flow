"""Compile identical bounded fixture against the two actual release archives.
No archive is modified, no provenance identity is normalized."""
from tooling import *
resource_guard(build=True)
base=RT/'native-full-time-candidate008-auditfix'
pins=json.loads((P/'SUCCESSOR-BASE.json').read_text())
baseline_identity='5f7d171c62708cf3e2e03f0fd843e7c906ed65a5db7d5955624773dee3516f2b'
baseline_library='098c1aaa14a533877ba4474f1ffae432994329db809b12b2f1391c6e2a65f9d7'
need(sha256(base/'build/libpic_native_v8.a')==baseline_library,'qualified baseline archive changed')
inputs=sealed_inputs()
manifest=json.loads((P/'SOURCE-SHA256.json').read_text())
out=OUT/'build-differential'; r=Runner(out,'build'); B=P/'build'
inc=['-I'+str(P/'source/include'),'-I'+str(P/'forms/include'),'-I'+str(P/'forms/gpu'),
     '-I'+str(GPU_COMMON),'-isystem',str(P/'inputs'),'-I'+str(VENDOR),'-isystem',str(CHEM/'include/pf_chem_v3'),
     '-isystem',str(ARRIVAL/'include/pf_arrival'),'-isystem',str(MERGE_SPLIT/'include'),
     '-I'+str(P/'source/src'),'-isystem',str(HALL),'-I'+str(P/'tests'),'-I'+str(HALL_V2/'source'),
     '-I'+str(ARRIVAL/'include'),'-I'+str(P/'source/adapters'),'-I'+str(V3_COMMIT/'source'),
     '-DPF_ARRIVAL_EXACT_COUNT_HEADER="raikou/exact_particle_count.h"','-Wno-c++26-extensions',
     '-isystem',str(CUDA_INC),'-isystem',str(NVRTC_DIR/'include')]+SYSTEM
for label, root, identity in [('baseline',base,baseline_identity),('candidate',P,manifest['source_identity'])]:
    obj=B/(label+'-birth-memo-diff.o')
    r.run(['/usr/bin/mpic++','-std=c++20','-O2','-g0','-ffp-contract=off','-fno-fast-math',
           '-ffunction-sections','-fdata-sections','-Wall','-Wextra','-Wpedantic','-Werror',
           '-DPIC_REACTION_NATIVE_SOURCE_SHA256="'+identity+'"',*inc,'-c',P/'tests/test_birth_memo_diff.cpp','-o',obj],
          label+'-compile.log',timeout=180)
    r.run(['/usr/bin/mpic++','-Wl,--gc-sections',obj,base/'build/kinetic-boundary.o',root/'build/libpic_native_v8.a',
           *LIBS,'-o',B/(label+'-birth-memo-diff')],label+'-link.log',timeout=120)
need(sealed_inputs()==inputs,'frozen dependencies changed')
for name,digest in manifest['files'].items(): need(sha256(P/name)==digest,'candidate source changed')
need(sha256(base/'build/libpic_native_v8.a')==baseline_library,'baseline archive changed')
(out/'build-result.json').write_text(json.dumps({'baseline_identity':baseline_identity,'baseline_library':baseline_library,
    'candidate_identity':manifest['source_identity'],'candidate_library':sha256(B/'libpic_native_v8.a'),
    'fixture_sha256':sha256(P/'tests/test_birth_memo_diff.cpp'),'commands':r.records},indent=2)+'\n')
print('actual-archive differential fixture compiled')
