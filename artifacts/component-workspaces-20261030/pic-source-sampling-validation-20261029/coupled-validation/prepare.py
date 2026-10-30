#!/usr/bin/env python3
"""Create only an owned, reproducible candidate; verify frozen source identity."""
import hashlib
import json
import os
from pathlib import Path

ROOT=Path(__file__).resolve().parent
SOURCE=Path('/home/treefor/projects/runtimes/mfem-gpu-pic-20261029/source')
SEALED=SOURCE.parent/'sealed-files.json'
OWNED=ROOT/'frozen'

def sha(data):return hashlib.sha256(data).hexdigest()
def main():
    if os.sched_getaffinity(0)!={6} or os.getpriority(os.PRIO_PROCESS,0)!=15:
        raise RuntimeError('CPU6/nice15 required')
    seal={r['path']:r['sha256'] for r in json.loads(SEALED.read_text())}
    paths=['include/raikou/mfem_gap_pic.h','src/mfem_gap_pic.cpp','src/bench_profile.hpp',
           'src/mfem_pic_distributed.hpp','src/cudss_poisson.hpp',
           'examples/mfem_heating_history_pic.cpp','examples/native_heating_fields.hpp']
    original={}
    for name in paths:
        source=SOURCE/name;data=source.read_bytes()
        if sha(data)!=seal[str(source)]:raise RuntimeError('Frozen seal mismatch: '+name)
        target=OWNED/name;target.parent.mkdir(parents=True,exist_ok=True)
        target.write_bytes(data);original[name]=sha(data)
    cpp=(OWNED/'src/mfem_gap_pic.cpp').read_text()
    header=(OWNED/'include/raikou/mfem_gap_pic.h').read_text()
    # Passive observer sees the unchanged native particle/hit before deletion.
    a='    void release_field_solver();'
    assert header.count(a)==1
    header=header.replace(a,a+'''
    void set_absorption_observer(std::function<void(const GapParticle&,int)> f) { absorption_observer_=std::move(f); }
''')
    a='    std::unique_ptr<PicPoissonState> poisson_;'
    assert header.count(a)==1
    header=header.replace(a,a+'\n    std::function<void(const GapParticle&,int)> absorption_observer_;')
    a='            const double energy=pic_kinetic_energy_J(p)*p.weight;'
    assert cpp.count(a)==1
    cpp=cpp.replace(a,a+'\n            if(absorption_observer_)absorption_observer_(p,hit);')
    assert cpp.count('#include "cudss_poisson.hpp"')==1
    cpp=cpp.replace('#include "cudss_poisson.hpp"','#include "gpu_disabled.hpp"')
    (ROOT/'mfem_gap_pic_candidate.cpp').write_text(cpp)
    (OWNED/'include/raikou/mfem_gap_pic.h').write_text(header)
    driver=(SOURCE/'examples/mfem_heating_history_pic.cpp').read_text()
    def take(start,end):
        a=driver.index(start);b=driver.index(end,a);return driver[a:b]
    structures=take('constexpr double eps=','\nvoid snapshot(')
    # These exact statements remain unchanged, including all safety constants.
    step=take('            double dt=pic.particles().empty()?', '            bool injected=false;')
    source=take('            bool injected=false;', '\n            if(!injected && pic.particles().empty())')
    generated='''#pragma once
#include "native_heating_fields.hpp"
#include <numeric>
#include <limits>
namespace validation_native {
'''+structures+'''
template<class Pic,class Field,class Bias>
auto allowed_step(Pic& pic,Field& fields,Bias& voltage,const Thermal& thermal,
 const std::array<Reservoir,4>& reservoir,const GapPicOptions& o,double t,double end,int next_frame) {
'''+step+'''
return std::pair<double,std::array<std::vector<double>,4>>{dt,std::move(T)};
}
template<class Source>
bool accumulate_source(Source& pic,std::array<Reservoir,4>& reservoir,
 std::array<std::vector<double>,4>& T,double dt,const GapPicOptions& o,std::uint64_t step) {
'''+source+'''
return injected;
}
} // namespace validation_native
'''
    (ROOT/'native_driver_parts.hpp').write_text(generated)
    old=ROOT.parent/'native_slice.hpp'
    # A separate source-only executable owns these verbatim native routines.
    (ROOT/'native_source_only.hpp').write_bytes(old.read_bytes())
    eq=json.loads(Path('<workspace>/outputs/pic-profile/source-sampling-validation/native-equivalence.json').read_text())
    assert sha(old.read_bytes())==eq['generated_sha256']
    source_proof={str(SOURCE/name):digest for name,digest in original.items()}
    manifest=dict(gpu_execution_enabled=False,coupled_execution_enabled=False,
      frozen_source_sha256=source_proof,
      verbatim_driver_slices={key:dict(sha256=sha(value.encode()),bytes=len(value.encode())) for key,value in [('structures',structures),('step_rules',step),('source_rules',source)]},
      source_sampler_sha256=sha(old.read_bytes()),
      mutations=['one passive absorption observer, no serialized ABI change','replace GPU header include with an unconditional throwing stub'],
      guard_files_unchanged=['src/mfem_pic_distributed.hpp'],
      generated_files_sha256={str(p.relative_to(ROOT)):sha(p.read_bytes()) for p in [ROOT/'mfem_gap_pic_candidate.cpp',ROOT/'native_driver_parts.hpp',ROOT/'native_source_only.hpp',OWNED/'include/raikou/mfem_gap_pic.h']})
    (ROOT/'candidate-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Prepared isolated candidate; native step/source slices verified; GPU disabled.')
if __name__=='__main__':main()
