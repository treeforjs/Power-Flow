"""Read preserved electron-only snapshot; prioritize searches, never qualify chemistry."""
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
BASE=Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs')
DIRECTORY=BASE/'full60-gpu-manual-pause-checkpoint'
COMPARATOR=BASE/'pic-profile/cpu-independent-review/checkpoint-validation/compare_saved_pic.py'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify(pins):
    for section in ('source_files','dependencies'):
        for n,h in pins[section].items():assert sha(ROOT/n)==h,n
    assert sha(Path(sys.executable).resolve())==pins['frozen_python_sha256']
def main():
    job=os.environ.get('SIMQ_JOB_ID');assert job and job.isdigit() and len(os.sched_getaffinity(0))==1
    assert not os.environ.get('CUDA_VISIBLE_DEVICES')
    pins=json.loads((ROOT/'INPUTS.json').read_text());verify(pins)
    spec=importlib.util.spec_from_file_location('historical_comparator',COMPARATOR);m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m)
    spec_reader=importlib.util.spec_from_file_location('pause_snapshot_reader',ROOT/'tools/snapshot_reader.py');reader=importlib.util.module_from_spec(spec_reader);spec_reader.loader.exec_module(reader)
    seed=m.strict_json(DIRECTORY/'seed.json');reader.validate_seed(seed)
    bad=[dict(seed,scope='full physical evolution'),dict(seed,checkpoint='../checkpoint.bin'),dict(seed,start_ns=float('nan')),dict(seed,checkpoint_sha256='f'),dict(seed,unexpected=1),{k:v for k,v in seed.items() if k!='scope'}]
    for trial in bad:
        try:reader.validate_seed(trial)
        except m.Invalid:pass
        else:raise AssertionError('malformed seed admitted')
    state=reader.read_state(DIRECTORY,m.schema())
    pause=json.loads((BASE/'full60-manual-pause.json').read_text())
    assert sha(state.checkpoint)==pause['checkpoint']['sha256']
    start=pause['last_reported_progress']['start_time_s']
    assert state.history[0]['step']==0 and state.history[0]['time_s']==start
    # Original driver reports start + native relative clock; restart.meta stores the relative clock.
    accepted=start+state.ledger['time_s']
    assert accepted==pause['checkpoint']['accepted_time_s']
    assert abs(state.seed_info['start_ns']-accepted*1e9)<=math.ulp(state.seed_info['start_ns'])
    assert all(p['charge_C']==-m.E_CHARGE for p in state.particles)
    rows=[(m.kinetic(p)/m.E_CHARGE,p['weight'],p['id']) for p in state.particles]
    total=math.fsum(w for e,w,i in rows);assert total>0
    bounds=[0,7.3605,10,13.018813320400712,18.009473168203966,100,150,200,1000,5000,10000,110000,220000]
    bins=[]
    for a,b in zip(bounds,bounds[1:]):
        selected=[r for r in rows if a<=r[0]<b]
        bins.append({'lower_eV':a,'upper_eV':b,'upper_exclusive':True,'macroparticles':len(selected),'represented_fraction':math.fsum(r[1] for r in selected)/total})
    high=[r for r in rows if r[0]>=bounds[-1]]
    bins.append({'lower_eV':bounds[-1],'upper_eV':None,'macroparticles':len(high),'represented_fraction':math.fsum(r[1] for r in high)/total})
    assert sum(b['macroparticles'] for b in bins)==len(rows)
    assert abs(math.fsum(b['represented_fraction'] for b in bins)-1)<1e-14
    assert min(e for e,w,i in rows)>=0
    sorted_rows=sorted(rows);quantiles={};acc=0
    for e,w,i in sorted_rows:
        acc+=w
        for q in (0.5,0.9,0.99,0.999):
            if str(q) not in quantiles and acc>=q*total:quantiles[str(q)]=e
    report={'schema':'historical-electron-energy-search-priority-v1.3','seed_negative_controls':len(bad),'job':job,'allocated_cpu':sorted(os.sched_getaffinity(0)),
            'source_identity':pins['source_identity'],'input_checkpoint_sha256':sha(state.checkpoint),
            'seed_display_ns_rounding_difference':state.seed_info['start_ns']-accepted*1e9,'seed_display_tolerance_ulp':1,'native_relative_clock_s':state.ledger['time_s'],'heating_start_time_s':start,'accepted_snapshot_time_s':accepted,'accepted_snapshot_step':state.step,
            'macroparticles':len(rows),'represented_electrons':total,'energy_frame':'laboratory; equal to collision energy only for stationary targets',
            'minimum_eV':min(e for e,w,i in rows),'maximum_eV':max(e for e,w,i in rows),
            'represented_mean_eV':math.fsum(e*w for e,w,i in rows)/total,'represented_energy_quantiles_eV':quantiles,'bins':bins,
            'provisional_search_envelope_eV':220000,'envelope_basis':'conditional recalled 220 kV feed, 110 kV midpoint, full 25 mm target; single-transit field-work scale, not a proven trajectory energy bound',
            'limits':['one paused electron-only snapshot, not a populated coupled plasma EEDF','conditional voltage waveform and geometry; reconstructed anode temperature, coarse B and height-averaged space charge',
                      'no ion/gas/collision feedback; cannot establish ion-target reachability or future electron energy bounds','only provisional prioritization of source searches; never a substitute for current all-species/thermal-support guards'],
            'simulation_advanced':False,'production_resumed':False,'chemistry_reachability_qualified':False,
            'production_qualified':False,'physical_model_validated':False,'network_complete':False}
    verify(pins);out=ROOT/'evidence'/('job-'+job);out.mkdir(exist_ok=False)
    (out/'report.json').write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
    print('PASS preserved snapshot identity, particle/weight partition and read-only checks')
if __name__=='__main__':main()
