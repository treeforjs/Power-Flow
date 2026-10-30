#!/usr/bin/env python3
"""Independent read-only audit of completed parent-owned GPU runs. No CUDA."""
import hashlib
import json
import math
from pathlib import Path
import statistics
import struct

R = Path('/home/treefor/projects/runtimes/pic-cuda-gather-investigation-20261029')
V = R / 'backend-v2'
H = V / 'lifecycle-harness'
O = Path('<workspace>/outputs/pic-profile/cuda-gather-investigation')
G = O / 'gpu-lifecycle-v2'
F = O / 'actual-midpoint-fixture'
REVIEW = O.parent / 'cuda-gather-independent-review/backend-v2/lifecycle-harness/disposition.json'
HERE = Path(__file__).resolve().parent
checks = []

def require(condition, label):
    checks.append({'check': label, 'passed': bool(condition)})
    if not condition:
        raise AssertionError(label)

def load(path):
    return json.loads(Path(path).read_text())

def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

def hold(d, label):
    require(d['production_host_pid'] == 278902 and d['start_ticks'] == '5288339'
            and d['state'] == 'T' and d['cgroup_frozen']
            and d['production_namespace_pids'] == [278902, 3]
            and d['exe'] == '/home/treefor/projects/runtimes/mfem-gpu-pic-20261029/qualified/heating-history-pic-gpu'
            and d['cgroup'] == '/powerflow-pic-gpu-reserved', label)

def cleanup(d, label, active):
    if active:
        require(d['attempted'] and d['context_established'] and d['stream_drained']
                and d['resources_released'] and not d['resources_retained']
                and d['caller_context_restored'] and d['primary_reference_released']
                and d['issue_count'] == 0 and not d['issues'], label)
    else:
        require(not any(d[k] for k in d if k != 'issues') and not d['issues'], label + ' (unused sink stays empty)')

def stats(values):
    ordered = sorted(values)
    def percentile(fraction):
        at = (len(ordered)-1)*fraction
        i = int(at)
        return ordered[i] + (ordered[min(i+1, len(ordered)-1)]-ordered[i])*(at-i)
    return {'median': statistics.median(values), 'mean': statistics.mean(values),
            'min': min(values), 'max': max(values), 'p10': percentile(.1), 'p90': percentile(.9)}

def main():
    review = load(REVIEW)
    require(review['review_complete'] and not review['blocking_findings'], 'fresh harness review accepted')
    for key, path in [('artifact_manifest_sha256', H/'artifact-sha256.json'),
                      ('review_request_sha256', H/'REVIEW-REQUEST.json'),
                      ('parent_command_request_sha256', H/'parent-gpu-request.json'),
                      ('force_lifecycle_executable_sha256', H/'build/force_lifecycle')]:
        require(sha(path) == review[key], 'accepted artifact identity: '+path.name)
    manifest = load(H/'artifact-sha256.json')
    require(all(sha(k) == v for k, v in manifest.items()), 'all 42 reviewed artifacts/dependencies unchanged')
    require(all(sha(k) == v for k, v in load(V/'v1-baseline-sha256.json').items()), 'all v1 evidence unchanged')
    require(all(sha(F/k) == v for k, v in load(F/'sha256.json').items()), 'actual fixture evidence unchanged')
    core_review = load(REVIEW.parents[1]/'disposition.json')
    require(all(sha(V/k) == v for k,v in core_review['accepted_source_hashes'].items()), 'accepted v2 core source unchanged')
    dispatch, continuation = load(G/'PARENT-DISPATCH.json'), load(G/'PARENT-CONTINUATION.json')
    require(dispatch['request_sha256'] == review['parent_command_request_sha256']
            and dispatch['artifact_sha256'] == review['artifact_manifest_sha256'], 'parent dispatch bound to accepted bundle')
    cases = dispatch['cases'] + continuation['cases']
    require([x['phase'] for x in cases] == ['shared-lifecycle','explicit-standalone','reference-inclusive-timing']
            and all(x['returncode'] == 0 and 0 < x['wall_s'] < 60 for x in cases), 'three serial processes returned0 within bounds')
    require(not dispatch['production_resume_authorized'] and not dispatch['application_integration_authorized']
            and not continuation['production_resume_authorized'] and not continuation['application_integration_authorized'], 'parent kept resume/integration unauthorized')
    require(dispatch['started_unix_s'] < continuation['started_unix_s'] < continuation['finished_unix_s'], 'continuation follows original shared dispatch; no shared rerun required')
    requests = load(H/'parent-gpu-request.json')['commands']
    for case, request in zip(cases, requests):
        phase = case['phase']; out = G/phase
        guard, pre = load(out/'host-guard.json'), load(out/'preflight.json')
        argv=request['command_argv']
        require(guard['command'] == argv[argv.index('--')+1:], phase+' exact reviewed child argv')
        require(guard['logical_cpu'] == [7] and guard['nice'] >= 15
                and guard['production_pidfd_retained_until_child_exit']
                and not guard['production_signalled'] and not guard['gpu_actions_in_guard']
                and guard['child_returncode'] == 0 and guard['guard_error'] is None
                and guard['periodic_checks'] >= 2, phase+' CPU/hold/exit guard')
        for name in ['before','after']:
            hold(guard[name], phase+' '+name+' stopped/frozen identity')
            require(guard[name]['pidfd_retained'], phase+' '+name+' retained pidfd')
        require(guard['before']['checked_unix_s'] < guard['after']['checked_unix_s'], phase+' hold chronology')
        require(len(guard['child_requests']) == 1 and guard['child_requests'][0]['stage'] == 'before_cpp_gpu_exec', phase+' child requested live pre-CUDA guard')
        hold(guard['child_requests'][0]['check'], phase+' live child hold')
        require(pre['passed'] and not pre['production_signalled'], phase+' preflight passed')
        hold(pre['hold'], phase+' preflight hold')
        require(pre['hold']['host_proc_readonly'] and pre['hold']['distinct_pid_views']
                and pre['hold']['host_executable_checked_with_retained_pidfd'], phase+' namespace/executable boundary')
        require((G/(phase+'.stderr')).stat().st_size == 0 and (G/(phase+'.stdout')).stat().st_size == 0, phase+' no unexpected process output')
    shared, standalone = load(G/'shared-lifecycle/lifecycle.json'), load(G/'explicit-standalone/lifecycle.json')
    for d, label, active in [(shared,'shared',{'shared','wrong','lost'}), (standalone,'standalone',{'standalone'})]:
        require(d['passed'] and not d['original_error'] and d['particles'] == 46726 and d['forces_completed'] == 2
                and d['GPU_kernel_dispatched'] and not d['full_PIC_replay'] and not d['production_actions']
                and not d['device_fault_injection'], label+' successful actual-midpoint scope')
        require(d['reset_attempts'] == 0 and d['borrowed_destroy_attempts'] == 0, label+' no forbidden context operation')
        for name, sink in d['cleanup'].items():
            cleanup(sink, label+' cleanup '+name, name in active)
    require(shared['successful_own_retains'] == shared['successful_own_releases'] == 3
            and shared['owned_context_destroys'] == 1
            and shared['wrong_current_rejected_before_dispatch'] and shared['lost_current_rejected_before_dispatch'], 'shared reference balance and negative-context outcomes')
    expected = ['solve1_before_force','force1','solve2_between_forces','force2','shared_explicit_close','solve3_after_force_close',
                'wrong_current_selected','wrong_close_restored_owned_context','wrong_owned_context_destroyed_only',
                'solve4_after_wrong_context_probe','lost_current_selected','lost_close_restored_null','solve5_after_lost_context_probe']
    require([x['stage'] for x in shared['events']] == expected, 'exact shared solve/force/close/negative-probe sequence')
    events = {x['stage']: x for x in shared['events']}; borrowed=shared['borrowed_primary']
    require(borrowed != '0x0' and all(x['primary_active'] == 1 for x in shared['events']), 'shared primary stays active')
    require(all(events[k]['current'] == borrowed for k in expected if k not in {'wrong_current_selected','wrong_close_restored_owned_context','lost_current_selected','lost_close_restored_null'}), 'shared identity around all solves/forces/normal close')
    for before, after in [('wrong_current_selected','wrong_close_restored_owned_context'),('lost_current_selected','lost_close_restored_null')]:
        require(events[before]['current'] == events[after]['current'] and
                all(events[before][k] == events[after][k] for k in ['launches','H2D_dispatches','D2H_dispatches']), before+' rejected before dispatch and caller restored')
    require(events['wrong_current_selected']['current'] not in {borrowed,'0x0'}
            and events['lost_current_selected']['current'] == '0x0', 'distinct wrong and missing current contexts exercised')
    require(len(shared['solves']) == 5 and shared['manufactured_Poisson_rows'] == 130321, 'five manufactured actual-resolution cuDSS solves')
    for i, s in enumerate(shared['solves']):
        limit = max(1e-30, 2e-13*s['rhs_norm'])
        require(all(math.isfinite(s[k]) for k in ['rhs_norm','residual_norm','guard_limit','relative_residual'])
                and s['guard_limit'] == limit and 0 <= s['residual_norm'] <= limit
                and math.isclose(s['relative_residual'],s['residual_norm']/s['rhs_norm'],rel_tol=1e-15), 'unchanged finite/residual gate solve'+str(i+1))
    require(standalone['successful_own_retains'] == standalone['successful_own_releases'] == 1
            and standalone['standalone_primary_active_before'] == standalone['standalone_primary_active_after'] == 0
            and standalone['events'][0]['current'] == standalone['events'][-1]['current'] == '0x0'
            and standalone['owned_context_destroys'] == 0, 'explicit standalone reference balance/null restoration/activity lifetime')
    for phase, names in [('shared-lifecycle',['force1','force2']),('explicit-standalone',['standalone_force1','standalone_force2']),('reference-inclusive-timing',[''])]:
        for name in names:
            prefix = (name+'-') if name else ''
            comparison = load(G/phase/(prefix+'comparison.json'))
            require(comparison['passed'] and comparison['failures'] == 0
                    and all(v == 0 for k,v in comparison.items() if k not in {'passed','failures'}), phase+' '+name+' exact reference comparison')
            data=(G/phase/(prefix+'device.bin')).read_bytes(); ref=(F/'reference.bin').read_bytes()
            require(len(data) == len(ref) == 46726*80, phase+' '+name+' complete output shape')
            # Ignore only C++ trailing struct padding, independently check every
            # E/B/u FP64 bit and status rather than trust the emitted receipt.
            require(all(data[i:i+76] == ref[i:i+76] for i in range(0,len(data),80)), phase+' '+name+' every E/B/u/status bit independently equals captured original')
            require(all(all(math.isfinite(v) for v in struct.unpack_from('<9d',data,i))
                        and struct.unpack_from('<i',data,i+72)[0] == 0 for i in range(0,len(data),80)), phase+' '+name+' independently finite and status0')
    timing = load(G/'reference-inclusive-timing/reference-inclusive-timing.json')
    require(timing['passed'] and timing['particles'] == 46726 and not timing['application_integrated']
            and not timing['python_force_loop'] and all(timing[k] for k in ['fresh_independent_CPU_reference_every_batch','reference_generation_charged','reference_gates_charged','commit_charged'])
            and len(timing['samples_ms']) == 20, 'twenty measured complete reference-inclusive force samples')
    tc = load(G/'reference-inclusive-timing/cleanup.json')
    require(tc['issue_count'] == 0 and not tc['issues'] and tc['resources_released'] and not tc['resources_retained']
            and tc['caller_context_restored'], 'timing cleanup recorded fields pass (full drain/ref fields not emitted by adapter)')
    rows = timing['samples_ms']
    for i, s in enumerate(rows):
        require(all(math.isfinite(x) and x >= 0 for x in s.values())
                and abs(s['force_pipeline_total'] - s['native_reference_generation'] - s['GPU_API_including_pack_and_gates'] - s['commit']) < 1e-10,
                'measured total partitions exactly sample'+str(i+1))
    summary={key:stats([s[key] for s in rows]) for key in rows[0]}
    paired=stats([(s['force_pipeline_total']/s['native_reference_generation']-1)*100 for s in rows])
    all_receipts={str(p):sha(p) for p in sorted(G.rglob('*')) if p.is_file()}
    pid=Path('/proc/278902');stat=pid.joinpath('stat').read_text().split(') ',1)[1].split()
    require(stat[0]=='T' and stat[19]=='5288339' and str(pid.joinpath('exe').readlink())=='/home/treefor/projects/runtimes/mfem-gpu-pic-20261029/qualified/heating-history-pic-gpu'
            and 'frozen 1' in Path('/sys/fs/cgroup/powerflow-pic-gpu-reserved/cgroup.events').read_text(), 'production remains original stopped/frozen process at audit')
    result={'passed':True,'audit_checks':len(checks),'GPU_executed_by_auditor':False,'production_actions':False,
            'receipt_hashes':all_receipts,'checks':checks,'timing_ms':summary,'paired_total_over_fresh_reference_percent':paired,
            'setup_excluded_ms':timing['setup_excluded_ms'],
            'shared_max_relative_residual':max(s['relative_residual'] for s in shared['solves']),
            'qualification':'Real midpoint E/B/u plus shared-primary lifecycle passed for reviewed isolated backend. Full evolving-step integration remains unqualified.',
            'limits':['Poisson RHS is manufactured; it tests real solver/context lifetime, not live PIC charge deposition.',
                      'Timing adapter emits only subset of cleanup booleans; full drain/own-release evidence comes from separately reviewed lifecycle processes, not those missing fields.',
                      'Matched baseline here is fresh native reference generation within each CPU7 sample, not a native whole force-only update or full timestep.',
                      'Setup/static refresh amortization, evolving fields/capacities, both PIC drifts/solves/sources are excluded. No application speedup claim.']}
    (HERE/'audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'passed':True,'audit_checks':len(checks),'timing_medians_ms':{k:v['median'] for k,v in summary.items()},'paired_total_over_reference_percent':paired,'setup_excluded_ms':timing['setup_excluded_ms']}))

if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        (HERE/'audit-failure.json').write_text(json.dumps({'passed':False,'error':str(error),'checks':checks},indent=2)+'\n')
        raise
