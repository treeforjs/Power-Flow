"""Long resident runs and carrier scaling on one fixture, candidate003 vs v5.
Run under core-lease-20261030/any-core.sh --wait (one core, nice 15, one thread)
with CUDA_VISIBLE_DEVICES empty. Each process runs through the frozen launcher
under RLIMIT_AS 3 GiB. Requires build-components.py and build-baseline-c003.py.
Usage: long-run.py [steps-for-v5 (default 1000000)]"""
import sys
from tooling import *

resource_guard()
steps = int(sys.argv[1]) if len(sys.argv) > 1 else 1000000
out = OUT / 'long-run'
r = Runner(out, 'long-run')
v5 = P / 'build/test_long_run'
c3 = P / 'build/c003/long_run_c003'
pins = {'v5_binary_sha256': sha256(v5), 'c003_binary_sha256': sha256(c3),
        'v5_source_identity': json.loads((P / 'SOURCE-SHA256.json').read_text())['source_identity'],
        'c003_source_identity': json.loads((C003 / 'SOURCE-SHA256.json').read_text())['source_identity']}
runs = []
# Same fixture: 2 carriers, zero-event exposures, genuine Hall v2 motion at zero field.
for label, binary, n, sample, extra in [('c003-2-carriers', c3, 2000, 10, []),
                                        ('v5-2-carriers-818', v5, 818, 1, ['1024'])]:
    result = out / (label + '.json')
    r.run([binary, result, n, sample] + extra, label + '.log', timeout=3600)
    runs.append(label)
for n_carriers, n_steps in [(10000, 6), (100000, 3)]:
    for label, binary in [('c003', c3), ('v5', v5)]:
        name = f'{label}-{n_carriers}-carriers'
        r.run([binary, out / (name + '.json'), n_steps, 1, '1024', n_carriers], name + '.log', timeout=3600)
        runs.append(name)
big = out / 'long-run-v5-1e6.json'
reuse = False
if big.exists() and (out / 'binary-sha256.txt').exists():
    reuse = (out / 'binary-sha256.txt').read_text().split()[0] == pins['v5_binary_sha256'] and \
        json.loads(big.read_text())['completed_steps'] == steps
if not reuse:
    r.run([v5, big, steps, max(1, steps // 100), '1024'], 'long-run-v5-1e6.log', timeout=7000)
runs.append('long-run-v5-1e6')


def load(name):
    return json.loads((out / (name + '.json')).read_text())


def mean_at(samples, upto):
    rows = [s for s in samples if s['step'] <= upto]
    return rows[-1]['mean_step_s'] if rows else None


c = load('c003-2-carriers')
v = load('long-run-v5-1e6')
v818 = load('v5-2-carriers-818')
summary = {
    'fixture': c['fixture'], 'pins': pins, 'reused_existing_v5_long_run': reuse,
    'candidate003': {'completed_steps': c['completed_steps'], 'refusal': c['refusal'],
                     'mean_step_s': c['mean_step_s'], 'native_mean_step_s': c['native_mean_step_s'],
                     'hall_provider_mean_step_s': c['hall_provider_mean_step_s'],
                     'mean_step_s_steps_1_to_10': mean_at(c['samples'], 10),
                     'mean_step_s_steps_801_to_810': c['samples'][-1]['mean_step_s'] if c['samples'] else None,
                     'peak_VmHWM_KiB': c['peak_VmHWM_KiB']},
    'v5_first_818_steps': {'completed_steps': v818['completed_steps'], 'mean_step_s': v818['mean_step_s'],
                           'native_mean_step_s': v818['native_mean_step_s'],
                           'hall_provider_mean_step_s': v818['hall_provider_mean_step_s'],
                           'peak_VmHWM_KiB': v818['peak_VmHWM_KiB']},
    'v5_long': {k: v[k] for k in ['completed_steps', 'refusal', 'mean_step_s', 'native_mean_step_s',
                                   'hall_provider_mean_step_s', 'peak_VmHWM_KiB', 'ru_maxrss_KiB', 'wall_s',
                                   'history', 'anchor_chain_verified', 'anchor_chain_verification_s',
                                   'final_accepted_step', 'final_accepted_arena_revision', 'final_time_s',
                                   'forbidden_native_solver_or_advance_calls']},
    'v5_long_rss_KiB_first_last_sample': [v['samples'][0]['VmRSS_KiB'], v['samples'][-1]['VmRSS_KiB']],
    'v5_long_block_mean_step_s_first_last': [v['samples'][0]['mean_step_s'], v['samples'][-1]['mean_step_s']],
    'scaling': {}}
for n_carriers in [10000, 100000]:
    for label in ['c003', 'v5']:
        j = load(f'{label}-{n_carriers}-carriers')
        summary['scaling'][f'{label}-{n_carriers}'] = {
            k: j[k] for k in ['completed_steps', 'refusal', 'mean_step_s', 'native_mean_step_s',
                              'hall_provider_mean_step_s', 'peak_VmHWM_KiB', 'phase_totals_s']}
summary['utc'] = utc()
(out / 'long-run-summary.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps({k: summary[k] for k in ['candidate003', 'v5_first_818_steps']}, indent=1))
