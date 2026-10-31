#!/usr/bin/python -I
"""Print the key numbers of gpu-vacuum.json files (diagnostic helper)."""
import json, sys
for p in sys.argv[1:]:
    d = json.load(open(p))
    ph = d.get('phases', {})
    def med(k):
        return ph.get(k, {}).get('median')
    print(p.split('/')[-2], 'spec', d['spec'], 'done', d.get('steps_done'), 'refusal', d.get('refusal', '')[:100])
    print('  ms/step', round(d.get('ms_per_step', 0), 3), 'median total', med('step.total_ms'), 'prepare', med('step.prepare_vacuum_ms'),
          'preflight', med('step.preflight_vacuum_final_ms'), 'iters', ph.get('solve.iterations', {}).get('median'), ph.get('solve.iterations', {}).get('max'))
    print('  host:', {k[8:]: round(v['median'], 3) for k, v in ph.items() if k.startswith('prepare.')})
    print('  device:', {k[7:]: round(v['median'], 3) for k, v in ph.items() if k.startswith('device.')})
    print('  guards', d.get('guards'), 'control', d.get('control'))
    print('  memory peak MB', d['memory']['own_peak_device_bytes'] / 1e6, 'model MB', d.get('memory_model_bytes', 0) / 1e6, 'pinned MB',
          d['memory']['own_pinned_host_bytes'] / 1e6, 'rss MB', d.get('max_rss_kB', 0) / 1e3)
    print('  gauss', d.get('gauss_initial'), '->', d['gauss_series'][-1]['gauss'], 'rate/s', d.get('gauss_ratio_rate_per_s'),
          'proj refusal', d.get('gauss_ratio_projected_refusal_time_s'))
    print('  B', d.get('B_bounds_final'), 'I', d.get('final_current_A'), 'cpu model ms', round(d.get('cpu_model_construct_ms', 0)),
          'gpu hold ms', round(d.get('gpu_hold_ms', 0)))
