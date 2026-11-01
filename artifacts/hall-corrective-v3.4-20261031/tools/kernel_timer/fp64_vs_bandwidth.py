#!/usr/bin/python -I
"""FP64 throughput vs memory bandwidth for the main Hall v3.4 device kernels (measurement only).
Inputs: the kernel-timer tables of job 189 (run-N-2, run-N-5; steady step = (5-step - 2-step) / 3) and the static PTX
counts (ptx-counts.json). Bounds, not measurements of dynamic counts: per launch, the static instruction count of a
kernel is applied to every thread (thread-parallel kernels) or to every (node, carrier) record (node walks, 8 CIC records
per carrier), which over-counts the executed work (untaken branches, loop bodies counted once per record).
RTX 3080 Ti (GA102, 80 SMs, 1.665 GHz boost): FP64 pipe 2 instr/clk/SM = 266e9 instr/s; INT32 64 lanes/SM = 8.5e12/s;
DRAM 912e9 B/s; 1536 resident threads per SM = 122880 thread slots.
Usage: fp64_vs_bandwidth.py DIR OUT_JSON"""
import json, sys
from pathlib import Path

D = Path(sys.argv[1])
FP64_PEAK, INT_PEAK, BW_PEAK, SLOTS = 266e9, 8.5e12, 912e9, 122880
counts = {k: c for prog in json.loads((D / 'ptx-counts.json').read_text()).values() for k, c in prog.items()}
NODE_WALKS = {'hgp_reduce_end_nodes_v3', 'hgp_reduce_start_nodes', 'k34_charge_nodes', 'k34_support_nodes'}


def steady(n):
    a = {k['name']: k for k in json.loads((D / f'run-{n}-2/kernels.json').read_text())['kernels']}
    b = {k['name']: k for k in json.loads((D / f'run-{n}-5/kernels.json').read_text())['kernels']}
    z = {'launches': 0, 'device_ms': 0.0, 'threads': 0}
    return {k: {f: (v[f] - a.get(k, z)[f]) / 3 for f in z} for k, v in b.items()}


out = {}
for n in (100000, 1000000):
    rows, st = [], steady(n)
    total = sum(v['device_ms'] for v in st.values())
    for name, v in sorted(st.items(), key=lambda kv: -kv[1]['device_ms'])[:8]:
        c = counts.get(name)
        t = v['device_ms'] * 1e-3
        r = {'kernel': name, 'ms_per_step': round(v['device_ms'], 3), 'share': round(v['device_ms'] / total, 3),
             'launches_per_step': v['launches'], 'threads_per_launch': round(v['threads'] / max(v['launches'], 1))}
        if c and t > 0:
            f64 = (c['f64_add_sub_mul'] + c['f64_fma'] + c['f64_setp'] + c['f64_cvt'] + c['f64_other'] + c['f64_rcp']
                   + 10 * (c['f64_div'] + c['f64_sqrt']))          # IEEE div/sqrt: ~10 FP64-pipe instructions each
            items = 8 * n * v['launches'] if name in NODE_WALKS else v['threads']
            r.update({'model': 'node walk, 8 records per carrier' if name in NODE_WALKS else 'one item per thread',
                      'static_f64_instr': f64, 'static_int_instr': c['int_ops'],
                      'static_bytes': c['ld_global_bytes'] + c['st_global_bytes'],
                      'fp64_fraction_of_peak_upper': round(items * f64 / t / FP64_PEAK, 4),
                      'int_fraction_of_peak_upper': round(items * c['int_ops'] / t / INT_PEAK, 4),
                      'dram_fraction_of_peak_upper': round(items * (c['ld_global_bytes'] + c['st_global_bytes']) / t / BW_PEAK, 4),
                      'thread_slots_used': round(r['threads_per_launch'] / SLOTS, 3)})
        rows.append(r)
    out[str(n)] = {'device_ms_per_step': round(total, 2), 'top': rows}
Path(sys.argv[2]).write_text(json.dumps(out, indent=1) + '\n')
for n, o in out.items():
    print(f'== {n} carriers, 200 um, lane: device {o["device_ms_per_step"]} ms/step')
    for r in o['top']:
        print('  %-26s %7.2f ms %5.1f%%  fp64<=%6.2f%%  int<=%6.2f%%  dram<=%6.2f%%  slots %5.2f' % (
            r['kernel'], r['ms_per_step'], 100 * r['share'], 100 * r.get('fp64_fraction_of_peak_upper', float('nan')),
            100 * r.get('int_fraction_of_peak_upper', float('nan')), 100 * r.get('dram_fraction_of_peak_upper', float('nan')),
            r.get('thread_slots_used', float('nan'))))
