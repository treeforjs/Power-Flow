#!/usr/bin/python -I
"""Timing table (host and device per coupled step, v3.3 paths against v3.4) and the per-kernel resource split.
Host: bench tables of record (bench-tables/b3.json: owner = v3.3 host Owner, v33 = v3.3 host path of the resident owner;
b12.json: device, lane = v3.4), median ms per step over the steps after the first, host = thread CPU of prepare +
preflight + publish, wall = prepare + preflight wall. Device: kernel-timer runs (job 189: lane 200 um, run-N-S; m2/:
run-MODE-T-N-S), steady step = (5-step - 2-step) / 3, sum of kernel device times.
Split: per kernel, the use of the FP64 pipe, INT32 and DRAM implied by the static PTX counts against RTX 3080 Ti peaks
(266e9 FP64 instr/s, 8.5e12 INT32/s, 912e9 B/s; balance 0.29 FP64 instr per byte). Node walks (one thread per node over
its records in carrier order): the whole static count applied to each of the 8 records per carrier -- an upper bound;
below 50% of every peak the kernel is latency-bound. Other kernels: static count per thread -- an estimate (untaken
branches over-count, small adjacency loops under-count); the static intensity (FP64 instr / byte) gives the roofline
side: FP64 candidate above the balance, bandwidth candidate below it, latency/launch-bound when both estimated uses are
under 10%. Kernels with data-dependent per-thread loops (exact-sum chunks, edge walks) are not estimated.
Usage: device_table.py DIR OUT_JSON"""
import json, sys
from pathlib import Path

D = Path(sys.argv[1]); BT = D.parent / 'bench-tables'
FP64_PEAK, INT_PEAK, BW_PEAK, SLOTS = 266e9, 8.5e12, 912e9, 122880
NODE_WALKS = {'hgp_reduce_end_nodes_v3', 'hgp_reduce_start_nodes', 'k34_charge_nodes', 'k34_support_nodes', 'k34_census_nodes'}
LOOPED = {'k34_sum', 'k34_sum_finish', 'hgp_reduce_edges', 'k_res_sum'}
BALANCE = FP64_PEAK / BW_PEAK
counts = {k: c for prog in json.loads((D / 'ptx-counts.json').read_text()).values() for k, c in prog.items()}
host = {**json.loads((BT / 'b3.json').read_text()), **json.loads((BT / 'b12.json').read_text())}
UM = {45: '200um', 180: '50um'}


def kernels(p):
    j = json.loads(p.read_text())
    return {k['name']: k for k in j['kernels']}, j.get('undrained_launches', 0)


def steady(mode, t, n):
    if mode == 'lane' and t == 45:
        a, b = D / f'run-{n}-2/kernels.json', D / f'run-{n}-5/kernels.json'
    else:
        a, b = D / f'm2/run-{mode}-{t}-{n}-2/kernels.json', D / f'm2/run-{mode}-{t}-{n}-5/kernels.json'
    if not (a.exists() and b.exists()):
        return None
    (ka, ua), (kb, ub) = kernels(a), kernels(b)
    z = {'launches': 0, 'device_ms': 0.0, 'threads': 0}
    return {k: {f: (v[f] - ka.get(k, z)[f]) / 3 for f in z} for k, v in kb.items()}, ua + ub


def split(name, v, n):
    c = counts.get(name); t = v['device_ms'] * 1e-3
    r = {'kernel': name, 'ms_per_step': round(v['device_ms'], 3), 'launches_per_step': v['launches'],
         'threads_per_launch': round(v['threads'] / max(v['launches'], 1))}
    if not c or t <= 0:
        r['class'] = 'no static counts'; return r
    if name in LOOPED:
        r['class'] = 'data-dependent loop, not estimated'; return r
    f64 = (c['f64_add_sub_mul'] + c['f64_fma'] + c['f64_setp'] + c['f64_cvt'] + c['f64_other'] + c['f64_rcp']
           + 10 * (c['f64_div'] + c['f64_sqrt']))
    items = 8 * n * v['launches'] if name in NODE_WALKS else v['threads']
    u = {'fp64': items * f64 / t / FP64_PEAK, 'int32': items * c['int_ops'] / t / INT_PEAK,
         'dram': items * (c['ld_global_bytes'] + c['st_global_bytes']) / t / BW_PEAK}
    nbytes = c['ld_global_bytes'] + c['st_global_bytes']
    inten = f64 / nbytes if nbytes else float('inf')
    if name in NODE_WALKS:
        cls = 'latency-bound' if max(u.values()) < 0.5 else 'walk near a peak: ' + max(u, key=u.get)
    elif max(u['fp64'], u['dram']) < 0.10:
        cls = 'latency/launch-bound'
    else:
        cls = 'FP64 candidate' if inten > BALANCE else 'bandwidth candidate'
    r.update({'model': 'node walk (upper bound)' if name in NODE_WALKS else 'per thread (estimate)',
              **{k + '_use': round(x, 4) for k, x in u.items()}, 'static_fp64_per_byte': round(inten, 3),
              'thread_slots': round(r['threads_per_launch'] / SLOTS, 3), 'class': cls})
    return r


table, splits = [], {}
for t in (45, 180):
    for n in (100000, 1000000):
        for mode in ('owner', 'v33', 'device', 'lane'):
            h = host.get(f'{mode}/{t}/{n}', {})
            s = steady(mode, t, n)
            row = {'mesh': UM[t], 'carriers': n, 'mode': mode, 'host_cpu_ms': h.get('host_cpu_ms'),
                   'wall_ms': (h.get('prepare_wall_ms') or 0) + (h.get('preflight_wall_ms') or 0) if h else None,
                   'iterations': h.get('iterations'), 'device_ms': None}
            if s:
                ks, und = s
                row['device_ms'] = round(sum(v['device_ms'] for v in ks.values()), 2)
                row['undrained_launches'] = und
                if mode == 'lane':
                    top = sorted(ks.items(), key=lambda kv: -kv[1]['device_ms'])[:12]
                    by_class = {}
                    for k, v in ks.items():
                        cl = split(k, v, n)['class']
                        by_class[cl] = round(by_class.get(cl, 0) + v['device_ms'] / row['device_ms'], 3)
                    splits[f'{UM[t]}/{n}'] = {'device_ms': row['device_ms'], 'share_by_class': by_class,
                                              'kernels': [dict(split(k, v, n), share=round(v['device_ms'] / row['device_ms'], 3)) for k, v in top]}
            table.append(row)
Path(sys.argv[2]).write_text(json.dumps({'table': table, 'split_lane': splits}, indent=1) + '\n')
print('per coupled step, ms (host = thread CPU; wall = prepare + preflight; device = sum of kernel times, kernel timer)')
print('  %-6s %-8s %-7s %10s %10s %10s %5s' % ('mesh', 'carriers', 'mode', 'host_cpu', 'wall', 'device', 'iter'))
for r in table:
    f = lambda x: '-' if x is None else ('%.1f' % x if x >= 10 else '%.2f' % x)
    print('  %-6s %-8.0e %-7s %10s %10s %10s %5s' % (r['mesh'], r['carriers'], r['mode'], f(r['host_cpu_ms']), f(r['wall_ms']),
                                                   f(r['device_ms']), '-' if r['iterations'] is None else '%.0f' % r['iterations']))
for key, s in splits.items():
    print(f'== lane {key}: device {s["device_ms"]} ms/step; share by class: ' +
          ', '.join('%s %.0f%%' % (k, 100 * v) for k, v in sorted(s['share_by_class'].items(), key=lambda kv: -kv[1])))
    for r in s['kernels']:
        print('  %-26s %8.2f ms %5.1f%%  fp64 %7.1f%% int %5.1f%% dram %6.1f%% I %6.2f slots %6.2f  %s' % (
            r['kernel'], r['ms_per_step'], 100 * r['share'], 100 * r.get('fp64_use', float('nan')), 100 * r.get('int32_use', float('nan')),
            100 * r.get('dram_use', float('nan')), r.get('static_fp64_per_byte', float('nan')), r.get('thread_slots', float('nan')), r['class']))
