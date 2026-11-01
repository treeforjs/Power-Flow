"""Timing tables of native v8 (host and device modes) next to the v7.1 baseline.
usage: python3 scripts/timing_table.py HOSTDIR DEVICEDIR"""
import json, sys, glob, os
host, dev = sys.argv[1], sys.argv[2]
# v7.1 baseline (coordinator v7.1 profile, HOST-WORK-REMAINING.txt): v2 chemistry subflow and native, ms per step
BASE = {(45, 10000): (19, 11), (180, 10000): (36, 25), (45, 100000): (297, 116), (180, 100000): (177, 112),
        (45, 1000000): (15895, 1219), (180, 1000000): (2639, 1169)}
def load(d, mode, mesh, n, w):
    f = os.path.join(d, 'v8-%s-%d-%d-w%s.json' % (mode, mesh, n, w))
    return json.load(open(f)) if os.path.exists(f) else None
print('%-6s %-8s %-6s %12s %12s %12s %10s %10s %8s %8s' % ('mesh', 'carriers', 'water', 'v7.1 chem+nat', 'v8 host', 'v8 device', 'dev wall', 'dev stage', 'events', 'rssMiB'))
for w in ['1e14', '1e6']:
    for mesh in [45, 180]:
        for n in [10000, 100000, 1000000]:
            h, d = load(host, 'host', mesh, n, w), load(dev, 'device', mesh, n, w)
            b = BASE[(mesh, n)]
            hv = '%.2f' % h['native_cpu_ms_per_step_median'] if h else '-'
            dv = '%.3f' % d['native_cpu_ms_per_step_median'] if d else '-'
            dw = '%.1f' % d['native_wall_ms_per_step_median'] if d else '-'
            ds = '-'
            ev = '-'
            rss = '-'
            if d:
                last = d['per_step'][1:] or d['per_step']
                ds = '%.1f' % (sorted(s['stage1']['device_ms'] + s['stage2']['device_ms'] for s in last)[len(last) // 2])
                ev = str(sorted(s['stage1']['events'] + s['stage2']['events'] for s in last)[len(last) // 2])
                rss = str(d['max_rss_kib'] // 1024)
            print('%-6s %-8d %-6s %12s %12s %12s %10s %10s %8s %8s' % ('200um' if mesh == 45 else '50um', n, w, '%d+%d' % b, hv, dv, dw, ds, ev, rss))
