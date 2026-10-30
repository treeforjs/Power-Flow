#!/usr/bin/env python3
"""Parent-controlled isolated benchmark. Never imports or initializes CUDA by default.

Requires a compiled PTX and fixture plus --parent-authorized. Caller must own the
production hold/resume policy; this script never signals or launches production.
"""
import argparse
from contextlib import ExitStack
import ctypes as c
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
from production_hold import paused_production
import statistics as stats
import struct
import subprocess
import time

class Grid(c.Structure):
    _fields_ = [('cells', c.c_uint64), ('nx', c.c_int), ('ny', c.c_int)]
class Tables(c.Structure):
    _fields_ = [('self', Grid), ('native', Grid)] + [(x, c.c_uint64) for x in (
        'phi', 'psi', 'az0', 'az1', 'rate0', 'rate1', 'xs', 'ys')] + [(x, c.c_double) for x in (
        'width', 'height', 'depth', 'short_z', 'fraction', 'voltage', 'dt')]
assert c.sizeof(Tables) == 152

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--fixture', type=Path, required=True)
p.add_argument('--ptx', type=Path)
p.add_argument('--out', type=Path, required=True)
p.add_argument('--parent-authorized', action='store_true')
p.add_argument('--validate-only', action='store_true', help='Validate input schema on CPU, without loading CUDA')
p.add_argument('--repeats', type=int, default=20)
p.add_argument('--warmup', type=int, default=3)
a = p.parse_args()
if not a.parent_authorized and not a.validate_only:
    p.error('GPU execution belongs to the parent; explicit --parent-authorized is required.')
if not a.validate_only and (os.sched_getaffinity(0) != {7} or os.getpriority(os.PRIO_PROCESS, 0) < 15):
    p.error('Parent GPU allocation requires CPU7 at nice15 or lower priority.')
if not 1 <= a.repeats <= 40 or not 0 <= a.warmup <= 5:
    p.error('Bounded request supports at most40 measured and5 warmup iterations.')
manifest = json.loads((a.fixture / 'manifest.json').read_text())
data = {x: (a.fixture / f'{x}.bin').read_bytes() for x in (
    'self_cells', 'native_cells', 'phi', 'psi', 'az0', 'az1', 'rate0', 'rate1', 'xs', 'ys', 'inputs')}
n = manifest['particles']
if not 0 < n <= 2000000 or len(data['inputs']) != n * 56:
    p.error('Invalid particle fixture length')
for x in ('self', 'native'):
    if not 1 <= manifest[f'{x}_nx'] <= 4096 or not 1 <= manifest[f'{x}_ny'] <= 4096:
        p.error('Grid dimensions exceed prototype bounds')
    if len(data[f'{x}_cells']) != manifest[f'{x}_nx'] * manifest[f'{x}_ny'] * 80:
        p.error('Invalid grid fixture length')
    field_length = len(data['phi' if x == 'self' else 'az0']) // 8
    for cell in struct.iter_unpack('<4i8d', data[f'{x}_cells']):
        if any(not 0 <= k < field_length for k in cell[:4]) or not all(map(math.isfinite, cell[4:])):
            p.error('Invalid connectivity or geometry')
for axis, dim in (('xs', 'native_nx'), ('ys', 'native_ny')):
    if len(data[axis]) != (manifest[dim] + 1) * 8:
        p.error('Invalid native axes length')
    axes = [x[0] for x in struct.iter_unpack('<d', data[axis])]
    if not all(map(math.isfinite, axes)) or any(y <= x for x, y in zip(axes, axes[1:])):
        p.error('Invalid native axis order')
if len(data['phi']) != len(data['psi']) or any(len(data[x]) != len(data['az0']) for x in ('az1', 'rate0', 'rate1')):
    p.error('Inconsistent nodal field lengths')
gpu_bytes = sum(map(len, data.values())) + n * 80
if gpu_bytes > 128 * 1024**2:
    p.error('Prototype memory budget128MiB exceeded')
for name in ('width', 'height', 'depth', 'dt'):
    if not math.isfinite(manifest[name]) or manifest[name] <= 0:
        p.error('Invalid positive geometry/timestep')
for name in ('short_z', 'fraction', 'voltage'):
    if not math.isfinite(manifest[name]):
        p.error('Nonfinite field configuration')
if not 0 <= manifest['fraction'] <= 1:
    p.error('Temporal extrapolation prohibited')
for name in ('phi', 'psi', 'az0', 'az1', 'rate0', 'rate1'):
    if len(data[name]) % 8 or any(not math.isfinite(x[0]) for x in struct.iter_unpack('<d', data[name])):
        p.error('Nonfinite nodal field')
for row in struct.iter_unpack('<7d', data['inputs']):
    if not all(map(math.isfinite, row)) or not math.isfinite(sum(x*x for x in row[3:6])):
        p.error('Invalid particle input')
    if not 0 <= row[0] <= manifest['width'] or not 0 <= row[1] <= manifest['height']:
        p.error('Particle outside PIC box')
if a.validate_only:
    print(json.dumps({'passed': True, 'gpu_executed': False, 'particles': n, 'device_allocation_bytes': gpu_bytes}))
    raise SystemExit(0)
if not a.ptx:
    p.error('Parent benchmark requires compiled --ptx')
ptx = a.ptx.read_bytes()
if not ptx or b'gather_push_kernel' not in ptx:
    p.error('PTX has no expected entry point')
build = json.loads(a.ptx.with_suffix('.build.json').read_text())
source = Path(__file__).resolve().parents[1] / 'source'
if build['core_sha256'] != hashlib.sha256((source / 'gather_push_core.hpp').read_bytes()).hexdigest() or build['kernel_sha256'] != hashlib.sha256((source / 'gather_push.cu').read_bytes()).hexdigest() or build['ptx_sha256'] != hashlib.sha256(ptx).hexdigest():
    p.error('Compiled PTX is stale or changed; rebuild before benchmarking')
result_root = Path('/home/treefor/sim-workspace/powerflow/outputs')
profile_root = result_root / 'pic-profile'

with ExitStack() as leases:
    # Same lease order as the parent replay harness: GPU, then physical test core.
    for name in ('testing-gpu.lock', 'testing-core7.lock'):
        lease = leases.enter_context((profile_root / name).open('a+'))
        fcntl.flock(lease, fcntl.LOCK_EX)
    hold_before = paused_production('before_cuda_load')
    # This is the first line that loads the CUDA driver, after every CPU gate.
    lib = c.CDLL('libcuda.so.1')
    def api(name, args):
        f = getattr(lib, name)
        f.argtypes, f.restype = args, c.c_int
        def call(*values):
            status = f(*values)
            if status:
                raise RuntimeError(f'{name}: CUDA driver status {status}')
        return call
    ptr, u64, size = c.c_void_p, c.c_uint64, c.c_size_t
    init = api('cuInit', [c.c_uint])
    device_get = api('cuDeviceGet', [c.POINTER(c.c_int), c.c_int])
    create = api('cuCtxCreate_v2', [c.POINTER(ptr), c.c_uint, c.c_int])
    destroy = api('cuCtxDestroy_v2', [ptr])
    alloc = api('cuMemAlloc_v2', [c.POINTER(u64), size])
    free = api('cuMemFree_v2', [u64])
    host_alloc = api('cuMemHostAlloc', [c.POINTER(ptr), size, c.c_uint])
    host_free = api('cuMemFreeHost', [ptr])
    h2d = api('cuMemcpyHtoD_v2', [u64, ptr, size])
    h2da = api('cuMemcpyHtoDAsync_v2', [u64, ptr, size, ptr])
    d2ha = api('cuMemcpyDtoHAsync_v2', [ptr, u64, size, ptr])
    stream_create = api('cuStreamCreate', [c.POINTER(ptr), c.c_uint])
    stream_destroy = api('cuStreamDestroy_v2', [ptr])
    module_load = api('cuModuleLoadDataEx', [c.POINTER(ptr), ptr, c.c_uint, ptr, ptr])
    module_unload = api('cuModuleUnload', [ptr])
    function_get = api('cuModuleGetFunction', [c.POINTER(ptr), ptr, c.c_char_p])
    launch = api('cuLaunchKernel', [ptr] + [c.c_uint]*7 + [ptr, c.POINTER(ptr), ptr])
    event_create = api('cuEventCreate', [c.POINTER(ptr), c.c_uint])
    event_record = api('cuEventRecord', [ptr, ptr])
    event_sync = api('cuEventSynchronize', [ptr])
    event_elapsed = api('cuEventElapsedTime', [c.POINTER(c.c_float), ptr, ptr])
    event_destroy = api('cuEventDestroy_v2', [ptr])
    context, stream, module, kernel = ptr(), ptr(), ptr(), ptr()
    device, buffers, pinned, events = c.c_int(), {}, {}, []
    started = time.perf_counter()
    try:
        init(0)
        device_get(c.byref(device), 0)
        create(c.byref(context), 4, device)  # CU_CTX_SCHED_BLOCKING_SYNC
        stream_create(c.byref(stream), 1)
        module_load(c.byref(module), c.cast(c.c_char_p(ptx), ptr), 0, None, None)
        function_get(c.byref(kernel), module, b'gather_push_kernel')
        for name, blob in data.items():
            address = u64()
            alloc(c.byref(address), len(blob))
            buffers[name] = address
            h2d(address, c.cast(c.c_char_p(blob), ptr), len(blob))
        buffers['outputs'] = u64()
        alloc(c.byref(buffers['outputs']), n * 80)
        for name, length in (('inputs', len(data['inputs'])), ('phi', len(data['phi'])), ('outputs', n * 80)):
            address = ptr()
            host_alloc(c.byref(address), length, 0)
            pinned[name] = address
        tables = Tables(Grid(buffers['self_cells'].value, manifest['self_nx'], manifest['self_ny']),
                        Grid(buffers['native_cells'].value, manifest['native_nx'], manifest['native_ny']),
                        *[buffers[x].value for x in ('phi', 'psi', 'az0', 'az1', 'rate0', 'rate1', 'xs', 'ys')],
                        *[manifest[x] for x in ('width', 'height', 'depth', 'short_z', 'fraction', 'voltage', 'dt')])
        count = c.c_uint(n)
        arguments = (ptr * 4)(c.addressof(tables), c.addressof(buffers['inputs']), c.addressof(buffers['outputs']), c.addressof(count))
        for _ in range(4):
            event = ptr()
            event_create(c.byref(event), 0)
            events.append(event)
        setup_ms = (time.perf_counter() - started) * 1000
        samples, first_output = [], None
        loop_start = time.perf_counter()
        for repetition in range(a.warmup + a.repeats):
            if time.perf_counter() - loop_start > 5:
                raise RuntimeError('Five-second steady benchmark budget exceeded')
            before = time.perf_counter()
            for name in ('inputs', 'phi'):
                c.memmove(pinned[name], data[name], len(data[name]))
            staged = time.perf_counter()
            event_record(events[0], stream)
            for name in ('inputs', 'phi'):
                h2da(buffers[name], pinned[name], len(data[name]), stream)
            event_record(events[1], stream)
            launch(kernel, (n+127)//128, 1, 1, 128, 1, 1, 0, stream, arguments, None)
            event_record(events[2], stream)
            d2ha(pinned['outputs'], buffers['outputs'], n*80, stream)
            event_record(events[3], stream)
            event_sync(events[3])
            done = time.perf_counter()
            parts = []
            for left, right in zip(events, events[1:]):
                elapsed = c.c_float()
                event_elapsed(c.byref(elapsed), left, right)
                parts.append(elapsed.value)
            blob = c.string_at(pinned['outputs'], n*80)
            if first_output is None:
                first_output = blob
            elif blob != first_output:
                raise RuntimeError('Read-only repeated device result changed')
            if repetition >= a.warmup:
                samples.append({'host_stage_ms': (staged-before)*1000, 'h2d_event_ms': parts[0],
                    'kernel_event_ms': parts[1], 'd2h_event_ms': parts[2], 'total_stage_through_sync_ms': (done-before)*1000})
        a.out.mkdir(parents=True, exist_ok=True)
        output_path = a.out / 'device.bin'
        output_path.write_bytes(first_output)
        checker = Path(__file__).resolve().parents[1] / 'build/gather_push_fixture'
        checked = subprocess.run([str(checker), 'check', str(a.fixture), str(output_path)], capture_output=True, text=True)
        (a.out / 'comparison.json').write_text(checked.stdout)
        if checked.returncode:
            raise RuntimeError(f'CPU-reference gate failed: {checked.stdout} {checked.stderr}')
        summary = {key: {'median': stats.median(s[key] for s in samples),
                         'min': min(s[key] for s in samples), 'max': max(s[key] for s in samples)} for key in samples[0]}
        (a.out / 'timing.json').write_text(json.dumps({'gpu_executed': True, 'particles': n,
            'setup_ms': setup_ms, 'gpu_allocated_bytes': gpu_bytes,
            'per_step_h2d_bytes': len(data['inputs']) + len(data['phi']), 'per_step_d2h_bytes': n*80,
            'ptx_sha256': hashlib.sha256(ptx).hexdigest(), 'samples': samples, 'summary_ms': summary,
            'logical_cpu': sorted(os.sched_getaffinity(0)), 'nice': os.getpriority(os.PRIO_PROCESS, 0),
            'production_hold_before': hold_before,
            'limitations': 'Python driver harness overhead included; particle packing, CPU output commit, field snapshot changes and full PIC advance are NOT included. Synthetic self phi/psi fixture.'}, indent=2)+'\n')
        print(json.dumps(summary, indent=2))
    finally:
        for event in events:
            event_destroy(event)
        for address in pinned.values():
            host_free(address)
        for address in buffers.values():
            if address.value:
                free(address)
        if module.value:
            module_unload(module)
        if stream.value:
            stream_destroy(stream)
        if context.value:
            destroy(context)
        hold_after = paused_production('after_cuda_cleanup')
        a.out.mkdir(parents=True, exist_ok=True)
        (a.out / 'production-hold.json').write_text(json.dumps({'before': hold_before, 'after': hold_after, 'production_signalled': False}, indent=2)+'\n')
