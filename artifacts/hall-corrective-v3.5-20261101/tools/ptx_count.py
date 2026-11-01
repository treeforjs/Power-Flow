#!/usr/bin/python -I
"""Static instruction counts of the Hall v3.4 device kernels (PTX from NVRTC with the session's options: compute_86,
--fmad=false --prec-div=true --prec-sqrt=true --ftz=false --std=c++17). Per kernel entry: FP64 arithmetic (add/sub/mul,
fma, div, sqrt, rcp, min/max/abs/neg, setp, cvt), integer, global loads/stores (count and bytes), branches. Static counts
of the PTX text; called device functions are added once per call site; loops are not unrolled into dynamic counts.
Usage: ptx_count.py NVRTC_LIBRARY OUT_JSON"""
import ctypes, json, re, sys
from pathlib import Path

R = Path('/home/treefor/projects/runtimes')
H = R / 'hall-corrective-v3.4-20261031/source'
G = R / 'hall-gpu-resident-v1-20261030'
F = R / 'pic-reaction-transactions-20261029/native-full-time-candidate007.1/forms/include/pf_v7/forms_core.h'
lib = ctypes.CDLL(sys.argv[1])
OPTS = [b'--gpu-architecture=compute_86', b'--fmad=false', b'--prec-div=true', b'--prec-sqrt=true', b'--ftz=false', b'--std=c++17']


def ptx_of(name, main, headers):
    prog = ctypes.c_void_p()
    names = [n.encode() for n in headers]
    srcs = [Path(p).read_bytes() for p in headers.values()]
    H_ = (ctypes.c_char_p * len(names))(*srcs); N_ = (ctypes.c_char_p * len(names))(*names)
    rc = lib.nvrtcCreateProgram(ctypes.byref(prog), Path(main).read_bytes(), name.encode(), len(names), H_, N_)
    assert rc == 0, rc
    O = (ctypes.c_char_p * len(OPTS))(*OPTS)
    rc = lib.nvrtcCompileProgram(prog, len(OPTS), O)
    if rc != 0:
        size = ctypes.c_size_t(); lib.nvrtcGetProgramLogSize(prog, ctypes.byref(size))
        log = ctypes.create_string_buffer(size.value); lib.nvrtcGetProgramLog(prog, log)
        raise SystemExit(name + ' compile failed:\n' + log.value.decode()[-4000:])
    size = ctypes.c_size_t(); lib.nvrtcGetPTXSize(prog, ctypes.byref(size))
    buf = ctypes.create_string_buffer(size.value); lib.nvrtcGetPTX(prog, buf)
    return buf.value.decode()


def classify(ptx):
    out, funcs, callees = {}, {}, {}
    blocks = re.split(r'\n(?=(?:\.(?:visible|weak|extern)\s+)?\.(?:entry|func)\s)', ptx)
    for blk in blocks:
        m = re.match(r'(?:\.(?:visible|weak|extern)\s+)?\.(entry|func)\s+(?:\([^)]*\)\s*)?(\w+)', blk)
        if not m or '{' not in blk:
            continue
        kind, name = m.group(1), m.group(2)
        if name in out or name in funcs:
            name = name + '#dup'
        callees[name] = re.findall(r'call(?:\.uni)?\s+(?:\([^)]*\),\s*)?(\w+)', blk)
        c = {'f64_add_sub_mul': 0, 'f64_fma': 0, 'f64_div': 0, 'f64_sqrt': 0, 'f64_rcp': 0, 'f64_other': 0, 'f64_setp': 0, 'f64_cvt': 0,
             'int_ops': 0, 'ld_global': 0, 'ld_global_bytes': 0, 'st_global': 0, 'st_global_bytes': 0, 'branches': 0, 'calls': 0, 'lines': 0}
        for line in blk.splitlines():
            t = line.strip()
            if not t or t.startswith('//') or t.startswith('.') or t.endswith(':') or t.startswith('{') or t.startswith('}'):
                continue
            c['lines'] += 1
            op = t.split()[0].rstrip(';')
            parts = op.split('.')
            base = parts[0].lstrip('@!%p0123456789 ')
            if t.startswith('@'):
                op2 = t.split()[1]; parts = op2.split('.'); base = parts[0]
            if base in ('bra', 'ret', 'exit'):
                c['branches'] += base == 'bra'; continue
            if base == 'call':
                c['calls'] += 1; continue
            if base in ('ld', 'st') and not any(sp in parts for sp in ('param', 'local', 'shared', 'const')):   # global or generic
                width = 8 if any(p in ('f64', 'u64', 's64', 'b64') for p in parts) else 4 if any(p in ('f32', 'u32', 's32', 'b32') for p in parts) else 1
                vec = 4 if 'v4' in parts else 2 if 'v2' in parts else 1
                key = 'ld_global' if base == 'ld' else 'st_global'
                c[key] += 1; c[key + '_bytes'] += width * vec; continue
            if 'f64' in parts:
                if base in ('add', 'sub', 'mul'): c['f64_add_sub_mul'] += 1
                elif base == 'fma': c['f64_fma'] += 1
                elif base == 'div': c['f64_div'] += 1
                elif base == 'sqrt': c['f64_sqrt'] += 1
                elif base == 'rcp': c['f64_rcp'] += 1
                elif base == 'setp': c['f64_setp'] += 1
                elif base == 'cvt': c['f64_cvt'] += 1
                else: c['f64_other'] += 1
            elif base in ('add', 'sub', 'mul', 'mad', 'shl', 'shr', 'and', 'or', 'xor', 'setp', 'selp', 'min', 'max', 'not', 'cvt', 'mov', 'clz', 'popc', 'bfe', 'prmt'):
                c['int_ops'] += 1
        # FP64 "flops": add/sub/mul 1, fma 2; div and sqrt with --prec-div/--prec-sqrt expand to ~20 FP64 operations each in SASS
        (out if kind == 'entry' else funcs)[name] = c
    # Called device functions (not inlined by NVRTC) counted once per call site, transitively (static, loops not unrolled).
    def total(name, depth=0):
        c = dict(out.get(name) or funcs[name])
        if depth < 8:
            for callee in callees.get(name, []):
                if callee in funcs:
                    for k, v in total(callee, depth + 1).items():
                        if k != 'f64_flop_static':
                            c[k] += v
        c['f64_flop_static'] = c['f64_add_sub_mul'] + 2 * c['f64_fma'] + 20 * (c['f64_div'] + c['f64_sqrt']) + 4 * c['f64_rcp'] + c['f64_other']
        return c
    return {name: total(name) for name in out}


programs = {
    'push_v3': (H / 'gpu/push_kernels_v3.cu', {'push_core.hpp': G / 'push/push_core.hpp', 'push_abi.hpp': G / 'push/push_abi.hpp',
                                                'push_bodies.hpp': G / 'push/push_bodies.hpp', 'sink_core.hpp': H / 'sink_core.hpp',
                                                'push_abi_v3.hpp': H / 'gpu/push_abi_v3.hpp', 'push_bodies_v3.hpp': H / 'gpu/push_bodies_v3.hpp',
                                                'deterministic_primitives.cuh': G / 'common/device/deterministic_primitives.cuh'}),
    'step_v34': (H / 'gpu/step_kernels_v34.cu', {'step_abi_v34.h': H / 'gpu/step_abi_v34.h', 'exact_sum_v34.cuh': H / 'gpu/exact_sum_v34.cuh',
                                                  'x87_v34.hpp': H / 'gpu/x87_v34.hpp', 'support_v34.hpp': H / 'gpu/support_v34.hpp',
                                                  'sink_core.hpp': H / 'sink_core.hpp', 'push_core.hpp': G / 'push/push_core.hpp', 'forms_core.h': F}),
    'resident': (H / 'gpu/resident_kernels.cu', {'resident_abi.h': H / 'gpu/resident_abi.h', 'forms_core.h': F,
                                                  'exact_sum_v34.cuh': H / 'gpu/exact_sum_v34.cuh'}),
}
res = {}
for name, (main, headers) in programs.items():
    res[name] = classify(ptx_of(name, main, headers))
# CG kernels: the kernel text is a C++ raw string in cg_kernels.hpp
src = (H / 'gpu/cg_kernels.hpp').read_text()
cg = src[src.index('R"CUDA(') + 7: src.index(')CUDA"')]
tmp = Path(sys.argv[2]).with_suffix('.cg.cu'); tmp.write_text(cg)
res['cg'] = classify(ptx_of('cg', tmp, {}))
tmp.unlink()
Path(sys.argv[2]).write_text(json.dumps(res, indent=1) + '\n')
for prog, ks in res.items():
    for k, c in ks.items():
        print('%-9s %-28s flop_static=%5d f64[asm=%d fma=%d div=%d sqrt=%d] ld=%d(%dB) st=%d(%dB) br=%d' % (
            prog, k, c['f64_flop_static'], c['f64_add_sub_mul'], c['f64_fma'], c['f64_div'], c['f64_sqrt'], c['ld_global'], c['ld_global_bytes'],
            c['st_global'], c['st_global_bytes'], c['branches']))
