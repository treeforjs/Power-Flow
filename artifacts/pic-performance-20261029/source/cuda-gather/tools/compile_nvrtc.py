#!/usr/bin/env python3
"""CPU-only CUDA C++ -> PTX. Requires an explicitly supplied, coordinated NVRTC.

No package download/install or CUDA driver/runtime initialization occurs here.
"""
import argparse
import ctypes as c
import hashlib
import json
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--nvrtc', type=Path, required=True)
p.add_argument('--out', type=Path, required=True)
a = p.parse_args()
if not a.nvrtc.is_file():
    p.error('NVRTC is absent. See toolchain-proposal.json; coordinate with parent before downloading.')
lib = c.CDLL(str(a.nvrtc.resolve()))
lib.nvrtcCreateProgram.argtypes = [c.POINTER(c.c_void_p), c.c_char_p, c.c_char_p, c.c_int, c.POINTER(c.c_char_p), c.POINTER(c.c_char_p)]
lib.nvrtcCompileProgram.argtypes = [c.c_void_p, c.c_int, c.POINTER(c.c_char_p)]
for name in ('nvrtcGetProgramLogSize', 'nvrtcGetPTXSize'):
    getattr(lib, name).argtypes = [c.c_void_p, c.POINTER(c.c_size_t)]
for name in ('nvrtcGetProgramLog', 'nvrtcGetPTX'):
    getattr(lib, name).argtypes = [c.c_void_p, c.c_void_p]
lib.nvrtcDestroyProgram.argtypes = [c.POINTER(c.c_void_p)]
def ok(result):
    if result:
        raise RuntimeError(f'NVRTC status {result}')

source = Path(__file__).resolve().parents[1] / 'source'
kernel = (source / 'gather_push.cu').read_bytes()
core = (source / 'gather_push_core.hpp').read_bytes()
headers = (c.c_char_p * 1)(core)
names = (c.c_char_p * 1)(b'gather_push_core.hpp')
program = c.c_void_p()
ok(lib.nvrtcCreateProgram(c.byref(program), kernel, b'gather_push.cu', 1, headers, names))
options = [b'--std=c++17', b'--gpu-architecture=compute_86', b'--fmad=false',
           b'--ftz=false', b'--prec-div=true', b'--prec-sqrt=true', b'--restrict']
try:
    result = lib.nvrtcCompileProgram(program, len(options), (c.c_char_p * len(options))(*options))
    size = c.c_size_t()
    ok(lib.nvrtcGetProgramLogSize(program, c.byref(size)))
    log = c.create_string_buffer(size.value)
    ok(lib.nvrtcGetProgramLog(program, log))
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.with_suffix('.nvrtc.log').write_bytes(log.value)
    ok(result)
    ok(lib.nvrtcGetPTXSize(program, c.byref(size)))
    output = c.create_string_buffer(size.value)
    ok(lib.nvrtcGetPTX(program, output))
    a.out.write_bytes(output.value)
    a.out.with_suffix('.build.json').write_text(json.dumps({
        'gpu_executed': False, 'options': [x.decode() for x in options],
        'core_sha256': hashlib.sha256(core).hexdigest(),
        'kernel_sha256': hashlib.sha256(kernel).hexdigest(),
        'ptx_sha256': hashlib.sha256(output.value).hexdigest(),
        'nvrtc_library': str(a.nvrtc.resolve()),
    }, indent=2) + '\n')
    print(f'PTX written without GPU execution: {a.out}')
finally:
    lib.nvrtcDestroyProgram(c.byref(program))
