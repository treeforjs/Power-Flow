#!/usr/bin/python -I
"""Forwarding stubs for every exported function of the real driver library except the intercepted ones.
Usage: gen_stubs.py SYMBOL_LIST OUT_S OUT_H"""
import sys
skip = {'cuLaunchKernel', 'cuStreamSynchronize', 'cuCtxSynchronize'}
names = [n.strip() for n in open(sys.argv[1]) if n.strip() and n.strip() not in skip]
with open(sys.argv[2], 'w') as s:
    s.write('    .text\n')
    for n in names:
        s.write(f'    .globl {n}\n    .type {n},@function\n{n}:\n    jmp *kt_real_{n}(%rip)\n')
    s.write('    .section .note.GNU-stack,"",@progbits\n')
with open(sys.argv[3], 'w') as h:
    h.write('#define KT_FORWARDED(X) \\\n')
    h.write(' \\\n'.join(f'    X({n})' for n in names) + '\n')
