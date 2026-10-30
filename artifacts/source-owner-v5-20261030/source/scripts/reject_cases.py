#!/usr/bin/env python3
"""Run compile-reject cases 0..11 (-fsyntax-only, -DREJECT_CASE=N) with the exact flags
of the source_v5_port_normal target from compile_commands.json, two at a time.
usage: reject_cases.py <build-dir> <compile_reject.cpp> <out-dir>"""
import json, shlex, subprocess, sys
from pathlib import Path
build, source, out = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
out.mkdir(parents=True, exist_ok=True)
cc = json.loads((build / 'compile_commands.json').read_text())
entry = next(e for e in cc if e['file'].endswith('audit-port/source_v5_port.cpp') and 'source_v5_port_normal.dir' in e['command'])
argv = shlex.split(entry['command'])
flags, skip = [], False
for a in argv[1:]:
    if skip:
        skip = False
        continue
    if a in ('-o', '-MF', '-MT'):
        skip = True
        continue
    if a in ('-c', '-MD') or a == entry['file']:
        continue
    flags.append(a)
(out / 'flags.json').write_text(json.dumps({'compiler': argv[0], 'flags': flags}, indent=1) + '\n')
for pair in ((0, 1), (2, 3), (4, 5), (6, 7), (8, 9), (10, 11)):
    procs = []
    for case in pair:
        cmd = [argv[0], *flags, '-fsyntax-only', f'-DREJECT_CASE={case}', source]
        log = open(out / f'case-{case:02d}.log', 'w')
        log.write('$ ' + shlex.join(cmd) + '\n')
        log.flush()
        procs.append((case, subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT), log))
    for case, p, log in procs:
        code = p.wait()
        log.close()
        (out / f'case-{case:02d}.exit').write_text(f'{code}\n')
print('compile-reject cases done')
