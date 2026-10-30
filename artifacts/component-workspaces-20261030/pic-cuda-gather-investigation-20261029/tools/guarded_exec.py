#!/usr/bin/env python3
"""One host-guard check then exec a C++ diagnostic/benchmark; no Python step loop."""
import json
import os
from pathlib import Path
import sys
from production_hold import paused_production

command=sys.argv[1:]
if command[:1]==['--']:
    command=command[1:]
if not command or os.sched_getaffinity(0)!={7} or os.getpriority(os.PRIO_PROCESS,0)<15:
    raise RuntimeError('Parent C++ GPU diagnostic requires a command and CPU7/nice15')
checked=paused_production('before_cpp_gpu_exec')
if receipt:=os.environ.get('PIC_CPP_PREFLIGHT_RECEIPT'):
    Path(receipt).write_text(json.dumps({'passed':True,'production_signalled':False,'hold':checked},indent=2)+'\n')
os.execvpe(command[0],command,os.environ)
