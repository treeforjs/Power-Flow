#!/usr/bin/python -I
from pathlib import Path
import runpy
P=Path(__file__).resolve().parent
m=runpy.run_path(str(P/'run_source_v7.py'),run_name='source_qualification')
for name,args in [('build',('c1','normal')),('build',('u1','ubsan')),('run',('c1',['u1']))]:
 result=m[name](*args)
 if result:raise SystemExit(result)
