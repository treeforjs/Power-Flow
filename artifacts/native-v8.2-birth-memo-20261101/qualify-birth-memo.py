"""Sequential owner/borrowable one-core functional qualification, no unleased work."""
import subprocess,sys
from tooling import *
resource_guard(build=True)
commands=[
 [sys.executable,'build-components.py','normal','core','hall','v8-sinks','step-timing','audit-target-guard',
  'audit-parcels-codec','v8-chemistry','v8-restore','flows','work-gross','v8-stage-b','v8-chemistry-gpu','v8-flows-gpu'],
 [sys.executable,'build-differential.py'],
 [sys.executable,'audit-tests.py','cpu','final-cpu'],
 [sys.executable,'run-differential.py','functional-differential',10000,8]]
if len(sys.argv)>1:
    need(sys.argv[1]=='checks','unknown qualification selector')
    commands=commands[1:]
for command in commands:
    subprocess.run([str(x) for x in command],cwd=P,check=True)
print('birth-memo functional CPU qualification passed')
