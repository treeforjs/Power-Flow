import subprocess,sys
subprocess.run([sys.executable,"build-components.py","normal","v83-export","v8-sinks","audit-target-guard","audit-parcels-codec","v8-chemistry","flows","work-gross","v8-stage-b","v8-restore","step-timing","v8-flows-gpu","v8-chemistry-gpu"],check=True)
subprocess.run([sys.executable,"check-export.py"],check=True)
subprocess.run([sys.executable,"audit-tests.py","cpu","first-cpu"],check=True)
