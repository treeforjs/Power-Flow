import subprocess,sys
subprocess.run([sys.executable,"build-components.py","normal","v83-export"],check=True)
subprocess.run([sys.executable,"check-export.py","sealed-export"],check=True)
subprocess.run([sys.executable,"audit-tests.py","cpu","sealed-cpu"],check=True)
