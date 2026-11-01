import subprocess,sys
subprocess.run([sys.executable,"build-components.py","normal","v83-export","v8-restore","v8-chemistry-gpu"],check=True)
subprocess.run([sys.executable,"check-export.py","final-export"],check=True)
subprocess.run([sys.executable,"audit-tests.py","cpu","final-cpu"],check=True)
