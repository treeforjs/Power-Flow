import subprocess,sys
from tooling import *
subprocess.run([sys.executable,"build-components.py","normal","v83-export","v8-restore"],check=True)
subprocess.run([sys.executable,"check-export.py","binding-export"],check=True)
subprocess.run([sys.executable,"audit-tests.py","cpu","binding-cpu", "restore-2-unbroken","restore-2-save","restore-2-restore","restore-3-unbroken","restore-3-save","restore-3-restore","restore-3-tamper"],check=True)
