import subprocess, sys
for command in ([sys.executable, "build-differential.py"], [sys.executable, "run-differential.py", "functional-differential-c2", "10000", "8"]):
    subprocess.run(command, check=True)
