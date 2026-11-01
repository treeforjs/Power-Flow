import subprocess, sys
for command in ([sys.executable, "build-differential.py"], [sys.executable, "run-differential.py", "functional-differential-c3", "10000", "8"], [sys.executable, "run-differential.py", "functional-anchor-differential-c3", "10000", "17"]):
    subprocess.run(command, check=True)
