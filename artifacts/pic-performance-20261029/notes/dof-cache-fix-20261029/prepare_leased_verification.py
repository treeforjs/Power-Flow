import pathlib
out=pathlib.Path(__file__).parent
p=out/"run_cpu_checks.py"
s=p.read_text()
s=s.replace("import subprocess,pathlib,os,json,hashlib,datetime,shutil,time", "import subprocess,pathlib,os,json,hashlib,datetime,shutil,time,fcntl",1)
needle="commands=[];results=[]\n"
replacement="""commands=[];results=[]
lease_path=review.parent/\"testing-core.lock\"
lease=lease_path.open(\"a+\")
print(\"Waiting for testing-core lease\",flush=True)
fcntl.flock(lease,fcntl.LOCK_EX)
print(\"Testing-core lease acquired\",flush=True)
"""
assert needle in s
s=s.replace(needle,replacement,1)
s=s.replace("\"compiler_jobs\":1,\"CUDA_VISIBLE_DEVICES\"", "\"compiler_jobs\":1,\"testing_core_lease\":str(lease_path),\"lease_held\":True,\"CUDA_VISIBLE_DEVICES\"",1)
p.write_text(s)
