"""Record bounded upstream verification commands without shell interpretation."""
import datetime, hashlib, json, pathlib, subprocess, sys, time
root=pathlib.Path(__file__).resolve().parent
label, cwd, *cmd=sys.argv[1:]
started=datetime.datetime.now(datetime.timezone.utc).isoformat()
t0=time.monotonic()
log=root/(label+'.log')
with log.open('w') as f:
 p=subprocess.run(cmd,cwd=cwd,stdout=f,stderr=subprocess.STDOUT,text=True)
entry={'command':cmd,'cwd':cwd,'start_utc':started,'end_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'elapsed_s':time.monotonic()-t0,'exit_code':p.returncode,'log':str(log),'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest()}
with (root/'commands.jsonl').open('a') as f: f.write(json.dumps(entry)+'\n')
print(json.dumps(entry,indent=2))
print(log.read_text()[-9000:])
raise SystemExit(p.returncode)
