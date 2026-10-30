from pathlib import Path
import argparse,json,os,time
p=argparse.ArgumentParser();p.add_argument('--memory-bytes',type=int,required=True);p.add_argument('--swap-bytes',type=int,required=True);p.add_argument('--unit',required=True);p.add_argument('--receipt',required=True);p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args()
command=a.command[1:] if a.command[:1]==['--'] else a.command
if not command:raise RuntimeError('No scoped command')
group=Path('/proc/self/cgroup').read_text().split(':',2)[2].strip()
if Path(group).name!=a.unit:raise RuntimeError('Command is not in the dedicated job scope')
cg=Path('/sys/fs/cgroup')/group.lstrip('/')
limits={key:(cg/key).read_text().strip() for key in ['memory.max','memory.swap.max','memory.oom.group']}
if limits!={'memory.max':str(a.memory_bytes),'memory.swap.max':str(a.swap_bytes),'memory.oom.group':'1'}:raise RuntimeError('Memory scope limits are not active')
cpus=sorted(os.sched_getaffinity(0))
if cpus!=list(range(4)) or os.environ.get('SIMQ_CPU_COUNT')!='4':raise RuntimeError('Expected four-core queue allocation')
receipt={'verified_unix_s':time.time(),'pid':os.getpid(),'cgroup':group,'cpu_affinity':cpus,'limits':limits,'command':command}
path=Path(a.receipt);tmp=path.with_suffix('.json.tmp');tmp.write_text(json.dumps(receipt,indent=2)+'\n');tmp.replace(path)
os.execvpe(command[0],command,os.environ)
