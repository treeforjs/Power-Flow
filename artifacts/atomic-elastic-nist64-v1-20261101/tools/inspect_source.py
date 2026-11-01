"""Bounded inspection of the already recovered primary atomic DCS exports."""
from pathlib import Path
import csv,hashlib,io,json,math,os,re,sys
R=Path(__file__).resolve().parents[1]
def need(x,s):
 if not x:raise RuntimeError(s)
need(os.environ.get('SIMQ_JOB_ID') and os.environ.get('SIMQ_LAUNCH_TOKEN') and os.environ.get('SIMQ_PARTITION')=='powerflow','SimQueue powerflow required')
need(len(os.sched_getaffinity(0))==1 and os.environ['SIMQ_UNIT'] in Path('/proc/self/cgroup').read_text(),'allocated CPU/cgroup')
need(os.environ.get('CUDA_VISIBLE_DEVICES')=='','CPU only')
need(all(os.environ.get(k)=='1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']),'one thread')
need(hashlib.sha256(Path('/usr/bin/python3').resolve().read_bytes()).hexdigest()=='375532f83be620283ab780c14ca6bcc5907dc379d2b621de0c5d0a112cf5458b','frozen Python')
def read(p,cap):
 with open(p,'rb') as f:b=f.read(cap+1)
 need(len(b)<=cap,'bounded read');return b
raw=read(R/'source-pins.json',32768)
need(hashlib.sha256(raw).hexdigest()=='85deb17fd131a99b6a2472f922852ad9ede15e8f093a0606a5ba92b583c3f11d','input manifest')
d=json.loads(raw);rows=[]
for name,pin in sorted(d['files'].items()):
 b=read(Path(d['source_dir'])/name,32768);need(hashlib.sha256(b).hexdigest()==pin['sha256'],'source identity '+name)
 t=b.decode('utf-8-sig').replace('\r\n','\n');need('RELATIVISTIC\nELASTIC SCATTERING CROSS SECTION\nd sigma/d omega vs theta\n' in t,'source process/frame')
 need(re.search(r'^Atomic number: '+str(pin['atomic_number'])+r'$',t,re.M),'species');need(re.search(r'^Energy: '+str(pin['energy_eV'])+r' eV$',t,re.M),'energy')
 need('(a0**2 = 2.8002852E-21 m**2)' in t,'source units')
 total=float(re.search(r'^Total cross section: ([0-9.E+-]+) a0\*\*2$',t,re.M).group(1))
 data=list(csv.reader(io.StringIO(t.split('Theta (deg), d sigma/d omega (a0**2/sr)\n',1)[1])));data=[r for r in data if r]
 need(len(data)==181 and all(len(x)==2 for x in data),'angular cardinality')
 need([float(x[0]) for x in data]==list(range(181)),'full source angular grid')
 cs=[float(x[1]) for x in data];need(all(math.isfinite(x) and x>0 for x in cs),'DCS positive finite')
 mu=[math.cos(math.radians(a)) for a in reversed(range(181))];cs.reverse()
 integral=math.fsum((b-a)*(x+y)/2 for a,b,x,y in zip(mu,mu[1:],cs,cs[1:]))*2*math.pi
 mt=math.fsum((b-a)/6*((1-a)*x+4*(1-(a+b)/2)*(x+y)/2+(1-b)*y) for a,b,x,y in zip(mu,mu[1:],cs,cs[1:]))*2*math.pi
 rows.append({'file':name,'energy_eV':pin['energy_eV'],'Z':pin['atomic_number'],'reported_integral_a0_squared':total,'piecewise_mu_integral_a0_squared':integral,'relative_integral_difference':integral/total-1,'piecewise_mu_transport_a0_squared':mt,'mean_mu':1-mt/integral,'raw_sha256':pin['sha256']})
out=R/'evidence'/sys.argv[1];need(not out.exists(),'preserve evidence');out.mkdir()
result={'passed':True,'tables':rows,'source_pins_sha256':hashlib.sha256(raw).hexdigest(),'frozen_python_sha256':hashlib.sha256(Path('/usr/bin/python3').resolve().read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'job':os.environ['SIMQ_JOB_ID'],'cpus':sorted(os.sched_getaffinity(0)),'model_qualified':False,'production_qualified':False}
(out/'result.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'tables':len(rows),'integral_relative_difference_min':min(x['relative_integral_difference'] for x in rows),'integral_relative_difference_max':max(x['relative_integral_difference'] for x in rows)}))
