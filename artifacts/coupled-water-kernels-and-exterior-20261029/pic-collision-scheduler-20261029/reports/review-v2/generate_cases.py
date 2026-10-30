from pathlib import Path
import math,csv,json,random
p=Path(__file__).parent;c=299792458.;me=9.1093837139e-31;qe=1.602176634e-19
cases=[]
def add(name,u,v):
 assert all(math.isfinite(x) for x in (*u,*v));cases.append({'label':name,'u':list(u),'v':list(v)})
for i,u in enumerate([0,1e-200,1e-100,1e-4,1,1e7,1e12,1e20,1e100,1e300]):add('stationary_'+str(i),(u,0,0),(0,0,0))
for i,u in enumerate([1,1e8,1e12,1e20,1e100,1e300]):add('comoving_'+str(i),(u,.37*u,-.21*u),(u,.37*u,-.21*u))
rest=me*c*c/qe;u0=c*math.sqrt((8000/rest)*(8000/rest+2));g0=math.sqrt(1+u0*u0/(c*c))
for i,g in enumerate([1.014770,2,1e3,1e6,1e12,1e100,1e160,1e200]):
 beta=math.sqrt(1-1/(g*g)) if g<1e154 else 1.;v=beta*c
 add('boost8keV_'+str(i),(g*(u0+v*g0),0,0),(g*v,0,0))
for i,u in enumerate([1e8,1e20,1e100,1e300]):
 add('opposite_'+str(i),(u,0,0),(-u,0,0));add('skew_'+str(i),(u,.3*u,-.4*u),(-.2*u,.7*u,.6*u))
for i,u in enumerate([1e8,1e20,1e100,1e300]):
 add('parallel_ulp_'+str(i),(u,0,0),(math.nextafter(u,math.inf),0,0))
 add('skew_ulp_'+str(i),(u,.37*u,-.21*u),(math.nextafter(u,math.inf),.37*u,-.21*u))
r=random.Random(82651)
for i in range(20):
 scale=10**r.uniform(-10,120);add('random_'+str(i),tuple(r.uniform(-1,1)*scale for _ in range(3)),tuple(r.uniform(-1,1)*scale for _ in range(3)))
with (p/'cases.csv').open('w',newline='') as f:
 w=csv.writer(f)
 for x in cases:w.writerow([x['label'],*(repr(v) for v in x['u']),*(repr(v) for v in x['v'])])
(p/'cases.json').write_text(json.dumps(cases,indent=2)+'\n');print(len(cases),'oracle cases')
