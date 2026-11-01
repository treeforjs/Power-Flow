"""Independent libm/Decimal controls; never imports the component implementation."""
from pathlib import Path
from decimal import Decimal,localcontext,Context,ROUND_HALF_EVEN
import bisect,csv,hashlib,io,json,math,sys
ROOT=Path(__file__).resolve().parents[1]
CTX=Context(prec=60,rounding=ROUND_HALF_EVEN,Emin=-999999,Emax=999999)
D=Decimal.from_float
C=Decimal(299792458)
def need(x,s):
 if not x:raise RuntimeError(s)
def close(a,b,tol=1e-11):return abs(a-b)<=tol*max(abs(a),abs(b),1e-300)
def vec(x):return [D(v) for v in x]
def dot(a,b):return sum((x*y for x,y in zip(a,b)),Decimal(0))
def momentum(m,u):return [m*x/C for x in u]
def kinetic(m,p):
 p2=dot(p,p);return p2/((m*m+p2).sqrt()+m)
def main():
 records=[json.loads(s) for s in Path(sys.argv[1]).read_text().splitlines()]
 identity=records[0];need(identity['type']=='identity','identity first')
 inspection=json.loads((ROOT/'evidence/inspect-c1/result.json').read_text())
 rows={(r['Z'],r['energy_eV']):r for r in inspection['tables']}
 pins=json.loads((ROOT/'source-pins.json').read_text());dcs={}
 for name,pin in pins['files'].items():
  raw=(Path(pins['source_dir'])/name).read_bytes();need(hashlib.sha256(raw).hexdigest()==pin['sha256'],'primary SHA')
  data=list(csv.reader(io.StringIO(raw.decode().split('Theta (deg), d sigma/d omega (a0**2/sr)\r\n',1)[1])))
  mu=[math.cos(math.radians(i)) for i in reversed(range(181))];sigma=[float(x[1]) for x in reversed(data)]
  accum=[0.]
  for a,b,x,y in zip(mu,mu[1:],sigma,sigma[1:]):accum.append(accum[-1]+(b-a)*(x+y)/2)
  dcs[(pin['atomic_number'],pin['energy_eV'])]=(mu,sigma,accum)
 summaries=events=analytic=0;max_energy=max_momentum=max_angle=max_recoil=0.;angles={1:[],8:[]}
 with localcontext(CTX):
  me=D(identity['electron_mass_eV'])
  for r in records:
   if r['type']=='summary':
    o=rows[(r['Z'],r['E'])];need(r['reported']==o['reported_integral_a0_squared'],'reported total unchanged')
    need(close(r['integrated'],o['piecewise_mu_integral_a0_squared']),'independent angular integral')
    need(close(r['transport'],o['piecewise_mu_transport_a0_squared']),'independent momentum integral');summaries+=1
   elif r['type']=='cross_section':
    E=r['E'];lo=min(10000,1000*int(E/1000));hi=min(10000,lo+1000);a=(E-lo)/1000
    expected=((1-a)*rows[(r['Z'],lo)]['reported_integral_a0_squared']+a*rows[(r['Z'],hi)]['reported_integral_a0_squared'])*2.8002852e-21
    need(close(r['sigma_m2'],expected,5e-15),'source total/interpolation')
   elif r['type'] in ('event','analytic'):
    mass=[me,D(r['M'])];U0=[vec(v) for v in r['incoming']];U1=[vec(v) for v in r['outgoing']]
    P0=[momentum(m,u) for m,u in zip(mass,U0)];P1=[momentum(m,u) for m,u in zip(mass,U1)]
    K0=[kinetic(m,p) for m,p in zip(mass,P0)];K1=[kinetic(m,p) for m,p in zip(mass,P1)]
    er=abs(sum(K1)-sum(K0))/max(Decimal(1),sum(K1)+sum(K0));max_energy=max(max_energy,float(er));need(er<Decimal('2e-13'),'independent energy closure')
    residual=[P1[0][j]+P1[1][j]-P0[0][j]-P0[1][j] for j in range(3)]
    scale=sum(dot(p,p).sqrt() for p in P0+P1)
    pr=dot(residual,residual).sqrt()/max(Decimal(1),scale);max_momentum=max(max_momentum,float(pr));need(pr<Decimal('2e-13'),'independent momentum closure')
    for direct,ke in ((K0,r['K0']),(K1,r['K1'])):
     for a,b in zip(direct,ke):need(abs(a-D(b))<Decimal('3e-13')*max(Decimal(1),abs(a)),'canonical stored kinetic')
    w=[x/C for x in U0[1]];g=(1+dot(w,w)).sqrt()
    def unboost(p,m):
     energy=m+kinetic(m,p);a=dot(w,p)/(g+1)-energy;return [x+a*y for x,y in zip(p,w)]
    before=unboost(P0[0],me);after=unboost(P1[0],me)
    recovered=dot(before,after)/(dot(before,before)*dot(after,after)).sqrt()
    ae=abs(float(recovered)-r['mu']);max_angle=max(max_angle,ae);need(ae<2e-12,'target-rest laboratory angle retained')
    recoil=kinetic(mass[1],unboost(P1[1],mass[1]));re=abs(float(recoil)-r['recoil'])/max(1.,float(recoil));max_recoil=max(max_recoil,re);need(re<2e-10,'recoil target-rest consistency')
    if r['type']=='analytic':
     analytic+=1;T=D(r['E']);M=mass[1]
     if r['mu']==-1:
      expected=2*M*T*(T+2*me)/(M*M+me*me+2*M*(me+T))
      need(abs(recoil-expected)<Decimal('3e-13')*max(Decimal(1),expected),'independent maximum recoil formula')
     if r['mu']==1:need(recoil<Decimal('1e-20'),'zero forward recoil')
    else:events+=1
   elif r['type']=='angle':angles[r['Z']].append(r['mu'])
 errors=[]
 for Z,values in angles.items():
  need(len(values)==8192,'statistical population')
  mu,sigma,accum=dcs[(Z,5000)]
  for x in [-.5,0,.5,.9,.95,.99,.999]:
   i=bisect.bisect_right(mu,x)-1;dx=x-mu[i];h=mu[i+1]-mu[i]
   expected=(accum[i]+dx*(sigma[i]+.5*(sigma[i+1]-sigma[i])*dx/h))/accum[-1]
   observed=sum(v<=x for v in values)/len(values);error=abs(expected-observed);errors.append(error);need(error<.025,'angular CDF population sanity')
 need(summaries==20 and events==8192 and analytic==190,'expected complete controls')
 need(records[-1]['events']==events+analytic,'fixture complete')
 print(json.dumps({'passed':True,'summaries':summaries,'random_events':events,'analytic_events':analytic,'assertions':records[-1]['assertions'],'max_relative_energy_error':max_energy,'max_relative_momentum_error':max_momentum,'max_angle_error':max_angle,'max_recoil_scaled_error':max_recoil,'max_empirical_cdf_error':max(errors),'independent_decimal_precision':60,'production_qualified':False,'physical_model_validated':False,'native_adopted':False},sort_keys=True))
if __name__=='__main__':main()
