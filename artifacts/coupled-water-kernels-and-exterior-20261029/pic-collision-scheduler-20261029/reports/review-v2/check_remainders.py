from pathlib import Path
import json,csv,decimal,os
p=Path(__file__).parent
assert sorted(os.sched_getaffinity(0))==[7] and os.getpriority(os.PRIO_PROCESS,0)>=15
with decimal.localcontext() as ctx:
 ctx.prec=900;r=[]
 for x in csv.DictReader((p/'remainders.csv').open()):
  a,b,h,l=[decimal.Decimal.from_float(float(x[k])) for k in ['a','b','high','low']]
  r.append({'label':x['label'],'exact':a-b==h+l,'low_nonzero':l!=0,'single_conversion_rejected':bool(int(x['single_conversion_rejected']))})
 assert all(x['exact'] and x['low_nonzero']==x['single_conversion_rejected'] for x in r)
 out={'status':'passed','oracle':'900-decimal-digit exact input/output binary64 values','cases':len(r),'nonzero_low_conversion_rejections':sum(x['low_nonzero'] for x in r),'results':r,'cpu_affinity':sorted(os.sched_getaffinity(0)),'nice':os.getpriority(os.PRIO_PROCESS,0)}
 (p/'remainder-oracle.json').write_text(json.dumps(out,indent=2)+'\n');print('exact differences',len(r),'nonzero low conversion rejections',out['nonzero_low_conversion_rejections'])
