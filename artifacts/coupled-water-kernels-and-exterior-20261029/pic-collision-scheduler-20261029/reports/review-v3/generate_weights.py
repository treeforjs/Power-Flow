from pathlib import Path
import math,random,csv,sys
p=Path(__file__).parent
pairs=[(1e20,1),(16,8),(1,math.nextafter(1,0)),(sys.float_info.max,math.nextafter(sys.float_info.max,0)),(sys.float_info.max,math.ulp(0)),(1,math.ulp(0)),(math.ulp(0),math.ulp(0)),(2*math.ulp(0),math.ulp(0)),(sys.float_info.min,math.ulp(0)),(1e6,.031334406939830223)]
r=random.Random(4726)
for i in range(100):
 a=math.ldexp(r.uniform(1,2),r.randint(-1022,1022));b=a*r.random() if i%2 else math.ldexp(r.uniform(1,2),r.randint(-1074,1022));a,b=sorted((a,b),reverse=True);pairs.append((a,b))
with (p/'weights.csv').open('w',newline='') as f:
 w=csv.writer(f)
 for i,(a,b) in enumerate(pairs):w.writerow([f'w{i}',float(a).hex(),float(b).hex()])
print(len(pairs),'weight differences')
