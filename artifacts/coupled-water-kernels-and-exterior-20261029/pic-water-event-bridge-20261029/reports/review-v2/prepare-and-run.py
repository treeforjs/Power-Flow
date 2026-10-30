from pathlib import Path
import os,subprocess,hashlib,json
p=Path(__file__).parent;old=p.parent/'review'
assert sorted(os.sched_getaffinity(0))==[7]
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
for f,r in json.loads((old/'output-manifest.json').read_text()).items():assert sha(old/f)==r['sha256'],f
source=(old/'probe.cpp').read_text()
assert source.count('ParentBridge bridge(data,parent_model,argv[2],argv[3]);')==1
source=source.replace('ParentBridge bridge(data,parent_model,argv[2],argv[3]);','ParentBridge bridge(data,parent_model,argv[2],argv[3],2.9915e-26);')
source=source.replace('<<p.energy_defect_J;}catch','<<p.energy_defect_J<<",\\\"momentum_bound\\\":"<<p.momentum_roundoff_limit_kg_m_s<<",\\\"momentum_defect\\\":";vec(o,p.momentum_defect_kg_m_s);}catch')
anchor='  axis_probe("near_unit_axis_tiny_weight",event,e2,n2);'
assert anchor in source
source=source.replace(anchor,anchor+'''\n  for(double weight:{1e-292,1e-289,1e-284,1e-280}){e2=e;n2=n;e2.weight=weight;n2.weight=64*weight;event=pc::candidates({e2},{n2},{{1,1e-12}},{bridge.channel()})[0];std::string name="valid_weight_"+std::to_string(static_cast<int>(std::log10(weight)));mutate(name.c_str(),e2,n2,event);}\n''')
(p/'replay-probe.cpp').write_text(source)
physical=source.replace('for(double boost:{0.,.8,1e2,1e3,1e4,1e5,1e6,1e7,1e8})','for(double boost:{0.,1e-6,1e-4,1e-3,.8,1e5})')
assert physical!=source
(p/'physical-probe.cpp').write_text(physical)
oracle=(old/'oracle.py').read_text()
for name in ('replay','physical'):
 text=oracle.replace("p/'probes.jsonl'",f"p/'{name}-probes.jsonl'").replace("p/'independent-oracle.json'",f"p/'{name}-oracle.json'")
 (p/f'{name}-oracle.py').write_text(text)
(p/'build').mkdir(exist_ok=True)
for name in ('replay','physical'):
 cmd=['g++','-std=c++20','-O2','-Wall','-Wextra','-Wpedantic','-Werror','-I'+str(p),'-I'+str(p/'inputs'),'-I'+str(p/'inputs/include'),str(p/f'{name}-probe.cpp'),'-o',str(p/f'build/{name}-probe')]
 subprocess.run(cmd,check=True)
 subprocess.run([str(p/f'build/{name}-probe'),str(p/'inputs/h2o-gas.json'),sha(p/'inputs/h2o-gas.json'),sha(p/'inputs/include/raikou/h2o_conditional_parent.h'),str(p/f'{name}-probes.jsonl')],check=True)
 subprocess.run(['python3',str(p/f'{name}-oracle.py')],check=True,stdout=(p/f'{name}-oracle-summary.txt').open('w'))
for name in ('replay','physical'):
 s=json.loads((p/f'{name}-oracle.json').read_text())['summary']
 print(name,json.dumps({'cases':s['cases'],'accepted':s['accepted'],'rejected':s['rejected'],'max_direction_error_accepted':s['max_direction_error_accepted'],'mutation_cases':s['mutation_cases']}))
