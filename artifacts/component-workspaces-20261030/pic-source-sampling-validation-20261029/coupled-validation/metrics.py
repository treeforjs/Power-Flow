"""Physical-weight diagnostics and conservative decision gates; no solver calls."""
import math
import random
import statistics
import struct
from pathlib import Path

E=1.602176634e-19
C=299792458.
REGIONS=('all','gap','lower_back','upper_back','left','right')
LEDGER=('time_s','particle_charge_C','deposited_charge_C','lower_field_charge_C','upper_field_charge_C',
        'gauss_defect_C','particle_charge_defect_C','field_energy_J','particle_kinetic_energy_J',
        'lower_absorbed_charge_C','upper_absorbed_charge_C','lower_absorbed_energy_J','upper_absorbed_energy_J',
        'lower_injected_charge_C','upper_injected_charge_C','injected_energy_J','authored_charge_C',
        'authored_kinetic_energy_J','maximum_nodal_continuity_defect_C','escaped_charge_C',
        'escaped_kinetic_energy_J','short_absorbed_charge_C','short_absorbed_energy_J',
        'axial_escaped_charge_C','poisson_residual_C_per_m')

def finite(x):
    if not math.isfinite(x):raise ValueError('Nonfinite diagnostic')
    return x

def weighted_quantile(samples,q):
    if not 0<=q<=1:raise ValueError('Invalid quantile')
    total=math.fsum(w for x,w in samples)
    if not total:return None
    run=0
    for x,w in sorted(samples):
        if w<0:raise ValueError('Negative physical weight')
        run+=w
        if run>=q*total:return x
    return max(x for x,w in samples)

def weighted_cdf_distance(a,b):
    na=math.fsum(w for x,w in a);nb=math.fsum(w for x,w in b)
    if na<=0 or nb<=0:return None
    changes={}
    for x,w in a:changes[x]=changes.get(x,0)+w/na
    for x,w in b:changes[x]=changes.get(x,0)-w/nb
    value=0.;worst=0.
    for x in sorted(changes):value+=changes[x];worst=max(worst,abs(value))
    return min(1.,worst)

def region(x,y):
    if x<=.002:return 'left'
    if x>=.007:return 'right'
    if y<=.0018:return 'lower_back'
    if y>=.0072:return 'upper_back'
    if .002<=y<=.007:return 'gap'
    raise ValueError('Particle inside excluded foil')

def kinetic(p):
    # Native particle tuple: id,xyz,proper-velocity xyz,q,m,physical weight.
    u2=sum(x*x for x in p[4:7])
    return p[8]*u2/(math.sqrt(1+u2/(C*C))+1)/E

def read_frame(path):
    data=Path(path).read_bytes();offset=0
    def take(fmt):
        nonlocal offset
        size=struct.calcsize(fmt)
        if offset+size>len(data):raise ValueError('Truncated frame')
        value=struct.unpack_from(fmt,data,offset);offset+=size
        return value
    magic,time,n,fresh=take('<QdQQ')
    if magic!=0x43564652414d4531 or n>2000000:raise ValueError('Frame header invalid')
    ledger=dict(zip(LEDGER,take('<25d')))
    particles=[take('<Q9d') for i in range(n)]
    charges=take('<5d');energies=take('<5d');supplied,pending=take('<2d')
    cells,=take('<Q')
    if cells!=360*360:raise ValueError('Field mesh changed')
    fields=[take('<4d') for i in range(cells)]
    if offset!=len(data):raise ValueError('Trailing frame bytes')
    for x in [time,supplied,pending,*ledger.values(),*charges,*energies]:finite(x)
    ids=set()
    for p in particles:
        if p[0] in ids or p[7]!=-E or p[8]<=0 or p[9]<=0:raise ValueError('Invalid particle/id/weight')
        ids.add(p[0]);[finite(x) for x in p[1:]]
    for row in fields:
        [finite(x) for x in row]
        if row[3]<0:raise ValueError('Negative vacuum volume')
    return dict(time_s=time,fresh_start=fresh,ledger=ledger,particles=particles,
                fresh_sink_charge_C=charges,fresh_sink_energy_J=energies,
                source_electrons=supplied,pending_electrons=pending,fields=fields)

def populations(frame):
    rows={(part,pop,r):dict(samples=[],bins=[0.]*324,N=0.,sum_w2=0.,kinetic_eV=0.)
          for part in ['full','slab'] for pop in ['all','fresh'] for r in REGIONS}
    for p in frame['particles']:
        ident,x,y,z,*rest=p;w=p[9];energy=kinetic(p);r=region(x,y)
        ix=min(17,max(0,int(x/.0005)));iy=min(17,max(0,int(y/.0005)))
        for part in (['full','slab'] if abs(z)<.0005 else ['full']):
            for pop in (['all','fresh'] if ident>=frame['fresh_start'] else ['all']):
                for area in ['all',r]:
                    v=rows[part,pop,area];v['N']+=w;v['sum_w2']+=w*w;v['kinetic_eV']+=w*energy
                    v['samples'].append((energy,w));v['bins'][iy*18+ix]+=w
    for v in rows.values():
        v['N_eff']=v['N']**2/v['sum_w2'] if v['sum_w2'] else 0.
        v['mean_energy_eV']=v['kinetic_eV']/v['N'] if v['N'] else None
        for q in [.1,.5,.9]:v[f'q{int(q*100)}_energy_eV']=weighted_quantile(v['samples'],q)
    for (part,pop,r),v in rows.items():
        depth=.025 if part=='full' else .001
        areas=dict(all=.009*.009-2*.005*.0002,gap=.005*.005,lower_back=.005*.0018,upper_back=.005*.0018,left=.002*.009,right=.002*.009)
        v['volume_m3']=areas[r]*depth;v['density_cm3']=v['N']/v['volume_m3']/1e6
    return rows

def frame_observables(frame,initial):
    """Terminal/native frame -> scalar observables plus weighted distributions."""
    pop=populations(frame);values={}
    for (part,group,r),v in pop.items():
        tag=f'{part}:{group}:{r}';values['N:'+tag]=v['N'];values['density_cm3:'+tag]=v['density_cm3']
        if r=='all':
            for name in ['kinetic_eV','mean_energy_eV','q10_energy_eV','q50_energy_eV','q90_energy_eV']:values[name+':'+tag]=v[name]
    for name in ['lower_absorbed_charge_C','upper_absorbed_charge_C','short_absorbed_charge_C','escaped_charge_C','lower_injected_charge_C','upper_injected_charge_C']:
        values['delta_electrons:'+name]=-(frame['ledger'][name]-initial['ledger'][name])/E
    for name in ['lower_absorbed_energy_J','upper_absorbed_energy_J','short_absorbed_energy_J','escaped_kinetic_energy_J','injected_energy_J']:
        values['delta_energy_eV:'+name]=(frame['ledger'][name]-initial['ledger'][name])/E
    for name,index in [('lower',0),('upper',1),('short',3)]:
        values['fresh_sink_electrons:'+name]=-frame['fresh_sink_charge_C'][index]/E
        values['fresh_sink_energy_eV:'+name]=frame['fresh_sink_energy_J'][index]/E
    # hit4 is already an escape. Add side escape2 and axial escape4 once.
    values['fresh_sink_electrons:escape']=-(frame['fresh_sink_charge_C'][2]+frame['fresh_sink_charge_C'][4])/E
    values['fresh_sink_energy_eV:escape']=(frame['fresh_sink_energy_J'][2]+frame['fresh_sink_energy_J'][4])/E
    values['source_electrons']=frame['source_electrons'];values['pending_electrons']=frame['pending_electrons']
    return dict(scalars=values,populations=pop)

def paired_distances(a,b):
    result={};pa=populations(a);pb=populations(b)
    for part in ['full','slab']:
        for group in ['all','fresh']:
            x,y=pa[part,group,'all'],pb[part,group,'all'];tag=part+':'+group
            result['density_L1:'+tag]=math.fsum(abs(i-j) for i,j in zip(x['bins'],y['bins']))/y['N'] if y['N'] else None
            result['energy_CDF:'+tag]=weighted_cdf_distance(x['samples'],y['samples']) if min(x['N_eff'],y['N_eff'])>=50 else None
    fields=field_distances(a['fields'],b['fields'])
    result['field_potential']=fields['potential_normalized'];result['field_electric']=fields['electric_normalized']
    return result

def field_distances(a,b):
    if len(a)!=len(b):raise ValueError('Mismatched field meshes')
    sums=[0.]*5;volume=0.
    for x,y in zip(a,b):
        if x[3]!=y[3]:raise ValueError('Mismatched physical volumes')
        v=x[3];volume+=v
        for k,value in enumerate([(x[0]-y[0])**2,(x[1]-y[1])**2+(x[2]-y[2])**2,y[0]**2,y[1]**2+y[2]**2]):sums[k]+=v*value
    if volume<=0:raise ValueError('Empty field domain')
    phi,field,phi_ref,field_ref=[math.sqrt(v/volume) for v in sums[:4]]
    return dict(potential_RMS_V=phi,electric_RMS_V_m=field,
                potential_normalized=phi/max(phi_ref,.001),electric_normalized=field/max(field_ref,.2))

def trapezoid(times,values):
    if len(times)!=len(values) or len(times)<2 or any(b<=a for a,b in zip(times,times[1:])):raise ValueError('Bad averaging grid')
    if any(v is None for v in values):return None
    return math.fsum((b-a)*(x+y)/2 for a,b,x,y in zip(times,times[1:],values,values[1:]))/(times[-1]-times[0])

def scalar_gate(a,r,absolute_margin=None):
    if len(a)!=len(r) or len(a)!=6 or any(x is None for x in a+r):return dict(status='inconclusive',reason='six complete paired seeds required')
    if any(not math.isfinite(x) for x in a+r):raise ValueError('Nonfinite confidence input')
    # t(0.975,df5), also the one-sided97.5% bound requested by the proposal.
    t=2.570581835636305
    if absolute_margin is None:
        hi=[x-1.05*y for x,y in zip(a,r)];lo=[.95*y-x for x,y in zip(a,r)]
    else:
        hi=[x-y-absolute_margin for x,y in zip(a,r)];lo=[y-x-absolute_margin for x,y in zip(a,r)]
    def bounds(v):
        m=statistics.mean(v);h=t*statistics.stdev(v)/math.sqrt(6);return [m-h,m+h]
    upper,lower=bounds(hi),bounds(lo)
    status='pass' if upper[1]<0 and lower[1]<0 else 'fail' if upper[0]>0 or lower[0]>0 else 'inconclusive'
    return dict(status=status,upper_contrast_ci95=upper,lower_contrast_ci95=lower,
                paired_mean_difference=statistics.mean(x-y for x,y in zip(a,r)))

def distance_gate(values,limit,seed=1973):
    if len(values)!=6 or any(v is None for v in values):return dict(status='inconclusive',reason='six resolved paired distances required')
    if any(not math.isfinite(v) or v<0 for v in values):raise ValueError('Bad distance')
    # Deterministic resampling of SEED metrics only, never macros or times.
    rng=random.Random(seed);means=sorted(statistics.mean(rng.choices(values,k=6)) for i in range(5000))
    lower,upper=means[125],means[4874]
    mean=statistics.mean(values);sd=statistics.stdev(values)
    skew=statistics.mean((x-mean)**3 for x in values)/(sd**3) if sd else 0.
    unstable=abs(skew)>1.5
    status='inconclusive' if unstable else 'pass' if upper<=limit else 'fail' if lower>limit else 'inconclusive'
    return dict(status=status,bootstrap_ci95_mean=[lower,upper],approximate_coverage=True,
                seed_count=6,skewness=skew,unstable=unstable)

def evaluate_scalar_table(records):
    """Intake contract for paired terminal/averaged metrics, fixed before results.

    Each record has count, seed, stage, endpoint, metric, value and an optional
    absolute_margin. Missing seed/metric values remain inconclusive.
    """
    groups={}
    for row in records:
        if row['count'] not in [64,32,16]:raise ValueError('Unsupported count')
        key=(row['stage'],row['endpoint'],row['metric']);pairs=groups.setdefault(key,{})
        slot=(row['count'],row['seed'])
        if slot in pairs:raise ValueError('Duplicate paired metric')
        pairs[slot]=row
    results=[]
    for key,pairs in sorted(groups.items()):
        seeds=range(6,12) if key[0]=='native_adaptive_confirmation' else range(6)
        for candidate,reference in [(32,64),(16,64),(16,32)]:
            a=[pairs.get((candidate,s),{}).get('value') for s in seeds]
            r=[pairs.get((reference,s),{}).get('value') for s in seeds]
            margins={row.get('absolute_margin') for row in pairs.values()}
            if len(margins)>1:raise ValueError('Margin changed inside a paired family')
            result=scalar_gate(a,r,next(iter(margins)))
            results.append(dict(stage=key[0],endpoint=key[1],metric=key[2],candidate=candidate,reference=reference,**result))
    return dict(status='inconclusive' if not results else 'fail' if any(r['status']=='fail' for r in results) else 'inconclusive' if any(r['status']!='pass' for r in results) else 'pass',
                scalar_results=results,full_campaign_acceptance=False,
                remaining_requirements=['complete primary metric inventory','distance/field confidence gates','fresh-region20% gates','N_eff50 for distributions','native numerical validity','both stages complete'])

def primary_names():
    names={f'N:{part}:{pop}:{r}' for part in ['full','slab'] for pop in ['all','fresh'] for r in REGIONS}
    names|={f'{k}:{part}:{pop}:all' for part in ['full','slab'] for pop in ['all','fresh'] for k in ['kinetic_eV','mean_energy_eV','q10_energy_eV','q50_energy_eV','q90_energy_eV']}
    names|={'delta_electrons:'+k for k in ['lower_absorbed_charge_C','upper_absorbed_charge_C','short_absorbed_charge_C','escaped_charge_C','lower_injected_charge_C','upper_injected_charge_C']}
    names|={'delta_energy_eV:'+k for k in ['lower_absorbed_energy_J','upper_absorbed_energy_J','short_absorbed_energy_J','escaped_kinetic_energy_J','injected_energy_J']}
    names|={f'fresh_sink_{kind}:{r}' for kind in ['electrons','energy_eV'] for r in ['lower','upper','short','escape']}
    names|={'source_electrons','pending_electrons'}
    return names

def evaluate_bundle(bundle):
    """Fail-closed campaign intake; produces a review result, never deployment."""
    stages=['common_trace_count_isolation','native_adaptive_confirmation'];endpoints=['terminal_100ps','trapezoidal_mean_50_to_100ps']
    problems=[];rows=bundle.get('scalar_records',[]);keys={(r['stage'],r['endpoint'],r['metric']) for r in rows}
    for stage in stages:
        for endpoint in endpoints:
            for metric in primary_names():
                if (stage,endpoint,metric) not in keys:problems.append('missing_primary:'+stage+':'+endpoint+':'+metric)
    valid_cases=set()
    for run in bundle.get('runs',[]):
        key=(run['stage'],run['count'],run['seed'])
        if key in valid_cases:raise ValueError('Duplicate campaign run')
        if run.get('complete') and run.get('native_guards_passed') and run.get('accepted_duration_ps',0)>=100:valid_cases.add(key)
    for stage in stages:
        for seed in (range(6) if stage==stages[0] else range(6,12)):
            for count in [64,32,16]:
                if (stage,count,seed) not in valid_cases:problems.append(f'missing_or_invalid_run:{stage}:{count}:{seed}')
    if problems:return dict(status='inconclusive',candidate_accepted_for_review=False,missing=problems)
    scalar=evaluate_scalar_table(rows);distances=[]
    distance_names={f'{k}:{part}:{pop}' for k in ['density_L1','energy_CDF'] for part in ['full','slab'] for pop in ['all','fresh']}|{'field_potential','field_electric'}
    supplied=bundle.get('distance_records',[])
    for stage in stages:
        seeds=list(range(6) if stage==stages[0] else range(6,12))
        for endpoint in endpoints:
            for name in sorted(distance_names):
                for a,b in [(32,64),(16,64),(16,32)]:
                    matching=[r for r in supplied if (r['stage'],r['endpoint'],r['metric'],r['candidate'],r['reference'])==(stage,endpoint,name,a,b)]
                    by_seed={r['seed']:r['value'] for r in matching}
                    if len(by_seed)!=len(matching):raise ValueError('Duplicate distance seed')
                    limit=.1 if name.startswith('density_L1') else .05
                    distances.append(dict(stage=stage,endpoint=endpoint,metric=name,candidate=a,reference=b,**distance_gate([by_seed.get(s) for s in seeds],limit)))
    infos=bundle.get('information_records',[])
    for stage,count,seed in valid_cases:
        for part in ['full','slab']:
            for r in REGIONS:
                found=[v for v in infos if (v['stage'],v['count'],v['seed'],v['partition'],v['region'])==(stage,count,seed,part,r)]
                if len(found)!=1 or ((r=='all' or not found[0].get('near_zero',False)) and found[0].get('fresh_fraction',0)<.2):problems.append(f'fresh_region_unresolved:{stage}:{count}:{seed}:{part}:{r}')
    if bundle.get('data_status')!='owner_supplied':problems.append('manufactured_or_unattributed_data')
    statuses=[scalar['status']]+[r['status'] for r in distances]
    status='fail' if 'fail' in statuses else 'inconclusive' if problems or any(s!='pass' for s in statuses) else 'pass_for_scientific_review'
    return dict(status=status,candidate_accepted_for_review=status=='pass_for_scientific_review',scalar=scalar,distance_results=distances,
                unresolved=problems,deployment_authorized=False,whole_history_convergence=False,independent_execution_verified=False)

def self_test(where=None):
    assert weighted_quantile([(1,1),(100,9)],.5)==100
    assert weighted_cdf_distance([(1,.5),(2,.5)],[(1,1),(2,1)])<1e-14
    assert weighted_cdf_distance([(1,1)],[(2,1)])==1
    assert weighted_cdf_distance([],[(1,1)]) is None
    assert trapezoid([0,1,2],[0,2,4])==2
    assert scalar_gate([100]*6,[100]*6)['status']=='pass'
    assert scalar_gate([120]*6,[100]*6)['status']=='fail'
    assert scalar_gate([100]*3,[100]*3)['status']=='inconclusive'
    assert scalar_gate([None]*6,[100]*6)['status']=='inconclusive'
    assert distance_gate([.01]*6,.05)['status']=='pass'
    assert distance_gate([.1]*6,.05)['status']=='fail'
    assert evaluate_scalar_table([])['status']=='inconclusive'
    assert evaluate_bundle({})['status']=='inconclusive'
    assert field_distances([(1,2,3,4)],[(1,2,3,4)])['potential_normalized']==0
    if where is not None:
        sample=(1,.003,.003,0.,0.,0.,0.,-E,9.1093837139e-31,1.)
        fresh=(2,.001,.003,.001,0.,0.,0.,-E,9.1093837139e-31,.5)
        data=struct.pack('<QdQQ',0x43564652414d4531,1e-9,2,2)+struct.pack('<25d',*([0.]*25))
        data+=struct.pack('<Q9d',*sample)+struct.pack('<Q9d',*fresh)+struct.pack('<12d',*([0.]*12))+struct.pack('<Q',360*360)
        data+=struct.pack('<4d',0.,0.,0.,1e-12)*(360*360)
        path=Path(where)/'manufactured-frame.bin';path.write_bytes(data);frame=read_frame(path);pop=populations(frame)
        assert pop['full','all','all']['N']==1.5 and pop['slab','all','all']['N']==1
        assert pop['full','fresh','left']['N']==.5 and pop['slab','fresh','all']['N']==0
        assert frame_observables(frame,frame)['scalars']['fresh_sink_electrons:escape']==0
        assert paired_distances(frame,frame)['density_L1:full:all']==0
        assert paired_distances(frame,frame)['energy_CDF:full:all'] is None
        bad=Path(where)/'truncated-frame.bin';bad.write_bytes(data[:-1])
        try:read_frame(bad)
        except ValueError:pass
        else:raise AssertionError('Truncated diagnostic frame accepted')
    return dict(passed=True,weighted_metrics=True,paired_gate_screen_never_passes=True,empty_inputs_inconclusive=True,
                frame_codec_fixture=where is not None,physical_data_used=False)
