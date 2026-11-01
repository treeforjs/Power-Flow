"""Independent Decimal invariant/partition checks; no shared solver calls."""
from bisect import bisect_right
from decimal import Decimal as D, Context, ROUND_HALF_EVEN, localcontext
import importlib.util
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
CTX=Context(prec=70,rounding=ROUND_HALF_EVEN,Emin=-999999,Emax=999999)
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def add(a,b):return [x+y for x,y in zip(a,b)]
def scale(a,s):return [x*s for x in a]
def kinetic(m,p):return dot(p,p)/((m*m+dot(p,p)).sqrt()+m)
def boost_to_frame(m,p,w):
    gamma=(1+dot(w,w)).sqrt()
    return add(p,scale(w,dot(w,p)/(gamma+1)-(m+kinetic(m,p))))
def dec(x):return D.from_float(x) if isinstance(x,float) else D(x)
def main(path):
    meta=json.loads((ROOT/'build/model.json').read_text());sp=meta['species']
    masses=[D.from_float(s['mass_eV']) for s in sp]
    m,M,mh,mo=masses;loss=D.from_float(meta['solver_loss_binary64_eV']);Q=D.from_float(meta['reference_cost_eV'])
    # Independent exact two-body kinetic-energy division in the fragment rest frame.
    K=loss-Q;S=mh+mo+K
    KH=K*(2*mo+K)/(2*S);KO=K*(2*mh+K)/(2*S)
    assert abs(KH+KO-K)<D('1e-60')
    T0=loss*(2*(m+M)+loss)/(2*M)
    loader=Path('/home/treefor/projects/runtimes/oh-primary-supplement-v1.2-20261101/src/oh_primary.py')
    spec=importlib.util.spec_from_file_location('oh_primary',loader);mod=importlib.util.module_from_spec(spec);sys.modules['oh_primary']=mod;spec.loader.exec_module(mod)
    rows=mod.load_channel(Path('/home/treefor/projects/runtimes/water-fragment-source-inventory-20261101/references/PSSTab264csupp_e-OH-data.xlsx'),'excitation_C').rows
    rows=[(e,s) for e,s in rows if e>=D('7.3605')];es=[e for e,_ in rows]
    a0=D.from_float(meta['bohr_radius_m'])
    maxima={k:D(0) for k in ('chemical_energy_relative','momentum_relative','H_rest_partition_abs_eV','O_rest_partition_abs_eV','electron_CM_abs_eV','rate_relative')}
    mean_e=mean_f=D(0);count=rate_count=threshold_count=0
    for line in Path(path).read_text().splitlines():
        r=json.loads(line)
        if r['type']=='rate':
            # Source decimal abscissas are converted to binary64 by the pinned importer.
            E=dec(r['E']);points=[D.from_float(float(e)) for e in es];i=max(0,min(len(rows)-2,bisect_right(points,E)-1))
            e1,e2=points[i:i+2];s1,s2=[s*a0*a0 for _,s in rows[i:i+2]]
            expected=s1+(s2-s1)*(E-e1)/(e2-e1)
            err=abs(dec(r['sigma'])-expected)/expected
            assert err<D('1e-14');maxima['rate_relative']=max(maxima['rate_relative'],err);rate_count+=1;continue
        if r['type']=='threshold':
            assert bool(r['accepted'])==(dec(r['E'])>T0);threshold_count+=1;continue
        if r['type']=='controls':assert r=={'type':'controls','refusals':16,'events':2048};continue
        assert r['type']=='event';count+=1
        pin=[[dec(v) for v in p] for p in r['incoming']];pout=[[dec(v) for v in p] for p in r['outgoing']]
        kin=[kinetic(m,pin[0]),kinetic(M,pin[1])];kout=[kinetic(mm,p) for mm,p in zip((m,mh,mo),pout)]
        energy_scale=sum(kin)+sum(kout)+abs(Q);defect=sum(kin)-sum(kout)-Q
        relative=abs(defect)/energy_scale;assert relative<=D(2)**-40
        maxima['chemical_energy_relative']=max(maxima['chemical_energy_relative'],relative)
        pscale=sum(abs(x) for p in pin+pout for x in p)
        for j in range(3):
            err=abs(sum(p[j] for p in pin)-sum(p[j] for p in pout))/pscale
            assert err<D(2)**-40;maxima['momentum_relative']=max(maxima['momentum_relative'],err)
        assert r['internal']==0
        for reported,expected in zip(r['kinetic_in']+r['kinetic_out'],kin+kout):assert abs(dec(reported)-expected)<=D('5e-14')*max(D(1),expected)
        # Boost physical product momenta into their own exact fragment-pair COM.
        Pf=add(pout[1],pout[2]);Ef=mh+mo+kout[1]+kout[2]
        Mf=(Ef*Ef-dot(Pf,Pf)).sqrt();wf=scale(Pf,1/Mf)
        ph=boost_to_frame(mh,pout[1],wf);po=boost_to_frame(mo,pout[2],wf)
        for label,actual,expected in [('H',kinetic(mh,ph),KH),('O',kinetic(mo,po),KO)]:
            err=abs(actual-expected);assert err<D('2e-12')
            maxima[label+'_rest_partition_abs_eV']=max(maxima[label+'_rest_partition_abs_eV'],err)
        # Pair COM electron energy follows exact two-body recoil, not T-loss in the lab.
        Pi=add(pin[0],pin[1]);Ei=m+M+sum(kin);Mi=(Ei*Ei-dot(Pi,Pi)).sqrt();wi=scale(Pi,1/Mi)
        pe=boost_to_frame(m,pout[0],wi)
        T=dec(r['E']);kcm=2*M*T/(((m+M)**2+2*M*T).sqrt()+m+M);koutcm=kcm-loss
        expected_e=koutcm*(2*(M+loss)+koutcm)/(2*(m+M+loss+koutcm))
        err=abs(kinetic(m,pe)-expected_e);assert err<D('2e-12')
        maxima['electron_CM_abs_eV']=max(maxima['electron_CM_abs_eV'],err)
        mean_e+=pe[2]/dot(pe,pe).sqrt();mean_f+=ph[2]/dot(ph,ph).sqrt()
    assert count==2048 and rate_count==2*len(rows)-1 and threshold_count==2
    assert abs(mean_e/count)<D('0.06') and abs(mean_f/count)<D('0.06')
    return {'schema':'oh-diagnostic-independent-controls-v1','events':count,'rates':rate_count,
            'negative_controls':16,'solver_threshold_controls':2,'maxima':{k:str(v) for k,v in maxima.items()},
            'expected_H_fragment_rest_KE_eV':str(KH),'expected_O_fragment_rest_KE_eV':str(KO),
            'electron_CM_mean_cosine':str(mean_e/count),'fragment_mean_cosine':str(mean_f/count),
            'passed':True,'production_qualified':False,'physical_model_validated':False,'network_complete':False,'native_adopted':False}
if __name__=='__main__':
    with localcontext(CTX):print(json.dumps(main(sys.argv[1]),sort_keys=True,indent=2))
