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
    m,M,mh=masses;mo=mh;Q=D.from_float(meta['reference_cost_eV'])
    def load(name,path):
        spec=importlib.util.spec_from_file_location(name,path);mod=importlib.util.module_from_spec(spec);sys.modules[name]=mod;spec.loader.exec_module(mod);return mod
    h2=load('h2_primary',Path('/home/treefor/projects/runtimes/h2-dissociation-primary-v1.1-20261101/src/h2_primary.py'))
    cdfmod=load('primary_cdf',ROOT/'src/primary_cdf.py')
    refs=Path('/home/treefor/projects/runtimes/water-fragment-source-inventory-20261101/references')
    cdf=cdfmod.load(refs/'sartori2024-table-I.txt');quant=[D.from_float(float(q)) for q,e in cdf];frag=[D.from_float(float(e)) for q,e in cdf]
    lossmin=Q+2*frag[0];T0=lossmin*(2*(m+M)+lossmin)/(2*M)
    rows=[(e,s) for e,s in h2._parse_source((refs/'MCCC-el-H2-b3Su_DE.X1Sg_vi=0.txt').read_bytes()).rows if e>=D('4.5')];es=[e for e,_ in rows]
    def interpolate(xs,ys,x):
        i=max(0,min(len(xs)-2,bisect_right(xs,x)-1))
        return ys[i]+(ys[i+1]-ys[i])*(x-xs[i])/(xs[i+1]-xs[i])
    a0=D.from_float(meta['bohr_radius_m'])
    maxima={k:D(0) for k in ('chemical_energy_relative','momentum_relative','H_rest_partition_abs_eV','H2_rest_partition_abs_eV','electron_CM_abs_eV','rate_relative','CDF_interpolation_abs_eV','conditional_upper_quantile_abs')}
    mean_e=mean_f=D(0);count=rate_count=threshold_count=0
    uniforms=[];electron_cosines=[]
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
        if r['type']=='controls':assert r=={'type':'controls','refusals':20,'events':4096};continue
        assert r['type']=='event';count+=1
        E=dec(r['E']);Kcm=2*M*E/(((m+M)**2+2*M*E).sqrt()+m+M);cap=(Kcm-Q)/2
        upper=min(quant[-1],interpolate(frag,quant,cap))
        assert cap>frag[0] and upper>quant[0]
        uerr=abs(dec(r['quantile_hi'])-upper);assert uerr<D('2e-14');maxima['conditional_upper_quantile_abs']=max(maxima['conditional_upper_quantile_abs'],uerr)
        qdraw=dec(r['quantile']);assert quant[0]<qdraw<upper
        fragment=interpolate(quant,frag,qdraw);assert frag[0]<=fragment<cap
        loss=dec(r['loss']);K=loss-Q;KH=KO=K/2
        perr=abs(K/2-fragment);assert perr<D('1e-14');maxima['CDF_interpolation_abs_eV']=max(maxima['CDF_interpolation_abs_eV'],perr)
        assert abs(dec(r['KER'])-K)<D('2e-15')
        uniforms.append((qdraw-quant[0])/(upper-quant[0]))
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
        for label,actual,expected in [('H',kinetic(mh,ph),KH),('H2',kinetic(mo,po),KO)]:
            err=abs(actual-expected);assert err<D('2e-12')
            maxima[label+'_rest_partition_abs_eV']=max(maxima[label+'_rest_partition_abs_eV'],err)
        # Pair COM electron energy follows exact two-body recoil, not T-loss in the lab.
        Pi=add(pin[0],pin[1]);Ei=m+M+sum(kin);Mi=(Ei*Ei-dot(Pi,Pi)).sqrt();wi=scale(Pi,1/Mi)
        pe=boost_to_frame(m,pout[0],wi)
        T=dec(r['E']);kcm=2*M*T/(((m+M)**2+2*M*T).sqrt()+m+M);koutcm=kcm-loss
        expected_e=koutcm*(2*(M+loss)+koutcm)/(2*(m+M+loss+koutcm))
        err=abs(kinetic(m,pe)-expected_e);assert err<D('2e-12')
        maxima['electron_CM_abs_eV']=max(maxima['electron_CM_abs_eV'],err)
        cosine=pe[2]/dot(pe,pe).sqrt();mean_e+=cosine;mean_f+=ph[2]/dot(ph,ph).sqrt();electron_cosines.append(cosine)
    assert count==4096 and rate_count==2*len(rows)-1 and threshold_count==2
    assert abs(mean_e/count)<D('0.06') and abs(mean_f/count)<D('0.06')
    us=sorted(uniforms);n=D(count);KS=max(max(abs(u-D(i)/n),abs(u-D(i+1)/n)) for i,u in enumerate(us));assert KS<D('.03')
    mu=sum(uniforms)/n;mc=sum(electron_cosines)/n
    corr=sum((u-mu)*(c-mc) for u,c in zip(uniforms,electron_cosines))/(sum((u-mu)**2 for u in uniforms)*sum((c-mc)**2 for c in electron_cosines)).sqrt()
    assert abs(corr)<D('.06')
    return {'schema':'h2-cdf-diagnostic-independent-controls-v1' ,'events':count,'rates':rate_count,
            'negative_controls':20,'conditional_uniform_KS':str(KS),'partition_angle_correlation':str(corr),'solver_threshold_controls':2,'maxima':{k:str(v) for k,v in maxima.items()},
            'primary_CDF_knots':len(cdf),'omitted_probability':'0.002',
            'electron_CM_mean_cosine':str(mean_e/count),'fragment_mean_cosine':str(mean_f/count),
            'passed':True,'production_qualified':False,'physical_model_validated':False,'network_complete':False,'native_adopted':False}
if __name__=='__main__':
    with localcontext(CTX):print(json.dumps(main(sys.argv[1]),sort_keys=True,indent=2))
