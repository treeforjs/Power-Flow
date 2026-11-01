import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from decimal import Decimal as D, localcontext
ROOT=Path(__file__).resolve().parents[1]
REF=Path('/home/treefor/projects/runtimes/water-fragment-source-inventory-20261101/references')
DATA=REF/'alecs-06695e72'
CHEM=Path('/home/treefor/projects/runtimes/pic-water-chemistry-v3.1.1-20261031')
ARR=Path('/home/treefor/projects/runtimes/pic-arrival-owner-v1-20261030/include/pf_arrival')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify(pins):
    for section in ('source_files','dependencies'):
        for n,h in pins[section].items():assert sha(ROOT/n)==h,'input changed: '+n
    assert sha(Path(sys.executable).resolve())==pins['frozen_python_sha256']
def read_xs(path,columns):
    raw=path.read_bytes();assert len(raw)<16384
    lines=raw.decode('ascii').splitlines();assert lines[0]=='# OH cross section' and lines[1]=='256' and len(lines)==258
    rows=[tuple(map(D,line.split())) for line in lines[2:]]
    assert all(len(x)==columns and all(v.is_finite() and v>=0 for v in x) for x in rows)
    assert all(rows[i][0]>rows[i-1][0] for i in range(1,len(rows)))
    return rows
def half_print_quantum(v):return D(5).scaleb(v.as_tuple().exponent-1)
def model_bound(m,orbitals,E):
    # Conservative interval propagation of printed operand precision, not physical uncertainty.
    elo,ehi=E-half_print_quantum(E),E+half_print_quantum(E)
    lo=hi=D(0)
    for o in orbitals:
        bl,bh=o.B-half_print_quantum(o.B),o.B+half_print_quantum(o.B)
        ul,uh=o.U-half_print_quantum(o.U),o.U+half_print_quantum(o.U)
        if ehi<=bl:continue
        tl=max(D(1),elo/bh);th=max(D(1),ehi/bl)
        flo=tl.ln()/2*(1-1/tl**2)+1-1/tl-th.ln()/(tl+1)
        fhi=th.ln()/2*(1-1/th**2)+1-1/th-tl.ln()/(th+1)
        slo=4*m.PI*(m.RYD/bh)**2
        shi=4*m.PI*(m.RYD/bl)**2
        lo+=slo*max(D(0),flo)/(th+uh/bl+1)
        hi+=shi*max(D(0),fhi)/(tl+ul/bh+1)
    return lo,hi
def main():
    job=os.environ.get('SIMQ_JOB_ID');assert job and job.isdigit() and len(os.sched_getaffinity(0))==1
    assert not os.environ.get('PF_BUILD_ISA_FLAGS') and not os.environ.get('CUDA_VISIBLE_DEVICES')
    pins=json.loads((ROOT/'INPUTS.json').read_text());verify(pins)
    out=ROOT/'evidence'/('job-'+job);out.mkdir(exist_ok=False);start=time.monotonic()
    def run(args,label):
        p=subprocess.run(args,cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=90)
        (out/(label+'.stdout')).write_bytes(p.stdout);(out/(label+'.stderr')).write_bytes(p.stderr)
        assert p.returncode==0,p.stderr.decode();return p.stdout
    run([sys.executable,'-I',str(ROOT/'tests/test_model.py')],'tests')
    spec=importlib.util.spec_from_file_location('oh_beb',ROOT/'src/oh_beb.py');m=importlib.util.module_from_spec(spec);sys.modules['oh_beb']=m;spec.loader.exec_module(m)
    model=m.load(DATA/'NIST_orbitals/mp2/OH.norb')
    recommended=read_xs(DATA/'ion_xs/recommended/OH.xs',2);mp2=read_xs(DATA/'ion_xs/mp2/OH.xs',3)
    assert recommended==[(row[0],row[1]) for row in mp2],'recommended OH must equal raw MP2 BEB column, not damped column'
    tested=[row for row in recommended if row[0]<=5000]
    energies=[float(e) for e,_ in tested]+[70.,5000.]
    header='// Local source operands, not redistributed.\n'
    header+='const double orbitals[][3]={'+','.join('{'+','.join(float(x).hex() for x in (o.B,o.U,o.N))+'}' for o in model.orbitals)+'};\n'
    header+='const double energies[]={'+','.join(e.hex() for e in energies)+'};\n'
    (ROOT/'build/operands.hpp').write_text(header)
    run(['/usr/bin/c++','-std=c++20','-O2','-Wall','-Wextra','-Wpedantic','-Werror','-ffp-contract=off','-fno-fast-math',
         '-I'+str(CHEM/'include/pf_chem_v3'),'-I'+str(ARR),'-I'+str(ROOT/'build'),str(ROOT/'tests/kernel.cpp'),'-o',str(ROOT/'build/kernel')],'compile')
    # Keep source-derived numerical arrays local; public evidence contains only hashes and aggregates.
    raw=subprocess.check_output([str(ROOT/'build/kernel')]);(ROOT/'build/kernel-values.txt').write_bytes(raw)
    assert raw==subprocess.check_output([str(ROOT/'build/kernel')]),'kernel replay'
    kernel=[tuple(map(D,line.split())) for line in raw.decode().splitlines()]
    with localcontext(m.CTX):
        max_kernel=max_print=D(0);positive=0
        for e,s in kernel:
            calculated=model.diagnostic_total(m.Incident(e),mode='source-BEB-total-only').sigma_m2
            if calculated==0:assert s==0
            else:
                error=abs(s-calculated)/calculated;assert error<D('2e-14');max_kernel=max(max_kernel,error)
        for e,s in tested:
            lo,hi=model_bound(m,model.orbitals,e);half=half_print_quantum(s)
            assert lo<=s+half and hi>=s-half,(str(e),str(s),str(lo),str(hi))
            calculated=model.diagnostic_total(m.Incident(e),mode='source-BEB-total-only').sigma_m2/(m.A0_M*m.A0_M)
            if s>0:positive+=1;max_print=max(max_print,abs(calculated-s)/s)
        registry=json.loads((CHEM/'data/tables/chemistry-v3.1-manifest.json').read_text());species={s['name']:s for s in registry['species']}
        q=species['OH+']['ref_eV']-species['OH']['ref_eV']
        me=D.from_float(species['e']['mass_eV']);M=D.from_float(species['OH']['mass_eV']);Q=D.from_float(q)
        threshold=Q*(2*(me+M)+Q)/(2*M)
        source_lo=min(o.B for o in model.orbitals)
        above=next(e for e,s in tested if e>=threshold and s>0)
        fragment_costs={name:D.from_float(species[a]['ref_eV'])+D.from_float(species[b]['ref_eV'])-D.from_float(species['OH']['ref_eV']) for name,a,b in [('Oplus_H','O+','H'),('O_Hplus','O','H+')]}
        report={'schema':'oh-primary-beb-qualification-v1.1','primary_repository':'alecs-db/ALeCS-db','primary_commit':'06695e72aea538e82b896c14e5dad01900d38a0d',
                'source_identity':pins['source_identity'],'orbitals':len(model.orbitals),'occupied_electrons':sum(o.N for o in model.orbitals),
                'alpha_beta_rows_preserved':True,'source_BEB_Rydberg_eV':str(m.RYD),'source_Bohr_radius_m':str(m.A0_M),
                'source_calculated_onset_eV':str(source_lo),'registry_ground_ion_cost_eV':q,'finite_target_rest_threshold_eV':str(threshold),
                'first_positive_printed_source_knot_above_physical_threshold_eV':str(above),
                'minimum_dissociative_single_ion_costs_eV':{k:str(v) for k,v in fragment_costs.items()},
                'published_paper_scope_eV':[10,5000],'raw_table_extent_eV':[str(recommended[0][0]),str(recommended[-1][0])],
                'raw_high_energy_tail_admitted':False,'source_table_comparisons':len(tested),'positive_comparisons':positive,
                'printed_operand_interval_checks_passed':True,'max_table_relative_difference_at_printed_energy':str(max_print),
                'shared_kernel_comparisons':len(kernel),'shared_kernel_max_relative_error':str(max_kernel),
                'test_methods':8,'provenance_controls':10,'original_empty_model_gap_reproduced':True,'malformed_parser_cases':8,'kernel_replay_identical':True,
                'branch_assignment':'unresolved: source supplies total single ionization, not OH+ versus fragments',
                'production_qualified':False,'physical_model_validated':False,'network_complete':False,'native_adopted':False,'passed':True}
    (out/'report.json').write_text(json.dumps(report,sort_keys=True,indent=2)+'\n');verify(pins)
    receipt={'schema':'oh-primary-beb-receipt-v1.1','job':job,'allocated_cpu':sorted(os.sched_getaffinity(0)),
             'elapsed_seconds':time.monotonic()-start,'inputs_sha256':sha(ROOT/'INPUTS.json'),'source_identity':pins['source_identity'],
             'binary_sha256':sha(ROOT/'build/kernel'),'local_operand_header_sha256':sha(ROOT/'build/operands.hpp'),
             'local_kernel_values_sha256':sha(ROOT/'build/kernel-values.txt'),'files':{p.name:sha(p) for p in sorted(out.iterdir())},
             'frozen_python_sha256':pins['frozen_python_sha256'],'passed':True,'production_qualified':False,'physical_model_validated':False,
             'network_complete':False,'native_adopted':False}
    (out/'receipt.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n')
    print('PASS primary OH total single-ionization BEB checks; event branches and production remain refused')
if __name__=='__main__':main()
