import csv,json,math,os
from pathlib import Path
from dividing_surface_bound import ClassicalDividingSurfaceBound
from phase_map import sm
assert os.sched_getaffinity(0)=={7} and os.getpriority(os.PRIO_PROCESS,0)>=15
out=Path('/home/treefor/sim-workspace/powerflow/outputs/pic-oblique-phase-map-20261029')
inp=Path('/home/treefor/sim-workspace/powerflow/outputs/pic-emission-applicability-v2-20261029/frozen-probes-with-historical-temperature.csv')
rows=list(csv.DictReader(inp.open()));checks=0;states=0;products=0;max_inverse=max_energy=max_reverse=max_numeric=0.
def test(b,s):
    global checks
    assert b,s
    checks+=1
with (out/'conditional-classical-dividing-surface-bounds.csv').open('w',newline='') as f:
    w=csv.writer(f);w.writerow(['legacy_row','face','historical_T_K','F_V_m','Bn_T','Bt_T','Us_eV','conditional_classical_J_lower_A_m2','conditional_classical_J_upper_A_m2','log_gross_J_A_m2','mean_Fermi_energy_eV','conditional_classical_absolute_Fermi_energy_flux_upper_W_m2','full_quantum_source_qualified'])
    for i,row in enumerate(rows):
        n=[float(row[k]) for k in ('nx','ny','nz')];E=[float(row[k]) for k in ('Ex_V_m','Ey_V_m','Ez_V_m')];B=[float(row[k]) for k in ('Bx_T','By_T','Bz_T')]
        F=-sum(a*b for a,b in zip(n,E))
        if F<=0:continue
        Bn=sum(a*b for a,b in zip(n,B));Bt=math.sqrt(sum((a-Bn*b)**2 for a,b in zip(B,n)));T=float(row['historical_T_K'])
        d=ClassicalDividingSurfaceBound(T,F,Bn,Bt);states+=1;b=d.bounds()
        test(d.Us>0 and d.gross_J>0 and d.absolute_Fermi_energy_flux_upper_W_m2>0,'finite positive classical bound')
        for u in (.001,.01,.1,.25,.5,.75,.9,.99,.999):
            p=d.crossing(u,.371,.193,physical_weight=3.);inverse=abs((1-d.total_K_tail(p.normal_K_eV+p.transverse_K_eV))-u);max_inverse=max(max_inverse,inverse);test(inverse<2e-14,'exact FD total energy inverse')
            K=sum(x*x for x in p.momentum_kg_m_s)/(2*sm.me*sm.qe);error=abs(K-(p.total_Fermi_energy_eV-d.Us));max_energy=max(max_energy,error);test(error<1e-14,'physical phase energy matches same FD current/energy measure')
            rev=d.time_reversed(p);test(rev.B_T==tuple(-x for x in p.B_T) and rev.momentum_kg_m_s==tuple(-x for x in p.momentum_kg_m_s),'Hamiltonian magnetic time reversal')
            test(p.signed_electron_number_crossing+rev.signed_electron_number_crossing==0 and p.signed_charge_crossing_C+rev.signed_charge_crossing_C==0 and p.signed_Fermi_energy_crossing_J+rev.signed_Fermi_energy_crossing_J==0,'reciprocal product/charge/total Fermi-energy crossing once')
            test(not p.material_origin_qualified and not p.native_birth_qualified,'proposal cannot silently create a native source particle');products+=1
        # Independent direct two-dimensional radial FD integral, with s=K/kT.
        if i%53==0:
            N=16384;h=64/N;I0=I1=0.
            for k in range(N+1):
                x=k*h;v=math.exp(-x)/(1+d.z*math.exp(-x));weight=1 if k in (0,N) else (4 if k%2 else 2)
                I0+=weight*x*v;I1+=weight*x*x*v
            I0*=h/3;I1*=h/3
            je=abs(I0/d.r1-1);ee=abs((d.Us+d.t*I1/I0)-d.mean_Fermi_energy_eV);max_numeric=max(max_numeric,je,ee);test(je<1e-10 and ee<1e-10,'independent radial current and energy integration')
        try:d.emit()
        except RuntimeError:checks+=1
        else:raise AssertionError('unqualified wall source enabled')
        w.writerow([i,row['name'],T,F,Bn,Bt,d.Us,0,d.gross_J,d.log_gross_J,d.mean_Fermi_energy_eV,d.absolute_Fermi_energy_flux_upper_W_m2,False])
report={'schema':'conditional-classical-full-SN-dividing-surface-current-energy-product-bound-v1','passed':True,'checks':checks,'actual_extraction_parameter_states':states,'physical_gross_crossing_proposals':products,'maximum_inverse_CDF_absolute_error':max_inverse,'maximum_phase_energy_error_eV':max_energy,'independent_radial_FD_integral_max_error':max_numeric,'current_energy_products_same_FD_measure':True,'B_reversal_reciprocal_crossings_exact':True,'retains_full_SN_saddle_energy':True,'quadratic_prefactor_or_unmagnetized_guard_resampling':False,'magnetic_dependence':'material connectivity and recrossing remain unresolved; interval is [0,gross crossing flux], not a B-independent actual emission law','bound_scope':'only declared classical continuum f(H) model; excludes quantum tunnelling/reflection/band/spin errors','actual_frozen_tangential_field_or_new_midpoint_admitted':False,'native_birth_or_wall_heat_current_admitted':False,'CPU':7,'nice':15,'threads':1,'GPU_work':False,'production_action':False}
(out/'dividing-surface-qualification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
