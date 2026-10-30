import csv, json, math, os
from pathlib import Path
from phase_map import ObliqueSNCoordinates
import saddle_modes as sm

assert os.sched_getaffinity(0)=={7} and os.getpriority(os.PRIO_PROCESS,0)>=15
out=Path('/home/treefor/sim-workspace/powerflow/outputs/pic-oblique-phase-map-20261029')
inp=Path('/home/treefor/sim-workspace/powerflow/outputs/pic-emission-applicability-v2-20261029/frozen-probes-with-historical-temperature.csv')
rows=list(csv.DictReader(inp.open()))
checks=0;max_inverse=max_energy=max_pair_energy=max_current_return=max_bound_violation=0.
max_normal_only_error=0.;max_nonlocal=0.;actual=0;total_pairs=0
def test(b,s):
    global checks
    assert b,s
    checks+=1
def relative(a,b):return max(abs(x-y) for x,y in zip(a,b))/max(1.,max(abs(x) for x in a))
with (out/'actual-regime-coordinate-tests.csv').open('w',newline='') as f:
    writer=csv.writer(f);writer.writerow(['legacy_row','T_K','F_V_m','Bn_T','Bt_T','SN_saddle_x_m','unstable_frequency_ratio','stable_frequency_ratio','mode_point','physical_x_m','maximum_inverse_error','full_SN_energy_error_eV','kinematic_pair_energy_residual_eV','kinematic_pair_displacement_z_m','physical_scattering_phase_qualified'])
    for index,row in enumerate(rows):
        n=[float(row[k]) for k in ('nx','ny','nz')];E=[float(row[k]) for k in ('Ex_V_m','Ey_V_m','Ez_V_m')];B=[float(row[k]) for k in ('Bx_T','By_T','Bz_T')]
        F=-sum(a*b for a,b in zip(n,E))
        if F<=0:continue
        actual+=1;Bn=sum(a*b for a,b in zip(n,B));Bt=math.sqrt(sum((a-Bn*b)**2 for a,b in zip(B,n)))
        # Coordinates use the normal E component as a deliberately isolated planar
        # model. The actual nonzero tangential E is NOT marked admitted below.
        for signed in (1.,-1.):
            m=ObliqueSNCoordinates(F,signed*Bn,Bt)
            for point in range(4):
                Z=(.02*(-1)**point,.015*(-1)**(point//2),.03,.2*(-1)**point,.13,.04)
                p=m.physical(Z);back=m.canonical(p);ie=relative(Z,back);max_inverse=max(max_inverse,ie);test(ie<1e-10,'physical/canonical inverse')
                ee=abs(m.physical_energy_eV(p)-m.transformed_full_SN_energy_eV(Z));max_energy=max(max_energy,ee);test(ee<1e-10,'full SN Hamiltonian, including anharmonic remainder')
                paired,receipt=m.reciprocal_kinematic_pair(p,-sm.qe,tangential_E_V_m=(0.,0.),uniform_static_field_admitted=True)
                original,reverse=m.reciprocal_kinematic_pair(paired,-sm.qe,tangential_E_V_m=(0.,0.),uniform_static_field_admitted=True)
                inv=relative(p.position_m,original.position_m);invp=relative(tuple(v/m.P for v in p.momentum_kg_m_s),tuple(v/m.P for v in original.momentum_kg_m_s));test(inv<1e-10 and invp<1e-10,'candidate pairing is an involution')
                er=abs(receipt['energy_residual_J_per_electron'])/sm.qe;max_pair_energy=max(max_pair_energy,er);test(er<1e-10,'full SN product energy receipt')
                cr=max(abs(a+b) for a,b in zip(receipt['charge_displacement_C_m'],reverse['charge_displacement_C_m']))/sm.qe;max_current_return=max(max_current_return,cr);test(cr<1e-15,'paired candidate current displacement cancels once')
                test(receipt['wall_heat_J']==0 and receipt['electron_reference_change_J']==0,'kinematic pairing cannot absorb at metal')
                dz=abs(receipt['displacement_m'][2]);max_nonlocal=max(max_nonlocal,dz)
                # Replacing this by a flip only along the material normal generally
                # changes the stable mode too, and fails to invert the reaction pair.
                wrong=list(p.momentum_kg_m_s);wrong[0]=-wrong[0]
                from phase_map import MechanicalPhase
                wrongZ=m.canonical(MechanicalPhase(p.position_m,tuple(wrong)))
                stable_error=abs(wrongZ[4]-Z[4]);max_normal_only_error=max(max_normal_only_error,stable_error)
                writer.writerow([index,row['historical_T_K'],F,signed*Bn,Bt,m.xs,m.modes['unstable'],m.modes['stable'],point,p.position_m[0],ie,ee,er,receipt['displacement_m'][2],False]);total_pairs+=1
            for a in (.1,.25,.5):
                bound=m.uniform_potential_error_bound_eV(a)
                for k in range(33):
                    x=m.xs*(1-a+2*a*k/32);xi=x-m.xs
                    exact=m.SN_potential_eV(x)-(m.Us-m.F*xi*xi/m.xs)
                    rem=m.anharmonic_remainder_eV(x)
                    test(abs(exact-rem)<2e-14,'exact cubic-over-linear SN remainder')
                    violation=abs(rem)-bound;max_bound_violation=max(max_bound_violation,violation);test(violation<1e-14,'uniform local SN potential bound')
            try:m.reciprocal_kinematic_pair(p,-sm.qe,tangential_E_V_m=(1.,0.),uniform_static_field_admitted=True)
            except ValueError:checks+=1
            else:raise AssertionError('nonzero tangential E was silently admitted')
            try:m.emit()
            except RuntimeError:checks+=1
            else:raise AssertionError('unqualified physical emitter was enabled')
test(max_nonlocal>0 and max_normal_only_error>1e-5,'oblique map genuinely requires correlated momentum and displacement')
report={'schema':'exact-oblique-full-SN-coordinate-and-kinematic-reciprocity-audit-v1','passed':True,'checks':checks,'actual_positive_F_regimes':actual,'paired_coordinate_cases':total_pairs,'both_signed_Bn_tested':True,'maximum_inverse_dimensionless_error':max_inverse,'maximum_full_SN_energy_error_eV':max_energy,'maximum_kinematic_pair_energy_residual_eV':max_pair_energy,'maximum_reverse_current_displacement_residual_m_per_electron':max_current_return,'maximum_potential_bound_violation_eV':max_bound_violation,'maximum_material_normal_only_flip_stable_mode_error':max_normal_only_error,'maximum_candidate_pair_lateral_displacement_m':max_nonlocal,'exact_full_SN_potential_retained':True,'nominal_rate_or_quadratic_prefactor_deployed':False,'physical_quantum_scattering_phase_qualified':False,'positive_quantum_product_phase_measure_qualified':False,'full_SN_transmission_current_or_heat_error_bound':None,'native_product_source_admitted':False,'actual_tangential_E_admitted':False,'uniform_planar_model_only':True,'CPU':7,'nice':15,'threads':1,'GPU_work':False,'production_action':False}
(out/'qualification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report))
