"""Conditional classical gross-crossing bound, not an accepted wall emitter.

Ideal continuum f(H) reservoir; full planar SN saddle; no Landau/spin/band
correction. Unknown material connection lies in [0,1] per phase. Quantum
tunnelling or quantum reflection is not bounded by this classical result.
"""
import math
from dataclasses import dataclass
import phase_map
sm=phase_map.sm

def fermi_ratio(z,power):
    # Fermi polylog sum / z; z<=exp(-5), with alternating-series remainder.
    value=0.;term=1.
    for k in range(1,128):
        value+=term/(k**power)
        term*=-z
        if abs(term)/((k+1)**power)<1e-18:return value
    raise RuntimeError('bounded alternating Fermi series did not converge')

@dataclass(frozen=True)
class GrossCrossing:
    position_m: tuple
    momentum_kg_m_s: tuple
    total_Fermi_energy_eV: float
    normal_K_eV: float
    transverse_K_eV: float
    signed_electron_number_crossing: float
    signed_charge_crossing_C: float
    signed_Fermi_energy_crossing_J: float
    direction: int
    B_T: tuple
    material_origin_qualified: bool=False
    native_birth_qualified: bool=False

class ClassicalDividingSurfaceBound:
    def __init__(self,T_K,F_V_m,B_normal_T,B_tangent_T,phi_eV=4.3):
        if not (1<=T_K<=3000 and 1<=F_V_m<=5e9 and 2<=phi_eV<=6 and all(math.isfinite(x) for x in (T_K,F_V_m,B_normal_T,B_tangent_T,phi_eV))):raise ValueError('finite explicit numerical benchmark required')
        self.T,self.F,self.Bn,self.Bt,self.phi=T_K,F_V_m,B_normal_T,B_tangent_T,phi_eV
        self.t=sm.kb*T_K/sm.qe;self.xs=math.sqrt(sm.c/F_V_m);self.Us=phi_eV-2*math.sqrt(sm.c*F_V_m)
        if self.Us/self.t<5:raise ValueError('this bounded-series implementation requires Us/kT>=5; no extrapolation')
        self.logz=-self.Us/self.t;self.z=math.exp(self.logz)
        self.r1=fermi_ratio(self.z,2);self.r2=fermi_ratio(self.z,3)
        self.log_gross_J=math.log(sm.C)+2*math.log(self.t)+self.logz+math.log(self.r1)
        self.gross_J=math.exp(self.log_gross_J)
        self.mean_Fermi_energy_eV=self.Us+2*self.t*self.r2/self.r1
        self.absolute_Fermi_energy_flux_upper_W_m2=self.gross_J*self.mean_Fermi_energy_eV

    def total_K_tail(self,K_eV):
        if K_eV<0:raise ValueError('K>=0 required')
        a=K_eV/self.t;decay=math.exp(-a);z=self.z*decay
        soft_ratio=math.log1p(z)/z if z else 1.
        return decay*(a*soft_ratio+fermi_ratio(z,2))/self.r1

    def total_K_quantile(self,u):
        if not math.isfinite(u) or not 0<u<1:raise ValueError('caller-owned open draw required')
        lo,hi=0.,8*self.t
        for _ in range(16):
            if self.total_K_tail(hi)<=1-u:break
            hi*=2
        else:raise RuntimeError('bounded energy inversion bracket failed')
        for _ in range(80):
            mid=(lo+hi)/2
            if self.total_K_tail(mid)>1-u:lo=mid
            else:hi=mid
        return (lo+hi)/2

    def crossing(self,energy_uniform,partition_uniform,azimuth_uniform,physical_weight=1.):
        if not all(math.isfinite(u) and 0<u<1 for u in (energy_uniform,partition_uniform,azimuth_uniform)) or not math.isfinite(physical_weight) or physical_weight<=0:raise ValueError('three open draws and positive physical weight required')
        K=self.total_K_quantile(energy_uniform);Kn=K*partition_uniform;Kt=K-Kn;angle=2*math.pi*azimuth_uniform
        pn=math.sqrt(2*sm.me*sm.qe*Kn);pt=math.sqrt(2*sm.me*sm.qe*Kt)
        return GrossCrossing((self.xs,0.,0.),(pn,pt*math.cos(angle),pt*math.sin(angle)),self.Us+K,Kn,Kt,physical_weight,-sm.qe*physical_weight,physical_weight*(self.Us+K)*sm.qe,1,(self.Bn,self.Bt,0.))

    @staticmethod
    def time_reversed(crossing):
        return GrossCrossing(crossing.position_m,tuple(-p for p in crossing.momentum_kg_m_s),crossing.total_Fermi_energy_eV,crossing.normal_K_eV,crossing.transverse_K_eV,-crossing.signed_electron_number_crossing,-crossing.signed_charge_crossing_C,-crossing.signed_Fermi_energy_crossing_J,-crossing.direction,tuple(-b for b in crossing.B_T))

    def bounds(self):
        return {'conditional_classical_material_connected_J_interval_A_m2':[0.,self.gross_J],'conditional_classical_absolute_Fermi_energy_flux_interval_W_m2':[0.,self.absolute_Fermi_energy_flux_upper_W_m2],'log_gross_crossing_J_A_m2':self.log_gross_J,'mean_Fermi_energy_eV':self.mean_Fermi_energy_eV,'scope':'classical, static full-SN, continuum f(H) occupation, elastic Hamiltonian motion, material connection subset of gross positive crossings; not a quantum-emission bound','material_connection_qualified':False,'full_quantum_tunnelling_or_reflection_bound':None,'accepted_wall_heat_or_birth_current':False}

    def emit(self,*args,**kwargs):
        raise RuntimeError('gross dividing-surface proposals are not wall births; material connection, return paths, quantum scattering and source ownership remain unqualified')
