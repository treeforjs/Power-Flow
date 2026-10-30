"""Exact classical coordinate maps for a uniform oblique B and planar SN U(x).

No transmission law, product sampler or physical reflection kernel is supplied.
The reciprocal involution is a kinematic test of a candidate reaction pairing.
"""
import math
import sys
from dataclasses import dataclass
sys.path.insert(0,'<workspace>/outputs/pic-oblique-emission-20261029/source-backup-local-saddle-diagnostic-v1')
import saddle_modes as sm

@dataclass(frozen=True)
class MechanicalPhase:
    position_m: tuple
    momentum_kg_m_s: tuple

class ObliqueSNCoordinates:
    def __init__(self,F_V_m,B_normal_T,B_tangent_T,phi_eV=4.3):
        if not (math.isfinite(F_V_m) and F_V_m>0 and math.isfinite(B_normal_T) and B_normal_T!=0 and math.isfinite(B_tangent_T) and B_tangent_T>=0 and math.isfinite(phi_eV)):
            raise ValueError('finite extracting uniform planar fields with nonzero Bn required; tangent source remains separate')
        self.F,self.Bn,self.Bt,self.phi=F_V_m,B_normal_T,B_tangent_T,phi_eV
        self.xs=math.sqrt(sm.c/self.F)
        self.Us=self.phi-2*math.sqrt(sm.c*self.F)
        self.omega0=math.sqrt(2*sm.qe*sm.c/(sm.me*self.xs**3))
        self.L=self.xs
        self.P=sm.me*self.omega0*self.L
        self.energy_unit_eV=sm.me*self.omega0**2*self.L**2/sm.qe
        self.hbar_effective=sm.hbar/(self.P*self.L)
        self.bn,self.bt=sm.qe*self.Bn/(sm.me*self.omega0),sm.qe*self.Bt/(sm.me*self.omega0)
        self.modes=sm.modes(self.bn,self.bt)
        self.S=self.modes['basis']
        self.inverse=sm.multiply(sm.J,sm.multiply(sm.transpose(self.S),sm.J))
        self.inverse=[[-x for x in r] for r in self.inverse]

    def physical(self,Z):
        if len(Z)!=6 or not all(math.isfinite(x) for x in Z):raise ValueError('six finite canonical mode coordinates required')
        z=[sum(a*b for a,b in zip(row,Z)) for row in self.S]
        r=(self.xs+self.L*z[0],self.L*z[1],self.L*z[2])
        # Gauge Az=Bn*y-Bt*(x-xs), electron q=-e. Canonical pz is not mechanical.
        p=(self.P*z[3],self.P*z[4],self.P*(z[5]+self.bn*z[1]-self.bt*z[0]))
        if r[0]<=0:raise ValueError('phase outside the exterior SN material half-space')
        return MechanicalPhase(r,p)

    def canonical(self,phase):
        r,p=phase.position_m,phase.momentum_kg_m_s
        if r[0]<=0:raise ValueError('phase outside the exterior SN material half-space')
        q=[(r[0]-self.xs)/self.L,r[1]/self.L,r[2]/self.L]
        z=q+[p[0]/self.P,p[1]/self.P,p[2]/self.P-self.bn*q[1]+self.bt*q[0]]
        return tuple(sum(a*b for a,b in zip(row,z)) for row in self.inverse)

    def SN_potential_eV(self,x):
        if x<=0:raise ValueError('SN potential needs x>0')
        return self.phi-self.F*x-sm.c/x

    def anharmonic_remainder_eV(self,x):
        if x<=0:raise ValueError('SN potential needs x>0')
        xi=x-self.xs
        return self.F*xi**3/(self.xs*x)

    def uniform_potential_error_bound_eV(self,relative_half_width):
        a=relative_half_width
        if not 0<a<1:raise ValueError('finite local window must satisfy 0<a<1')
        return self.F*self.xs*a**3/(1-a)

    def physical_energy_eV(self,phase):
        return sum(p*p for p in phase.momentum_kg_m_s)/(2*sm.me*sm.qe)+self.SN_potential_eV(phase.position_m[0])

    def transformed_full_SN_energy_eV(self,Z):
        Qu,Qs,_,Pu,Ps,_=Z
        x=self.physical(Z).position_m[0]
        lam,s=self.modes['unstable'],self.modes['stable']
        return self.Us+self.energy_unit_eV*.5*(Pu*Pu-lam*lam*Qu*Qu+Ps*Ps+s*s*Qs*Qs)+self.anharmonic_remainder_eV(x)

    def reciprocal_kinematic_pair(self,phase,charge_C,*,tangential_E_V_m=(0.,0.),uniform_static_field_admitted=False):
        if not uniform_static_field_admitted or tangential_E_V_m!=(0.,0.):
            raise ValueError('this candidate pairing requires explicitly uniform static B and exactly normal E; no local-field relabeling')
        if not math.isfinite(charge_C):raise ValueError('finite charge required')
        Z=list(self.canonical(phase));Z[3]=-Z[3]
        paired=self.physical(Z)
        D=tuple(b-a for a,b in zip(phase.position_m,paired.position_m))
        dK=sum(b*b-a*a for a,b in zip(phase.momentum_kg_m_s,paired.momentum_kg_m_s))/(2*sm.me)
        dImage=sm.qe*sm.c*(1/phase.position_m[0]-1/paired.position_m[0])
        W=sm.qe*self.F*D[0]
        return paired,{'charge_displacement_C_m':tuple(charge_C*x for x in D),'displacement_m':D,'particle_K_change_J_per_electron':dK,'image_change_J_per_electron':dImage,'normal_field_work_J_per_electron':W,'wall_heat_J':0.,'electron_reference_change_J':0.,'energy_residual_J_per_electron':dK+dImage-W,'physical_scattering_phase_qualified':False}

    def emit(self,*args,**kwargs):
        raise RuntimeError('no full-SN channel transmission, scattering phase, quantum-to-PIC projection or bounded physical sampler is qualified')
