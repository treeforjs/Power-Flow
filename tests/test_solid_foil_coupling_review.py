"""Independent tensor/first-law checks for the conditional cell corrector."""
import numpy as np
import pytest

from mhdlab.solid_foil_mechanics import structured_hex_mesh
from mhdlab.solid_foil_coupling import coupled_cell_heat_step
from mhdlab.thermal_expansion import load_ss304_material


class _Caloric:
    capacity = 3.8e6

    def __init__(self, rho):
        self.receipt = {"reference_density_kg_m3": rho,
                        "basis": "independent manufactured constant stress-free capacity"}

    def energy_J_m3(self, t):
        return self.capacity*(np.asarray(t)-300.)

    def capacity_J_m3_K(self, t):
        return np.zeros_like(np.asarray(t))+self.capacity

    def temperature_K(self, u):
        return np.asarray(u)/self.capacity+300.


def _setup(reference=300.):
    mesh=structured_hex_mesh([-.0025,0.,.0025],[.0025,.0027],[-.0125,0.,.0125])
    material=load_ss304_material()
    rho=material.density_at_stress_free_temperature(reference)
    caloric=_Caloric(rho)
    args=dict(material=material,stress_free_temperature_K=reference,poisson_ratio=.29,
        boundary='free',caloric=caloric,state_revision=0,mesh_revision='reviewmesh',
        material_sha256='a'*64,source_sha256='b'*64,transaction_id='review',
        proof_stress_screen_Pa=lambda t:np.full_like(t,190e6))
    return mesh,args


def _direct_total_energy(mesh, temperature, displacement, velocity, reference, nu, caloric):
    # Direct 3x3 tensor implementation, independent of Voigt/work-port helpers.
    gradient=np.einsum('eai,eqaj->eqij',displacement[mesh.cells],mesh.shape_gradients_m_inv)
    strain=(gradient+gradient.swapaxes(-1,-2))/2
    sref=16e-6*(reference-293.15)
    eth=(16e-6*(temperature-293.15)-sref)/(1+sref)
    alpha=16e-6/(1+sref)
    young=200e9-75e6*(temperature-293.15)
    elastic=strain-eth[:,None,None,None]*np.eye(3)
    trace=np.trace(elastic,axis1=-2,axis2=-1)
    norm=np.sum(elastic*elastic,axis=(-1,-2))
    mu=young/(2*(1+nu)); lam=young*nu/((1+nu)*(1-2*nu))
    mu_t=-75e6/(2*(1+nu)); lam_t=-75e6*nu/((1+nu)*(1-2*nu))
    phi=mu[:,None]*norm+.5*lam[:,None]*trace**2
    phi_t=mu_t*norm+.5*lam_t*trace**2-alpha*(3*lam+2*mu)[:,None]*trace
    internal=phi-temperature[:,None]*phi_t
    internal_total=np.sum(internal*mesh.gauss_weights_m3)
    base=np.sum(caloric.energy_J_m3(temperature)*mesh.reference_volumes_m3)
    rho=caloric.receipt['reference_density_kg_m3']
    nodal_mass=np.zeros(len(mesh.nodes_m))
    for cell,volume in zip(mesh.cells,mesh.reference_volumes_m3):
        nodal_mass[cell]+=rho*volume/8
    kinetic=.5*np.sum(nodal_mass[:,None]*velocity*velocity)
    return base+internal_total+kinetic


def test_global_first_law_from_independent_tensor_internal_energy():
    mesh,args=_setup()
    t0=np.array([300.,301.,302.,303.])
    heat=mesh.reference_volumes_m3*args['caloric'].capacity*np.array([1.,2.,3.,4.])
    result=coupled_cell_heat_step(mesh,t0,heat,dt_s=2e-9,**args)
    mech=result['mechanics']; packet=result['exchange']
    total=[]
    for i in (0,-1):
        total.append(_direct_total_energy(mesh,mech['cell_temperature_K'][i],
            mech['displacement_m'][i],mech['velocity_m_s'][i],300.,.29,args['caloric']))
    assert total[1]-total[0] == pytest.approx(float(heat.sum()),rel=2e-10,abs=1e-14)
    assert abs(packet['coupled_energy_residual_J']) < 1e-10*heat.sum()


def test_zero_external_heat_affine_initial_expansion_cools():
    mesh,args=_setup(reference=320.)
    initial_velocity=1000.*(mesh.nodes_m-mesh.nodes_m.mean(axis=0))
    result=coupled_cell_heat_step(mesh,np.full(len(mesh.cells),320.),np.zeros(len(mesh.cells)),
        dt_s=1e-9,initial_velocity_m_s=initial_velocity,**args)
    mech=result['mechanics']; packet=result['exchange']
    assert np.all(np.asarray(packet['end_temperature_K']) < 320.)
    total=[]
    for i in (0,-1):
        total.append(_direct_total_energy(mesh,mech['cell_temperature_K'][i],
            mech['displacement_m'][i],mech['velocity_m_s'][i],320.,.29,args['caloric']))
    assert total[1] == pytest.approx(total[0],rel=1e-12,abs=1e-14)
    assert sum(packet['cell_thermoelastic_debit_J']) > 0
