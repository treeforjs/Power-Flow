"""Independent tensor/source reconstruction for centrally run native-v4 tests.

No native build is started here. The actual-source fixture is bound below once
the v4 demo API is available. Helpers deliberately avoid the reported combined
energy residual and the thermoelastic-energy implementation under review.
"""
from decimal import Decimal, localcontext
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path

import numpy as np
import pytest


QE = 1.602176634e-19
ME = 9.1093837139e-31
EPS0 = 8.8541878188e-12
LIGHT = 299792458.
EPS = np.finfo(float).eps


def _kinetic(proper_velocity, mass, weight):
    u = np.asarray(proper_velocity, dtype=float)
    u2 = float(u @ u)
    return mass*weight*u2/(math.sqrt(1.+u2/LIGHT**2)+1.)


def _thermal_state(case, nodal_temperature):
    """Independent analytic solid h, native corner weights and cell inversion."""
    t = np.asarray(nodal_temperature, dtype=float)
    assert np.all((t >= 298.2) & (t <= 773.15))
    rho = case.caloric.rho
    h = rho*(t-300.)*(472.+.068*(t+300.)-2820000./(300.*t))
    mean_h = (h[:-1, :-1]+h[:-1, 1:]+h[1:, :-1]+h[1:, 1:])/4
    xy_volume = (np.diff(case.native.y_edges_m)[:, None]
                 * np.diff(case.native.x_edges_m)[None, :]*case.native.height_m)
    H = math.fsum((mean_h*xy_volume).ravel())
    cell_t_xy = case.caloric.temperature_K(mean_h)
    cell_t = np.broadcast_to(cell_t_xy, case.mesh.cell_shape).ravel()
    return H, cell_t


def _thermal_increment(case, t0, t1):
    """Stable exact integral, independent of native quadrature delta_h."""
    a = np.asarray(t0, dtype=np.longdouble)
    b = np.asarray(t1, dtype=np.longdouble)
    change = np.longdouble(case.caloric.rho)*(b-a)*(472.+np.longdouble('.068')*(a+b)-2820000./(a*b))
    mean = (change[:-1, :-1]+change[:-1, 1:]+change[1:, :-1]+change[1:, 1:])/4
    volume = (np.diff(case.native.y_edges_m)[:, None]
              * np.diff(case.native.x_edges_m)[None, :]*case.native.height_m)
    return float(np.sum(mean*volume, dtype=np.longdouble))


def _solid_energy(case, u, velocity, temperature):
    """Three-by-three tensor Lamé law; no Voigt/helper energy reuse."""
    mesh, m = case.mesh, case.material
    t = np.asarray(temperature, dtype=float)
    gradient = np.einsum('eai,eqaj->eqij', np.asarray(u)[mesh.cells], mesh.shape_gradients_m_inv)
    strain = .5*(gradient+np.swapaxes(gradient, -1, -2))
    j = np.searchsorted(m.temperature_K, t, side='right')-1
    k = np.searchsorted(m.modulus_temperature_K, t, side='right')-1
    assert np.all((j >= 0) & (j < len(m.temperature_K)-1))
    assert np.all((k >= 0) & (k < len(m.modulus_temperature_K)-1))
    assert not np.any(t[:, None] == m.temperature_K)
    assert not np.any(t[:, None] == m.modulus_temperature_K)
    stretch_ref = 1.+np.interp(case.tref, m.temperature_K, m.relative_length_change)
    s = np.interp(t, m.temperature_K, m.relative_length_change)
    s_ref = stretch_ref-1.
    theta = (s-s_ref)/stretch_ref
    theta_t = (np.diff(m.relative_length_change)/np.diff(m.temperature_K))[j]/stretch_ref
    E = np.interp(t, m.modulus_temperature_K, m.young_modulus_Pa)
    E_t = (np.diff(m.young_modulus_Pa)/np.diff(m.modulus_temperature_K))[k]
    mu = E/(2*(1+case.nu));lam = E*case.nu/((1+case.nu)*(1-2*case.nu))
    mu_t = E_t/(2*(1+case.nu));lam_t = E_t*case.nu/((1+case.nu)*(1-2*case.nu))
    elastic = strain-theta[:, None, None, None]*np.eye(3)
    trace = np.trace(elastic, axis1=-2, axis2=-1)
    norm2 = np.sum(elastic*elastic, axis=(-1, -2))
    phi = mu[:, None]*norm2+.5*lam[:, None]*trace**2
    phi_t = (mu_t[:, None]*norm2+.5*lam_t[:, None]*trace**2
             -(2*mu+3*lam)[:, None]*theta_t[:, None]*trace)
    Uel = math.fsum(((phi-t[:, None]*phi_t)*mesh.gauss_weights_m3).ravel())
    lumped_mass = np.zeros(len(mesh.nodes_m))
    np.add.at(lumped_mass, mesh.cells.ravel(), np.repeat(case.caloric.rho*mesh.reference_volumes_m3/8, 8))
    Ks = .5*math.fsum(lumped_mass*np.sum(np.asarray(velocity)**2, axis=1))
    return Uel, Ks


def _rows(source_state, primary, alias=None):
    if primary in source_state:
        return source_state[primary]
    if alias is not None and alias in source_state:
        return source_state[alias]
    raise AssertionError(f'missing actual native source rows: {primary}')


def _vec(row, prefix, suffix=''):
    return np.array([row[f'{prefix}{axis}{suffix}'] for axis in 'xyz'])


def _electron_energy(source_state):
    kinetic, reference, image = [], [], []
    for layer in _rows(source_state, 'layers'):
        if layer['status'] == 2:
            assert layer['current_weighted_K_J'] == 0
            continue
        assert layer['status'] in (0, 1)
        K = _kinetic(_vec(layer, 'proper_v', '_m_s'), ME, layer['weight'])
        assert K == pytest.approx(layer['current_weighted_K_J'], rel=5e-14, abs=0)
        kinetic.append(K)
        reference.append(layer['weight']*QE*layer['work_function_eV'])
        if layer['status'] == 0:
            distance = float((_vec(layer, 'position_', '_m')-_vec(layer, 'origin_', '_m')) @ _vec(layer, 'normal_'))
            assert distance > 0
            image.append(-layer['weight']*QE**2/(16*math.pi*EPS0*distance))
    return np.array([math.fsum(kinetic), math.fsum(reference), math.fsum(image)])


def _image_coordinate_bound(source_state):
    """Absolute global-coordinate export can round a microscopic wall distance."""
    result = 0.
    for layer in _rows(source_state, 'layers'):
        if layer['status'] != 0:
            continue
        point, origin = _vec(layer, 'position_', '_m'), _vec(layer, 'origin_', '_m')
        normal = _vec(layer, 'normal_')
        distance = float((point-origin) @ normal)
        uncertainty = 64*EPS*float((np.abs(point)+np.abs(origin)) @ np.abs(normal))
        result += layer['weight']*QE**2/(16*math.pi*EPS0*distance)*(uncertainty/distance)
    return result


def _gas_energy(source_state):
    kinetic, internal = [], []
    for packet in _rows(source_state, 'gas_packets', 'gas'):
        assert packet['charge_C'] == packet['reference_J'] == 0.
        velocity = _vec(packet, 'proper_v', '_m_s')
        K = _kinetic(velocity, packet['mass_kg'], packet['weight'])
        assert K == pytest.approx(packet['kinetic_J'], rel=5e-14, abs=0)
        np.testing.assert_allclose(_vec(packet, 'support_impulse_', '_N_s'),
            -packet['mass_kg']*packet['weight']*velocity, rtol=5e-15, atol=0)
        kinetic.append(K);internal.append(packet['internal_J'])
    return np.array([math.fsum(kinetic), math.fsum(internal)])


def _adsorbed_increment(case, source0, source1):
    """Decimal hi/lo count difference avoids erasing sub-ULP water release."""
    keys = ('adsorbed_N_hi','adsorbed_N_lo','gas_N_hi','gas_N_lo','adsorbed_reference_energy_J')
    a = np.array([[row[k] for k in keys] for row in source0['reservoirs']])
    b = np.array([[row[k] for k in keys] for row in source1['reservoirs']])
    patches = case.native.binding['patches']
    values = []
    with localcontext() as ctx:
        ctx.prec = 60
        D = Decimal.from_float
        for index, patch in enumerate(patches):
            n0 = D(a[index, 0])+D(a[index, 1])
            n1 = D(b[index, 0])+D(b[index, 1])
            ng0 = D(a[index, 2])+D(a[index, 3])
            ng1 = D(b[index, 2])+D(b[index, 3])
            released = n0-n1
            mismatch = released-(ng1-ng0)
            scale = abs(released)+abs(ng1-ng0)
            assert abs(mismatch) <= D(64*EPS)*scale
            if not patch['water_enabled']:
                assert released == 0
                continue
            width = case.native.binding['width_m'] if patch['face'] in (0, 1) else case.native.binding['thickness_m']
            area = (patch['end']-patch['begin'])*width*case.native.height_m
            nml = D(1.3e19)*D(area)
            delta = D(.87)*D(QE)*released*(1-D(.5)*(n0+n1)/(2*nml))
            values.append(float(delta))
    return math.fsum(values)


def audit_trial(case, trial, state0, u0, v0):
    """Audit a real prepared native trial from states, tensors, paths and packets."""
    r = trial.receipt;b = r['native_budget']
    source0, source1 = r['native_source_state0'], r['native_source_state1']
    t0 = state0['temperature_K'];t1 = np.asarray(r['native_end_temperature_K'])
    H0, cell0 = _thermal_state(case, t0);H1, cell1 = _thermal_state(case, t1)
    assert H0 == pytest.approx(state0['stored_energy_J'], rel=1e-12, abs=b['roundoff_bound_J'])
    assert H1 == pytest.approx(b['trial_reference_energy_J'], rel=1e-12, abs=b['roundoff_bound_J'])
    np.testing.assert_allclose(cell1, r['end_temperature_K'], rtol=0, atol=2e-12)
    dH = _thermal_increment(case, t0, t1)
    U0, K0 = _solid_energy(case, u0, v0, cell0)
    U1, K1 = _solid_energy(case, trial.end_displacement_m, trial.end_velocity_m_s, cell1)
    de = _electron_energy(source1)-_electron_energy(source0)
    gas0, gas1 = _rows(source0, 'gas_packets', 'gas'), _rows(source1, 'gas_packets', 'gas')
    assert gas1[:len(gas0)] == gas0
    _gas_energy(source0)  # Old accepted rows retain their own physical energy identities.
    dg = _gas_energy({'gas_packets':gas1[len(gas0):]})
    dads = _adsorbed_increment(case, source0, source1)
    old_paths = len(_rows(source0, 'paths'))
    field_work = [];path_coordinate_bound = 0.
    for path in _rows(source1, 'paths')[old_paths:]:
        before, after = _vec(path, 'from_', '_m'), _vec(path, 'to_', '_m')
        displacement = after-before
        field = _vec(r['source_stage'], 'E_', '_V_m')
        expected = path['charge_C']*float(field @ displacement)
        coordinate_bound = 64*EPS*abs(path['charge_C'])*float(np.abs(field) @ (np.abs(before)+np.abs(after)))
        bound = 2e-10*(abs(expected)+abs(path['field_work_J']))+coordinate_bound
        path_coordinate_bound += coordinate_bound
        assert abs(expected-path['field_work_J']) <= bound
        assert r['interval_s'][0] <= path['start_time_s'] <= path['end_time_s'] <= r['interval_s'][1]
        field_work.append(expected)
    W = math.fsum(field_work)
    source_scale = sum(abs(x) for x in (*de, *dg, dads, W, b['source_wall_J']))
    source_bound = (3e-10*source_scale+_image_coordinate_bound(source0)
                    +_image_coordinate_bound(source1)+path_coordinate_bound)
    np.testing.assert_allclose(de, [b['source_particle_K_change_J'], b['source_electron_reference_change_J'],
                                   b['source_image_energy_change_J']], rtol=3e-10, atol=source_bound)
    np.testing.assert_allclose(dg, [b['neutral_K_J'], b['neutral_internal_J']], rtol=3e-13, atol=0)
    assert dads == pytest.approx(b['adsorbed_reference_change_J'], rel=3e-12, abs=0)
    assert W == pytest.approx(b['source_field_work_J'], rel=3e-10, abs=source_bound)
    ledger = {key:np.asarray(value, dtype=np.longdouble) for key,value in b['nodal_ledger'].items()}
    lhs = (ledger['delta_H_J']-ledger['mechanical_load_J']-ledger['electron_load_J']
           -ledger['water_load_J']+ledger['conduction_out_J'])
    activity = sum(np.abs(ledger[key]) for key in ('delta_H_J','mechanical_load_J','electron_load_J',
                                                   'water_load_J','conduction_out_J','residual_J'))
    assert np.all(np.abs(lhs-ledger['residual_J']) <= 512*EPS*activity)
    assert np.max(np.abs(lhs)) <= b['roundoff_bound_J']+512*EPS*float(np.max(activity))
    assert math.fsum(ledger['electron_load_J'].ravel()) == pytest.approx(b['source_wall_J'], rel=3e-13, abs=0)
    assert math.fsum(ledger['water_load_J'].ravel()) == pytest.approx(-math.fsum([*dg, dads]), rel=3e-12, abs=0)
    assert math.fsum(ledger['mechanical_load_J'].ravel()) == pytest.approx(b['mechanical_exchange_J'],
        rel=3e-13, abs=512*EPS*b['mechanical_absolute_activity_J'])
    external_error = ledger['external_load_J']-ledger['electron_load_J']-ledger['water_load_J']
    external_scale = np.abs(ledger['external_load_J'])+np.abs(ledger['electron_load_J'])+np.abs(ledger['water_load_J'])
    assert np.all(np.abs(external_error) <= 64*EPS*external_scale)
    residual = math.fsum([dH, U1-U0, K1-K0, *de, *dg, dads, -W])
    assert abs(residual) <= r['energy_tolerance_J']+source_bound
    assert b['joule_J'] == b['authored_heat_J'] == b['gas_reference_J'] == 0.
    return {'independent_combined_residual_J':residual, 'independent_delta_H_J':dH,
            'independent_delta_Uel_J':U1-U0, 'independent_delta_Ksolid_J':K1-K0,
            'independent_electron_changes_J':de.tolist(), 'independent_water_changes_J':[*dg, dads],
            'independent_path_work_J':W}


@pytest.fixture
def actual_native_factory():
    executable = os.environ.get('SS304_NATIVE_V4_EXECUTABLE')
    if not executable or not Path(executable).is_file():
        pytest.skip('requires separately built sealed native-v4 worker; no build is attempted')
    assert sorted(os.sched_getaffinity(0)) == [7], 'central frozen CPU7 lease required'
    assert os.getpriority(os.PRIO_PROCESS, 0) >= 15
    path = Path(__file__).resolve().parents[1]/'tools/run_native_thermoelastic_v4.py'
    spec = importlib.util.spec_from_file_location('independent_native_v4_demo', path)
    demo = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(demo)
    cases = []
    def create(**kwargs):
        case = demo.make_case(Path(executable), **kwargs)
        cases.append(case)
        return case
    yield create, demo.inputs
    for case in cases:
        case.native.close()


@pytest.mark.parametrize('velocity_sign',[1., -1.])
def test_actual_sources_independent_tensor_and_packet_energy(actual_native_factory, velocity_sign):
    create, inputs = actual_native_factory
    case = create(temperature_K=700., velocity_sign=velocity_sign)
    original_digest = case.state_digest()
    packet_counts = []
    for step, dt in enumerate((1e-15, 1e-12, 1e-10)):
        state0 = case.native.state();u0 = case.u.copy();v0 = case.v.copy()
        trial = case.prepare(dt, transaction_id=f'independent-source-energy-{step}',
                             **inputs(case, dt, begin=(step == 0)))
        r = trial.receipt;b = r['native_budget']
        audit_trial(case, trial, state0, u0, v0)
        assert r['source_candidate_held_at_accepted_start'] is True
        assert r['moving_source_surface_coupled'] is False
        assert r['field_solve_performed'] is False
        source0, source1 = r['native_source_state0'], r['native_source_state1']
        assert all(packet['reference_J'] == 0 for packet in source1['gas_packets'])
        assert len(source1['gas_packets'])-len(source0['gas_packets']) == b['new_gas_packet_count'] > 0
        packet_counts.append(len(source1['gas_packets']))
        assert b['adsorbed_reference_change_J'] > 0
        if step == 0:
            assert case.state_digest() == original_digest
            assert len(source1['layers']) == 1
            assert source1['layers'][0]['status'] == 0  # Actual partial layer, not a handoff fixture.
            assert b['physical_electron_number_change'] == 1000.
            assert b['source_electron_reference_change_J'] > 0
        else:
            # Continuation and tail restoration do not repeat birth/reference cost.
            assert b['physical_electron_number_change'] == 0
            assert b['source_electron_reference_change_J'] == 0
            assert b['source_wall_J'] == 0
            assert len(source1['layers']) == len(source0['layers']) == 1
        case.commit(trial)
    assert packet_counts == sorted(set(packet_counts))


@pytest.mark.parametrize('admission',['material_phase_admitted','planar_barrier_admitted','electron_supply_admitted'])
def test_actual_source_requires_each_explicit_admission(actual_native_factory, admission):
    create, inputs = actual_native_factory
    case = create(temperature_K=700.)
    before = case.state_digest()
    request = inputs(case, 1e-15, begin=True)
    request['source_stage'][admission] = False
    canonical = {k:v for k,v in request['source_stage'].items() if k != 'field_sha256'}
    request['source_stage']['field_sha256'] = hashlib.sha256(json.dumps(
        canonical, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
    # Deliberately change an actual serialized admission, not an unused side field.
    with pytest.raises(ValueError):
        case.prepare(1e-15, transaction_id='missing-source-admission', **request)
    assert case.state_digest() == before


def test_reviewed_cold_water_rejected_and_conditional_choice_explicit(actual_native_factory):
    create, inputs = actual_native_factory
    # Construction may reject immediately; otherwise prepare must fail before acceptance.
    reviewed = None
    try:
        reviewed = create(temperature_K=299.5, water_range=0)
    except ValueError as error:
        assert 'water' in str(error).lower() or 'source' in str(error).lower()
    if reviewed is not None:
        before = reviewed.state_digest()
        with pytest.raises(ValueError, match='water|source'):
            reviewed.prepare(1e-15, transaction_id='reviewed-water-cold', **inputs(reviewed, 1e-15))
        assert reviewed.state_digest() == before
    conditional = create(temperature_K=299.5, water_range=1)
    state0 = conditional.native.state();u0 = conditional.u.copy();v0 = conditional.v.copy()
    trial = conditional.prepare(1e-15, transaction_id='explicit-conditional-cold-water',
                                **inputs(conditional, 1e-15))
    audit_trial(conditional, trial, state0, u0, v0)
    r = trial.receipt
    assert conditional.native.binding['patches'][0]['water_range'] == 1
    assert r['native_budget']['new_gas_packet_count'] > 0
    assert min(r['end_temperature_K']) < 300
    assert all(p['reference_J'] == 0 for p in r['native_source_state1']['gas_packets'])
