import csv
from dataclasses import replace
import hashlib
import json

import numpy as np
import pytest

from mhdlab.eos_tables import BoundedEOSTable, COLUMNS, SCHEMA, UNITS, load_eos_table


def authored(tmp_path):
    rho=np.array([100.,1000.,8000.]);temp=np.array([300.,1000.,5000.])
    pressure=rho[:,None]*100.*temp[None,:]
    energy=np.broadcast_to(250.*(temp-1000.),pressure.shape).copy()
    path=tmp_path/'eos.csv'
    with path.open('w',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(COLUMNS)
        for i,r in enumerate(rho):
            for j,t in enumerate(temp): writer.writerow([r,t,pressure[i,j],energy[i,j]])
    metadata={'schema':SCHEMA,'units':UNITS,'parameter_status':'manufactured',
        'material_grade':'authored ideal gas','table_version':'authored-v1',
        'source_reference':'Authored ideal-gas fixture, R=100 J/(kg K), Cv=250 J/(kg K)',
        'energy_reference':'u=0 at1000K','phase_policy':'single_phase','phase_reference':'Authored ideal gas',
        'thermodynamic_consistency_reference':'Analytic ideal gas P=rhoRT,u=Cv(T-1000K)',
        'interpolation':'bilinear_in_density_and_temperature','extrapolation':'error',
        'csv_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
    sidecar=tmp_path/'eos.json';sidecar.write_text(json.dumps(metadata))
    return path,sidecar


def test_bilinear_evaluation_and_energy_inverse_preserve_negative_reference_energy(tmp_path):
    csv_path,metadata=authored(tmp_path)
    table=load_eos_table(csv_path,metadata,expected_material_grade='authored ideal gas')
    for rho,temp in [(100.,300.),(550.,650.),(8000.,5000.)]:
        values=table.evaluate(rho,temp)
        assert values.pressure_pa==pytest.approx(rho*100*temp,rel=1e-15)
        assert values.internal_energy_j_kg==pytest.approx(250*(temp-1000),rel=1e-15)
        assert table.temperature_from_internal_energy_k(rho,values.internal_energy_j_kg)==pytest.approx(temp,rel=1e-15)


def test_no_out_of_support_temperature_density_or_energy_is_clamped(tmp_path):
    csv_path,metadata=authored(tmp_path)
    table=load_eos_table(csv_path,metadata,expected_material_grade='authored ideal gas')
    for rho,temp in [(np.nextafter(100.,0.),300.),(100.,np.nextafter(5000.,np.inf))]:
        with pytest.raises(ValueError,match='outside'): table.evaluate(rho,temp)
    with pytest.raises(ValueError,match='outside'):
        table.temperature_from_internal_energy_k(100.,np.nextafter(-175000.,-np.inf))
    for invalid in [True,'8000',float('nan')]:
        with pytest.raises(ValueError,match='scalar'): table.evaluate(invalid,300.)


@pytest.mark.parametrize('mutation',[
    lambda m:m['units'].update(internal_energy='eV'),
    lambda m:m.update(phase_policy='unspecified'),
    lambda m:m.update(energy_reference='unknown'),
    lambda m:m.update(extrapolation='clamp'),
])
def test_required_units_phase_energy_and_bounds_semantics(tmp_path,mutation):
    path,sidecar=authored(tmp_path);m=json.loads(sidecar.read_text());mutation(m);sidecar.write_text(json.dumps(m))
    with pytest.raises(ValueError):load_eos_table(path,sidecar,expected_material_grade='authored ideal gas')


def test_source_digest_material_selection_and_rectangular_coverage_are_checked(tmp_path):
    path,sidecar=authored(tmp_path)
    with pytest.raises(ValueError,match='material'):load_eos_table(path,sidecar,expected_material_grade='SS304')
    with path.open('a') as stream:stream.write('\n')
    with pytest.raises(ValueError,match='SHA256'):load_eos_table(path,sidecar,expected_material_grade='authored ideal gas')
    path,sidecar=authored(tmp_path)
    lines=path.read_text().splitlines();path.write_text('\n'.join(lines[:-1])+'\n')
    m=json.loads(sidecar.read_text());m['csv_sha256']=hashlib.sha256(path.read_bytes()).hexdigest();sidecar.write_text(json.dumps(m))
    with pytest.raises(ValueError,match='rectangular'):load_eos_table(path,sidecar,expected_material_grade='authored ideal gas')


def test_increasing_energy_and_immutable_fields_are_required(tmp_path):
    path,sidecar=authored(tmp_path);table=load_eos_table(path,sidecar,expected_material_grade='authored ideal gas')
    energy=table.internal_energy_j_kg.copy();energy[:,1]=energy[:,0]
    with pytest.raises(ValueError,match='increase'):replace(table,internal_energy_j_kg=energy)
    with pytest.raises(ValueError):table.pressure_pa.setflags(write=True)
    with pytest.raises(TypeError):table.metadata['energy_reference']='silently changed'
    pressure=table.pressure_pa.copy();copy=replace(table,pressure_pa=pressure);pressure[:]=0
    assert copy.evaluate(100.,300.).pressure_pa==3000000.
