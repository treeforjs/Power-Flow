"""New ledger export from a preserved two-state launch-packet association."""
import csv
from fractions import Fraction
import hashlib
import importlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import pytest

from mhdlab.constants import KB,QE
from mhdlab.flash_surface import FlashSurfaceHistory
from mhdlab.planar_release import PlanarSurfaceHistory
from mhdlab.temkin_planar import evolve_temkin_planar_release
from mhdlab.temkin_planar_io import write_temkin_planar_release_history
from mhdlab.temkin_release import TemkinReservoir

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/'tools/export_packet_source_ledger.py'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path,value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    return path


def command(directory,output):
    return [str(SCRIPT),'--packet-export',str(directory),'--output-dir',str(output)]


def run(directory,output):
    return subprocess.run([sys.executable,'-W','error',*command(directory,output)],cwd=ROOT,
                          text=True,capture_output=True,timeout=60)


@pytest.fixture(scope='module')
def preserved_packets(tmp_path_factory):
    directory=tmp_path_factory.mktemp('new-packet-ledger-source')
    provenance='Manufactured two-parent ledger reference; no physical water parameters'
    base=FlashSurfaceHistory([0.,1e-7],[0.,1.],[[300.],[300.]],'synthetic',
                            'new-packet-ledger-fixture','Manufactured prescribed temperature',
                            'Manufactured stationary face')
    surface=PlanarSurfaceHistory([0.,1e-7],[-.003,-.001,.002],[.001,.004],
                                [[[400.,900.]],[[400.,900.]]],base.attributes)
    pools=[TemkinReservoir('parent_a','M:a',1e7,0.,.2,.04,0.,1e18,1e18,'Manufactured constant rate','Manufactured supplied costs'),
           TemkinReservoir('parent_b','M:b',2e7,0.,.5,.12,0.,1e18,1e18,'Manufactured constant rate','Manufactured supplied costs')]
    initial=np.array([[[1.2e17,2.4e17]],[[.8e17,1.5e17]]])
    history=evolve_temkin_planar_release(surface,initial,pools,parameter_provenance=provenance)
    source=directory/'source.h5'
    write_temkin_planar_release_history(source,history)
    summary=write_json(directory/'summary.json',{'schema':'mhdlab-temkin-planar-run-v1',
        'parameter_status':'manufactured','parameter_provenance':provenance,'output_sha256':{source.name:sha(source)}})
    config=write_json(directory/'config.json',{'schema':'mhdlab-temkin-windowed-packet-config-v1',
        'parameter_status':'manufactured','conditioning_provenance':'Manufactured disjoint source selections',
        'mass_status':'manufactured','mass_provenance':'Two supplied manufactured parent masses',
        'state_masses_kg':{'M:a':3e-26,'M:b':6e-26},'selections':[
            {'interval_index':0,'start_s':1e-8,'stop_s':3e-8},
            {'interval_index':0,'start_s':6e-8,'stop_s':9e-8}],
        'log_mass_tolerance':1e-7,'face_normal_coordinate_m':0.,'outward_normal_sign':1,
        'packet_count':32,'seed':19431})
    result=directory/'preserved'
    process=subprocess.run([sys.executable,'-W','error',str(ROOT/'tools/sample_temkin_source_windows.py'),
        '--source',str(source),'--source-run-summary',str(summary),'--config',str(config),'--output-dir',str(result)],
        cwd=ROOT,text=True,capture_output=True,timeout=60)
    assert process.returncode==0,process.stderr
    return result


@pytest.fixture(scope='module')
def exported(preserved_packets,tmp_path_factory):
    output=tmp_path_factory.mktemp('new-packet-ledger-result')/'result'
    process=run(preserved_packets,output)
    assert process.returncode==0,process.stderr
    return output


def test_parent_counts_and_all_energy_channels_match_independent_rates(exported,preserved_packets):
    report=json.loads((exported/'report.json').read_text())
    ledger=json.loads((exported/'ledger.json').read_text())
    assert report['schema']=='mhdlab-packet-source-ledger-export-v1'
    assert report['parameter_status']=='manufactured'
    for flag in ('host_budget_accepted','source_evolution_repeated','packet_sampling_repeated','source_moments_recomputed'):
        assert report[flag] is False
    expected={name:0. for name in ('released_number','translational_energy_j','internal_energy_j','binding_energy_j','total_wall_energy_debit_j')}
    for pool,rate,initial,internal,removal in ((0,1e7,[1.2e17,2.4e17],.04,.2),(1,2e7,[.8e17,1.5e17],.12,.5)):
        for density,width,temperature in zip(initial,[.002,.003],[400.,900.]):
            for start,stop in ((1e-8,3e-8),(6e-8,9e-8)):
                number=density*math.exp(-rate*start)*(-math.expm1(-rate*(stop-start)))*width*.003
                costs={'released_number':1.,'translational_energy_j':2*KB*temperature,
                       'internal_energy_j':internal*QE,'binding_energy_j':removal*QE,
                       'total_wall_energy_debit_j':2*KB*temperature+(internal+removal)*QE}
                for name,cost in costs.items():expected[name]+=number*cost
    for name,value in expected.items():
        quantity=ledger['resolved']['total']['quantities'][name]
        assert quantity['available']
        assert quantity['value']==pytest.approx(value,rel=3e-13,abs=0.)
        assert ledger['matched_uniform_width']['total']['quantities'][name]['exact_sum']==quantity['exact_sum']
    assert {row['key']['state_id'] for row in ledger['resolved']['by_state']}=={'M:a','M:b'}
    with (exported/'ledger.csv').open(newline='') as stream:
        rows=list(csv.DictReader(stream))
    assert rows and {row['quantity'] for row in rows}==set(expected)
    assert all(row['window_index']=='' and row['window_start_s']=='' for row in rows)
    assert all('{' not in value for row in rows for value in row.values())
    for name,digest in report['input_export_sha256'].items():
        assert sha(preserved_packets/name)==sha(exported/'packet_export'/name)==digest
    for name,digest in report['output_sha256'].items():assert sha(exported/name)==digest


def test_export_does_not_regenerate_sources_packets_moments_or_randomness(preserved_packets,tmp_path,monkeypatch):
    module=importlib.import_module('tools.export_packet_source_ledger')
    def forbidden(*args,**kwargs):raise AssertionError('unchanged generation routine called')
    monkeypatch.setattr(importlib.import_module('mhdlab.temkin_planar'),'evolve_temkin_planar_release',forbidden)
    monkeypatch.setattr(importlib.import_module('mhdlab.temkin_windowed_boundary'),'sample_temkin_windowed_emission_packets',forbidden)
    monkeypatch.setattr(importlib.import_module('tools.inspect_windowed_source_moments'),'windowed_joint_source_moments',forbidden)
    monkeypatch.setattr(np.random,'default_rng',forbidden)
    output=tmp_path/'result'
    monkeypatch.setattr(sys,'argv',command(preserved_packets,output))
    module.main()
    assert (output/'ledger.json').is_file()


@pytest.mark.parametrize('kind',['directory','file','broken_symlink'])
def test_existing_output_is_preserved(preserved_packets,tmp_path,kind):
    output=tmp_path/'existing'
    if kind=='directory':output.mkdir();(output/'keep').write_bytes(b'preserved')
    elif kind=='file':output.write_bytes(b'preserved')
    else:output.symlink_to(tmp_path/'missing')
    process=run(preserved_packets,output)
    assert process.returncode!=0
    assert output.exists() or output.is_symlink()
    if kind=='directory':assert (output/'keep').read_bytes()==b'preserved'
    elif kind=='file':assert output.read_bytes()==b'preserved'


@pytest.mark.parametrize('name',['source_input.h5','report.json'])
def test_changed_original_input_during_csv_write_rolls_back(preserved_packets,tmp_path,monkeypatch,name):
    module=importlib.import_module('tools.export_packet_source_ledger')
    copied=tmp_path/'copied';shutil.copytree(preserved_packets,copied)
    output=tmp_path/'result';original=module.write_ledger_csv
    def mutate(path,records):
        original(path,records)
        with (copied/name).open('ab') as stream:stream.write(b' ')
    monkeypatch.setattr(module,'write_ledger_csv',mutate)
    monkeypatch.setattr(sys,'argv',command(copied,output))
    with pytest.raises(SystemExit) as error:module.main()
    assert error.value.code==2 and not output.exists()
    assert not list(tmp_path.glob('.packet-source-ledger-*'))


def test_partial_write_failure_rolls_back(preserved_packets,tmp_path,monkeypatch):
    module=importlib.import_module('tools.export_packet_source_ledger')
    def fail(path,records):path.write_bytes(b'incomplete');raise OSError('manufactured write failure')
    monkeypatch.setattr(module,'write_ledger_csv',fail)
    output=tmp_path/'result';monkeypatch.setattr(sys,'argv',command(preserved_packets,output))
    with pytest.raises(SystemExit) as error:module.main()
    assert error.value.code==2 and not output.exists()
    assert not list(tmp_path.glob('.packet-source-ledger-*'))


def test_code_pin_change_during_export_rolls_back(preserved_packets,tmp_path,monkeypatch):
    module=importlib.import_module('tools.export_packet_source_ledger')
    original=module.pin_implementation;calls=[]
    def changed():
        calls.append(True);pins=original()
        if len(calls)>1:pins['tools/export_packet_source_ledger.py']='0'*64
        return pins
    monkeypatch.setattr(module,'pin_implementation',changed)
    output=tmp_path/'result';monkeypatch.setattr(sys,'argv',command(preserved_packets,output))
    with pytest.raises(SystemExit) as error:module.main()
    assert error.value.code==2 and not output.exists()


def test_indirect_parent_input_is_refused(preserved_packets,tmp_path):
    link=tmp_path/'link';link.symlink_to(preserved_packets.parent,target_is_directory=True)
    output=tmp_path/'result';process=run(link/preserved_packets.name,output)
    assert process.returncode!=0 and not output.exists()


@pytest.mark.parametrize('location',['conditioned_source','packets_source','npz_metadata'])
def test_saved_source_metadata_types_are_replayed_strictly(preserved_packets,tmp_path,location):
    copied=tmp_path/'copied';shutil.copytree(preserved_packets,copied)
    report_path=copied/'report.json';report=json.loads(report_path.read_text())
    if location=='npz_metadata':
        path=copied/'boundary_packets.npz'
        with np.load(path,allow_pickle=False) as stored:
            arrays={name:stored[name] for name in stored.files}
        metadata=json.loads(str(arrays['metadata_json']))
        metadata['source']['components'][0]['source_window']['quantities']['released_number']['qualified']=1
        arrays['metadata_json']=np.asarray(json.dumps(metadata))
        np.savez_compressed(path,**arrays)
        report['output_sha256'][path.name]=sha(path)
    else:
        source=report['conditioned_source'] if location=='conditioned_source' else report['packets']['source']
        source['components'][0]['source_window']['quantities']['released_number']['qualified']=1
    write_json(report_path,report)
    output=tmp_path/'result';process=run(copied,output)
    assert process.returncode!=0 and not output.exists()
    assert 'JSON types differ' in process.stderr
