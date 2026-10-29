"""Saved kinetic estimates keep their support, physical axes and timestamps."""
import json

import h5py
import numpy as np
from PIL import Image
import pytest

from mhdlab.neutral_temperature_render import inspect_temperature_maps,render_temperature


def fixture(tmp_path):
    path=tmp_path/'kinetic.h5'
    temperature=np.full((3,4,3,4),np.nan)
    temperature[0,0,0,0]=100.
    temperature[1,0,0,3]=250.
    temperature[2,0,2,3]=700.
    temperature[2,1,0,1]=0.
    temperature[:,2,2,0]=550.
    valid=np.isfinite(temperature)
    neff=np.zeros_like(temperature);neff[valid]=100.
    density=np.zeros_like(temperature);density[valid]=1e15
    with h5py.File(path,'w') as hdf:
        hdf.create_dataset('time_s',data=[.75e-9,30e-9,60e-9])
        hdf.create_dataset('x_edges_m',data=[-.003,-.001,0.,.001,.003])
        hdf.create_dataset('y_edges_m',data=[-.002,-.0001,.0001,.002])
        hdf.create_dataset('foil_bounds_m',data=[-.001,.001,-.0001,.0001])
        hdf.create_dataset('species',data=np.asarray(['H','H2','H2O','OH'],dtype=object),dtype=h5py.string_dtype('utf-8'))
        hdf.create_dataset('translational_temperature_k',data=temperature)
        hdf.create_dataset('effective_sample_size',data=neff)
        hdf.create_dataset('number_density_m3',data=density)
        hdf.attrs['status']='manufactured'
        hdf.attrs['definition_metadata_json']=json.dumps({'formula':'m trace(centered velocity covariance)/(3 kB)',
            'scope':'Authored verification values, not physical gas temperatures'})
    metadata=tmp_path/'metadata.json'
    metadata.write_text(json.dumps({'status':'manufactured','model_description':'Authored kinetic-temperature render fixture',
        'time_registration':'Authored unshifted times','assumptions':['No physical kinetic model is established'],
        'invariant_current_length_m':.01,'min_Neff':50.,'min_expected_cell_particles':1.,
        'normal_axis_magnification':1.}))
    return path,metadata,temperature


def test_temperature_encoding_retains_nans_sub300_and_supported_zero(tmp_path):
    path,metadata,original=fixture(tmp_path)
    before=path.read_bytes()
    report=render_temperature(path,metadata,tmp_path/'render',dpi=60)
    assert report['frame_count']==3
    assert report['physical_time_s']==[.75e-9,30e-9,60e-9]
    assert report['color_scale']['limits_k']==[0.,700.]
    assert report['source_summary']['species_temperature_ranges'][0]['valid_min_k']==100.
    assert report['source_summary']['species_temperature_ranges'][1]['valid_min_k']==0.
    assert report['source_summary']['species_temperature_ranges'][3]['valid_min_k'] is None
    assert report['source_summary']['valid_cell_counts']==[[1,0,1,0],[1,0,1,0],[1,1,1,0]]
    assert report['temperature_values_changed'] is False
    assert report['unsupported_nans_replaced'] is False
    assert report['temperature_floor_applied'] is False
    assert path.read_bytes()==before
    with h5py.File(path) as hdf:
        np.testing.assert_array_equal(hdf['translational_temperature_k'][...],original)
    with Image.open(tmp_path/'render/neutral_temperature.gif') as encoded:
        assert encoded.n_frames==3
        assert encoded.info['loop']==0
    for name in ('preview_first.png','preview_middle.png','preview_last.png'):
        assert (tmp_path/'render'/name).is_file()


def test_equal_scale_crop_and_explicit_color_cap_change_display_only(tmp_path):
    path,metadata,_=fixture(tmp_path)
    report=render_temperature(path,metadata,tmp_path/'cropped',dpi=60,
        view_bounds_m=[-.0025,.0025,-.0015,.0015],temperature_max_k=600.)
    assert report['geometry_aspect']=='equal_millimetre_units'
    assert report['display_view_bounds_m']==[-.0025,.0025,-.0015,.0015]
    assert report['view_crops_axes_only'] is True
    assert report['color_scale']['maximum_supported_temperature_k']==700.
    assert sum(map(sum,report['color_scale']['bins_above_color_range']))==1
    with pytest.raises(ValueError,match='inside'):
        render_temperature(path,metadata,tmp_path/'outside',dpi=60,view_bounds_m=[-.004,.004,-.002,.002])


@pytest.mark.parametrize('value',[-1.,np.inf,-np.inf])
def test_negative_or_infinite_temperature_rejects(tmp_path,value):
    path,metadata,_=fixture(tmp_path)
    with h5py.File(path,'r+') as hdf: hdf['translational_temperature_k'][0,0,0,0]=value
    with pytest.raises(ValueError,match='temperature'):
        inspect_temperature_maps(path,metadata)


@pytest.mark.parametrize('dataset,value',[('effective_sample_size',49.),('number_density_m3',0.)])
def test_finite_temperature_without_declared_support_rejects(tmp_path,dataset,value):
    path,metadata,_=fixture(tmp_path)
    with h5py.File(path,'r+') as hdf: hdf[dataset][0,0,0,0]=value
    with pytest.raises(ValueError,match='support'):
        inspect_temperature_maps(path,metadata)


def test_all_unsupported_is_explicitly_empty_without_zero_replacement(tmp_path):
    path,metadata,_=fixture(tmp_path)
    with h5py.File(path,'r+') as hdf: hdf['translational_temperature_k'][...]=np.nan
    report=render_temperature(path,metadata,tmp_path/'empty',dpi=60)
    assert report['color_scale']['maximum_supported_temperature_k'] is None
    assert all(count==0 for row in report['source_summary']['valid_cell_counts'] for count in row)
    assert report['unsupported_nans_replaced'] is False


def test_metadata_status_and_physical_axis_mismatch_reject(tmp_path):
    path,metadata,_=fixture(tmp_path)
    with h5py.File(path,'r+') as hdf: hdf.attrs['status']='provisional'
    with pytest.raises(ValueError,match='status'):
        inspect_temperature_maps(path,metadata)
    with h5py.File(path,'r+') as hdf: hdf.attrs['status']='manufactured'
    record=json.loads(metadata.read_text());record['normal_axis_magnification']=4.
    metadata.write_text(json.dumps(record))
    with pytest.raises(ValueError,match='equal'):
        inspect_temperature_maps(path,metadata)


def test_overwrite_rejects_without_changing_previous_result(tmp_path):
    path,metadata,_=fixture(tmp_path)
    output=tmp_path/'previous';output.mkdir();marker=output/'keep.txt';marker.write_text('keep')
    with pytest.raises(ValueError,match='new render'):
        render_temperature(path,metadata,output,dpi=60)
    assert marker.read_text()=='keep'
