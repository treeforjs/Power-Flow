"""Source/display contracts; optional exact-PDF verification stays external."""

from collections import Counter
import copy
import csv
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("calibration_recovery", ROOT/"tools/recover_dissertation_calibration.py")
recovery = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(recovery)


def artifact(stem):
    directory = ROOT/"data/dissertation"
    with (directory/(stem+'.csv')).open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    metadata = json.loads((directory/(stem+'.json')).read_text())
    return rows, metadata, (directory/(stem+'.csv')).read_bytes()


def test_lamp_curve_retains_original_vertices_literal_units_and_source_identity():
    rows,meta,payload=artifact('figure_4_18_calibration')
    assert len(rows)==243 and list(rows[0])==list(recovery.LAMP_FIELDS)
    assert meta['output']['csv_sha256']==recovery.sha256(payload)=='20482064aa5bf2c49cc15b230e46498da0fc66361e18d06fad15af026a77481b'
    assert meta['source']['pdf_sha256']==recovery.PDF_SHA256
    assert meta['source']['printed_page']==80 and meta['source']['pdf_page_1based']==105
    assert meta['observable']['literal_ordinate_unit']=='µW/(sr·nm·mm)'
    assert not meta['observable']['absolute_factor_resolved']
    assert not meta['observable']['SI_radiance_conversion_applied']
    assert not meta['observable']['is_raw_detector_count']
    wavelength=np.array([float(row['wavelength_nm']) for row in rows])
    assert wavelength[0]==pytest.approx(110.09395710969116,abs=1e-10)
    assert wavelength[-1]==pytest.approx(169.81177290618862,abs=1e-10)
    assert np.flatnonzero(np.diff(wavelength)<=0).tolist()==[159]
    assert not meta['recovery']['wavelength_strictly_increasing']
    assert not meta['recovery']['interpolation_ready']
    assert meta['recovery']['nonincreasing_wavelength_segments'][0]['from_vertex']==159
    assert meta['recovery']['source_vertex_order_preserved']
    assert [int(row['source_vertex_index']) for row in rows]==list(range(243))
    assert all(float(row['intensity_display_lower'])<float(row['plotted_intensity'])<float(row['intensity_display_upper']) for row in rows)
    assert meta['source_context']['caption_wavelength_range_nm']==[100,200]
    assert meta['source_context']['plotted_tick_range_nm']==[110,170]


def test_gain_series_keep_marker_and_polygon_roles_without_legend_contamination():
    rows,meta,payload=artifact('figure_4_19_gain')
    assert list(rows[0])==list(recovery.GAIN_FIELDS)
    assert recovery.sha256(payload)==meta['output']['csv_sha256']=='b099c5558ec3e2d8a0c056e606e2c342b713de6f4c80f8b15840cc0f78fb3cc8'
    assert Counter(row['series'] for row in rows)=={'exponential_fit':101,'measured_marker':36,'fit_error_boundary':203}
    measured=[row for row in rows if row['series']=='measured_marker']
    assert [int(row['source_path_index']) for row in measured]==list(range(3,39))
    assert all(row['record_kind']=='marker_center' and row['source_element_index']=='-1' for row in measured)
    # Native path/tick rounding causes small unsnapped gain-coordinate offsets.
    np.testing.assert_allclose([float(row['gain']) for row in measured],np.arange(20,56),rtol=0,atol=.001)
    assert 107 not in {int(row['source_path_index']) for row in rows}
    boundary=[row for row in rows if row['series']=='fit_error_boundary']
    assert (boundary[0]['page_x_pt'],boundary[0]['page_y_pt'])==(boundary[-1]['page_x_pt'],boundary[-1]['page_y_pt'])
    assert all(row['record_kind']=='polygon_boundary_vertex' for row in boundary)
    assert meta['source_context']['fit_error_statistical_definition'] is None
    assert meta['source_context']['per_marker_gate_assignment'] is None
    assert meta['recovery']['legend_symbols_excluded']


def test_logarithmic_gain_axis_and_display_extents_are_not_physical_error_bars():
    rows,meta,_=artifact('figure_4_19_gain')
    for row in rows:
        assert float(row['plotted_intensity'])==pytest.approx(10**float(row['log10_plotted_intensity']),rel=3e-15)
        assert float(row['intensity_display_lower'])<float(row['plotted_intensity'])<float(row['intensity_display_upper'])
    assert 'not marker-centre error' in meta['calibration']['display_envelope']
    assert not meta['recovery']['rows_are_independent_measurements']
    fit=[row for row in rows if row['series']=='exponential_fit']
    assert fit[0]['gain_region']=='below_measured_gain_range'
    assert fit[-1]['gain_region']=='above_measured_gain_range'
    assert not meta['source_context']['gain_extrapolation_is_validation']


def test_native_geometry_hash_changes_with_coordinates_but_not_renderer_sequence():
    path=dict(type='s',color=(0.,0.,1.),fill=None,width=1.,closePath=False,
              items=[('l',(1.,2.),(3.,4.))],seqno=7)
    first=recovery.drawings_digest([path])
    path['seqno']=8
    assert recovery.drawings_digest([path])==first
    path['items']=[('l',(1.,2.),(3.,4.01))]
    assert recovery.drawings_digest([path])!=first


def test_polyline_requires_original_connected_straight_segments():
    points=recovery.polyline_vertices({'items':[('l',(0.,0.),(1.,2.)),('l',(1.,2.),(3.,4.))]})
    assert points==[(0.,0.),(1.,2.),(3.,4.)]
    with pytest.raises(ValueError):
        recovery.polyline_vertices({'items':[('l',(0.,0.),(1.,2.)),('l',(2.,2.),(3.,4.))]})
    with pytest.raises(ValueError):
        recovery.polyline_vertices({'items':[('c',(0.,0.),(1.,2.),(2.,2.),(3.,4.))]})


def test_marker_center_and_symbol_extent_are_separate_from_measured_value():
    path={'items':[('c',)]*8,'rect':(20.,40.,26.,46.)}
    assert recovery.marker_center(path)==(23.,43.)
    path['rect']=(20.,40.,27.,46.)
    with pytest.raises(ValueError):
        recovery.marker_center(path)


def test_axis_transform_has_independent_linear_and_log_tick_checks():
    xf,xres=recovery.axis_fit([10.,30.,50.],[110.,130.,150.])
    yf,yres=recovery.axis_fit([200.,150.,100.],[1.,2.,3.])
    np.testing.assert_allclose(xf,[1.,100.],rtol=0,atol=2e-13)
    np.testing.assert_allclose(yf,[-.02,5.],rtol=0,atol=2e-14)
    values=list(recovery.calibrated_points([(30.,150.)],xf,yf,xres,yres,.5,True))[0]
    assert values[:3]==pytest.approx((130.,2.,100.))
    assert values[3]<130<values[4] and values[5]<100<values[6]
    for coords,vals in (([1.,1.],[0.,1.]),([1.,np.nan],[0.,1.]),([1.,2.],[0.])):
        with pytest.raises(ValueError): recovery.axis_fit(coords,vals)


def test_wrong_source_is_rejected_before_pdf_import_or_output(tmp_path):
    wrong=tmp_path/'wrong.pdf'
    wrong.write_bytes(b'not the inspected source')
    with pytest.raises(ValueError,match='PDF SHA256'):
        recovery.recover(wrong)
    assert list(tmp_path.iterdir())==[wrong]


def test_qa_rejects_changed_source_before_producing_an_overlay(tmp_path):
    pytest.importorskip('pymupdf')
    wrong=tmp_path/'changed.pdf'
    wrong.write_bytes(b'changed source')
    destination=tmp_path/'overlay.png'
    with pytest.raises(ValueError,match='QA source PDF changed'):
        recovery.make_qa(wrong,'4.18',[],destination)
    assert not destination.exists()


def test_optional_exact_pdf_reproduces_artifacts_and_rejects_changed_path_guard(monkeypatch):
    import os
    source=Path(os.environ.get('DISSERTATION_SOURCE_PDF','/workspace/research_sources/smtrevor_1.pdf'))
    if not source.is_file():
        pytest.skip('source PDF is an external cache artifact')
    pytest.importorskip('pymupdf')
    results=recovery.recover(source)
    for figure,stem,fields in (('4.18','figure_4_18_calibration',recovery.LAMP_FIELDS),('4.19','figure_4_19_gain',recovery.GAIN_FIELDS)):
        rows,metadata=results[figure]
        _,committed,payload=artifact(stem)
        assert recovery.csv_bytes(rows,fields)==payload
        assert metadata['source']==committed['source']
    changed=copy.deepcopy(recovery.NATIVE_DRAWINGS_SHA256)
    changed[105]='0'*64
    monkeypatch.setattr(recovery,'NATIVE_DRAWINGS_SHA256',changed)
    with pytest.raises(ValueError,match='native drawing geometry'):
        recovery.recover(source)
