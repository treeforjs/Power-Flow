"""Graphics primitives must not become additional measured mass channels."""
import importlib.util
from pathlib import Path
import pytest

SPEC=importlib.util.spec_from_file_location('rga_recovery',Path(__file__).resolve().parents[1]/'tools/recover_dissertation_rga.py')
MODULE=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(MODULE)


def test_overlapping_pdf_rectangles_are_fifty_channel_readouts_not450_observations():
    rectangles=[]
    for mass in range(1,51):
        for sub in range(9):
            center=mass+(sub-4)*.04
            rectangles.append(dict(mass_to_charge_amu=center,bar_left_amu=center-.4,
                bar_right_amu=center+.4,plotted_channel_pressure_torr=mass*100+sub))
    readouts=MODULE.channel_readouts(rectangles)
    assert len(rectangles)==450 and len(readouts)==50
    assert [r['nominal_mass_channel_amu'] for r in readouts]==list(range(1,51))
    assert [r['published_plot_envelope_pressure_torr'] for r in readouts]==[mass*100+8 for mass in range(1,51)]
    assert not any(r['neutral_species_partial_pressure_established'] for r in readouts)


def test_unrepresented_channel_remains_unknown_instead_of_becoming_zero():
    rows=MODULE.channel_readouts([dict(mass_to_charge_amu=18,bar_left_amu=17.6,bar_right_amu=18.4,plotted_channel_pressure_torr=8e-5)])
    assert rows[17]['published_plot_envelope_pressure_torr']==8e-5
    assert rows[4]['published_plot_envelope_pressure_torr'] is None
    assert rows[4]['encoded_center_nearest_pressure_torr'] is None


def test_unverified_pdf_is_refused_before_geometry_is_interpreted(tmp_path):
    p=tmp_path/'wrong.pdf';p.write_bytes(b'not the source-bound dissertation')
    with pytest.raises(ValueError,match='source-bound'):
        MODULE.recover(p)
