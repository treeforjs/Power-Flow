"""Bind the literal dissertation arrays before admitting a physical ensemble."""
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("_pilot_bridge_test", ROOT / "tools/prepare_pilot_comparison.py")
bridge = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bridge)


def packet():
    manifest = json.loads((ROOT / "data/dissertation/pilot_comparison_manifest.json").read_bytes())
    return bridge.prepare_packet(manifest, "15750_4_20_wavelength")


def observation(value):
    curves = []
    for raw in value["curves"]:
        row = dict(raw)
        for key in ("coordinate", "center", "lower", "upper", "support_segment_ids"):
            row[key] = np.array(row[key])
        curves.append(SimpleNamespace(**row))
    return SimpleNamespace(curves=curves, status="reconstructed", mapping_status="qualified",
                           gate_start_s=-30e-9, gate_end_s=-27e-9, time_registration="invented")


def test_actual_pilot_packet_retains_samples_and_withholds_inference():
    value = packet()
    assert sum(len(curve["coordinate"]) for curve in value["curves"]) == 157
    assert value["allowed_C_abs"] == [1]
    assert not value["physical_comparison_ready"]
    assert "sourced monolayer convention for coverage" not in value["comparison_blockers"]
    assert value["coverage_only_blockers"] == ["sourced monolayer convention for coverage"]
    with pytest.raises(ValueError, match="withheld"):
        bridge.bind_literal_observation(value, observation(value))


@pytest.mark.parametrize("field,new", [("quantity", "spectral_energy_radiance"),
                                        ("signal_units", "W m^-2 sr^-1 nm^-1"),
                                        ("scale_role", "nominal_energy_radiance")])
def test_relabeling_derived_display_cannot_bypass_manifest(field, new):
    value = packet()
    obs = observation(value)
    setattr(obs.curves[0], field, new)
    with pytest.raises(ValueError, match="substituted " + field):
        bridge.bind_literal_observation(value, obs)


def test_changed_ordinate_or_omitted_support_rejects_before_readiness():
    value = packet()
    for field in ("center", "support_segment_ids"):
        obs = observation(value)
        getattr(obs.curves[0], field)[0] += 1
        with pytest.raises(ValueError, match="substituted " + field):
            bridge.bind_literal_observation(value, obs)


def test_gap_readout_sort_preserves_original_rows_and_infers_no_continuity():
    manifest = json.loads((ROOT / "data/dissertation/pilot_comparison_manifest.json").read_bytes())
    value = bridge.prepare_packet(manifest, "15750_4_20_gap")
    assert sum(len(curve["coordinate"]) for curve in value["curves"]) == 106
    for curve in value["curves"]:
        assert np.all(np.diff(curve["coordinate"]) > 0)
        assert len(set(curve["support_segment_ids"])) == len(curve["coordinate"])
