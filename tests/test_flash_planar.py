import h5py
import numpy as np
import pytest

from mhdlab.flash_extract import FlashLayout, FlashSnapshot, PlanarFace
from mhdlab.flash_planar import as_planar_boundary_temperature_proxy, extract_planar_surface_history
from mhdlab.planar_release import evolve_planar_release
from mhdlab.thermal_release import BindingReservoir


def fixture(tmp_path):
    layout = FlashLayout(bounding_box_dataset="bbox", node_type_dataset="nodes",
                         temperature_dataset="lattice", material_fraction_dataset="foil",
                         bounding_box_axes=("x", "y", "z"), scalar_axes=("block", "z", "y", "x"),
                         physical_axes={"x": "width", "y": "normal", "z": "current"},
                         leaf_node_value=1, nonleaf_node_values=(2,), coordinate_unit="m", coordinate_scale_to_m=1,
                         temperature_unit="K", temperature_scale_to_k=1, temperature_offset_k=0,
                         min_material_fraction=1)
    face = PlanarFace((0, .002, .005), 0, 1, 1e-5, (0, .025))
    snapshots = []
    for index in range(2):
        path = tmp_path / f"raw{index}.h5"
        with h5py.File(path, "w") as f:
            f["bbox"] = [[[0, .005], [-.0002, 0], [0, .01]], [[0, .005], [-.0002, 0], [.01, .025]]]
            f["nodes"] = [1, 1]
            f["lattice"] = np.array([300, 1200]).reshape(2, 1, 1, 1)
            f["foil"] = np.ones((2, 1, 1, 1))
        snapshots.append(FlashSnapshot(path, index * 1e-8, "s", 1))
    return layout, face, snapshots


def extract(tmp_path, z=(0, .01, .025)):
    layout, face, snapshots = fixture(tmp_path)
    return extract_planar_surface_history(snapshots, layout, face, z, source_kind="synthetic",
                                          run_id="manufactured", source_reference="two-temperature fixture",
                                          surface_region="single face")


def test_export_retains_current_patches_and_provenance_before_nonlinear_release(tmp_path):
    surface, manifest = extract(tmp_path)
    np.testing.assert_array_equal(surface.temperature_k, [[[300, 300], [1200, 1200]]] * 2)
    assert surface.source_attributes["temperature_location"] == "material_sample_plane"
    assert surface.source_attributes["material_sample_depth_m"] == 1e-5
    assert manifest["current_patches"] == 2
    nested = surface.attributes["resolved_patch_extraction"]["patch_extraction_manifests"]
    assert nested[0]["sampling"]["current_interval_m"] == [0, .01]
    assert nested[1]["sampling"]["current_interval_m"] == [.01, .025]
    pools = [BindingReservoir("sites", "M", 1e8, 1, 1)]
    with pytest.raises(ValueError, match="boundary"):
        evolve_planar_release(surface, np.full((1, 2, 2), 1e18), pools, parameter_provenance="manufactured")
    proxy = as_planar_boundary_temperature_proxy(surface, reference="manufactured field has no normal gradient")
    result = evolve_planar_release(proxy, np.full((1, 2, 2), 1e18), pools, parameter_provenance="manufactured")
    assert np.all(result.released_m2[0, 0, 1] > 1e10 * result.released_m2[0, 0, 0])
    assert result.surface_attributes["material_sample_depth_m"] == 1e-5
    assert np.all(result.extensive_totals() > 0)


@pytest.mark.parametrize("z", [(0, .02), (0, .02, .01, .025), (0, True, .025), (0, np.inf, .025)])
def test_invalid_partition_rejected(tmp_path, z):
    with pytest.raises(ValueError):
        extract(tmp_path, z)
