import h5py
import numpy as np
import pytest

from mhdlab.flash_surface import FlashSurfaceHistory, load_surface_history, write_surface_history
from mhdlab.hdf_integrity import local_dataset


@pytest.mark.parametrize("dependency", ["external_link", "soft_link", "virtual", "external_storage", "external_group"])
def test_canonical_temperatures_cannot_hide_unhashed_dependencies(tmp_path, dependency):
    canonical = tmp_path / "canonical.h5"
    other = tmp_path / "other.h5"
    write_surface_history(canonical, FlashSurfaceHistory([0, 1], [0, 1], [[300], [600]],
                                                        "synthetic", "test", "manufactured", "face"))
    with h5py.File(other, "w") as f:
        f["T"] = [[300.0], [600.0]]
        f["T"].attrs["units"] = "K"
        f["group/T"] = [[300.0], [600.0]]
    with h5py.File(canonical, "a") as f:
        del f["temperature_k"]
        if dependency == "external_link":
            f["temperature_k"] = h5py.ExternalLink(str(other), "T")
        elif dependency == "soft_link":
            f["local_T"] = [[300.0], [600.0]]
            f["local_T"].attrs["units"] = "K"
            f["temperature_k"] = h5py.SoftLink("/local_T")
        elif dependency == "virtual":
            layout = h5py.VirtualLayout((2, 1), dtype=float)
            layout[:] = h5py.VirtualSource(str(other), "T", shape=(2, 1))
            f.create_virtual_dataset("temperature_k", layout).attrs["units"] = "K"
        elif dependency == "external_storage":
            dataset = f.create_dataset("temperature_k", (2, 1), dtype=float,
                                       external=[(str(tmp_path / "external.bin"), 0, h5py.h5f.UNLIMITED)])
            dataset[:] = [[300], [600]]
            dataset.attrs["units"] = "K"
        else:
            f["external_group"] = h5py.ExternalLink(str(other), "group")
            with pytest.raises(ValueError):
                local_dataset(f, "external_group/T")
            return
    with pytest.raises(ValueError):
        load_surface_history(canonical)


def test_same_file_hard_link_and_nested_dataset_remain_valid(tmp_path):
    with h5py.File(tmp_path / "local.h5", "w") as f:
        f["group/value"] = np.arange(3)
        f["alias"] = f["group/value"]
        np.testing.assert_array_equal(local_dataset(f, "alias")[:], np.arange(3))
        np.testing.assert_array_equal(local_dataset(f, "group/value")[:], np.arange(3))
