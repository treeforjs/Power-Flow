import numpy as np
import pytest


def test_native_core_laplacian_if_built():
    core = pytest.importorskip("mhdlab._mhd_core")
    out = core.laplacian2d([0, 0, 0, 0, 1, 0, 0, 0, 0], 3, 3, 1.0, 1.0)
    assert np.isclose(out[4], -4.0)
