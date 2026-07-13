import numpy as np

from mhdlab.cross_sections import CrossSectionLibrary, collision_probability


def test_cross_section_manifest_loads_saved_nist_tables():
    library = CrossSectionLibrary.from_manifest("data/cross_sections/cross_section_manifest.yaml")
    assert "e_H2_ionization_ref75" in library.tables
    table = library.tables["e_H2_ionization_ref75"]
    assert table.sigma(15.43) == 0.0
    assert table.sigma(30.0) > 0.0
    assert table.source_file.endswith("e_H2_ionization_ref75_nist_beb.csv")


def test_collision_probability_is_bounded():
    p = collision_probability(
        target_density_m3=np.asarray([0.0, 1.0e20, 1.0e25]),
        sigma_m2=np.asarray([1.0e-20, 1.0e-20, 1.0e-20]),
        relative_speed_m_s=np.asarray([1.0e6, 1.0e6, 1.0e6]),
        dt_s=1.0e-9,
    )
    assert np.all(p >= 0.0)
    assert np.all(p <= 1.0)
    assert p[0] == 0.0
    assert p[-1] > p[1]
