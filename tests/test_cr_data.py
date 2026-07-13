from mhdlab import load_config
from mhdlab.config import resolve_path
from mhdlab.cr_data import CRDataLibrary


def test_cr_data_library_loads_yacora_references():
    cfg = load_config("configs/mykonos.yaml")
    library = CRDataLibrary.from_config(cfg, resolve_path)
    summary = library.summary()
    assert summary["reference_count"] >= 12
    assert summary["reaction_overview_count"] >= 12
    assert "h2_electron_vibrational_excitation" in summary["bolsig_candidate_sets"]
