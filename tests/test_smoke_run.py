from pathlib import Path
import json

from mhdlab import load_config, run_from_config


def test_small_end_to_end_run(tmp_path):
    cfg = load_config("configs/mykonos.yaml")
    cfg["output_root"] = str(tmp_path)
    cfg["grid"] = {"nx": 16, "ny": 12}
    cfg["time"] = {"dt_s": 1.0e-9, "end_s": 3.0e-9, "sample_every": 1}
    cfg["bolsig"]["enabled"] = False
    cfg["mhd"]["poisson_iterations"] = 10
    cfg["mhd"]["electrostatic_iterations"] = 10
    result = run_from_config(cfg)
    run_dir = Path(result["run_dir"])
    assert (run_dir / "fields.npz").exists()
    assert (run_dir / "summary.json").exists()
    assert (run_dir / "synthetic_spectrum.csv").exists()
    assert result["rl_fit"]["resistance_ohm"] > 0.0
    summary = json.loads((run_dir / "summary.json").read_text(encoding="utf-8"))
    assert summary["cross_sections"]["loaded_count"] >= 1
    assert summary["cr_model"]["reference_count"] >= 12
