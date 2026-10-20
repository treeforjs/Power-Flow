import json
from pathlib import Path
import subprocess
import sys

from mhdlab.heating_comparison import load_heating_profile_csv


ROOT = Path(__file__).resolve().parents[1]


def test_reference_cli_writes_typed_profiles_energy_budget_and_point_temperature(tmp_path):
    config = json.loads((ROOT / "configs/heating_diffusion_manufactured.json").read_text())
    config["depth_cells"] = 8
    path = tmp_path / "parameters.json"
    path.write_text(json.dumps(config))
    output = tmp_path / "reference"
    result = subprocess.run([sys.executable, str(ROOT / "tools/run_heating_benchmark.py"), "--config", str(path),
                             "--output-dir", str(output), "--plot"], check=True, capture_output=True, text=True)
    assert json.loads(result.stdout)["parameter_status"] == "manufactured"
    summary = json.loads((output / "summary.json").read_text())
    assert summary["is_flash_result"] is False
    assert len(summary["profiles"]) == 3
    for item in summary["profiles"]:
        profile = load_heating_profile_csv(output / item["csv"])
        assert profile.source_kind == "synthetic"
        assert len(profile.magnetic_field_t) == 8
    assert (output / "energy_budget.csv").is_file()
    assert (output / "adiabatic_point_temperature.csv").is_file()
    assert (output / "overview.png").is_file()
