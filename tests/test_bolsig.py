from pathlib import Path
import sys

import numpy as np

from mhdlab.bolsig import (
    _run_bolsig_process,
    approximate_bolsig_table,
    parse_bolsig_output,
    supplement_missing_rate_coefficients,
)


def test_parse_bolsig_mean_energy_block(tmp_path):
    output = tmp_path / "bolsig_out.dat"
    output.write_text(
        """BOLSIG+ version: 07/2024
Gas temperature (K)                        300.000

E/N (Td)      	Grid type
0.100000	0.00000
0.231013	3.00000

E/N (Td)	Mean energy (eV)
0.100000	0.432885E-01
0.231013	0.533381E-01
0.533670	0.813053E-01

E/N (Td)	Mobility *N
0.100000	1.0

C17   H2    Ionization    15.40 eV
E/N (Td)	Rate coefficient (m3/s)
0.100000	0.00000
0.231013	0.100000E-20
0.533670	0.200000E-20
""",
        encoding="utf-8",
    )
    table = parse_bolsig_output(output)
    assert table.reduced_field_td.tolist() == [0.1, 0.231013, 0.53367]
    assert np.isclose(table.mean_energy_ev[0], 0.0432885)
    assert table.source_file == str(output)
    assert "c17_h2_ionization_15_40_ev_m3_s" in table.rate_coefficients
    assert np.isclose(table.rate_coefficients["c17_h2_ionization_15_40_ev_m3_s"][-1], 2.0e-21)


def test_supplement_missing_rate_coefficients_adds_prototype_rates(tmp_path):
    output = tmp_path / "bolsig_out.dat"
    output.write_text(
        """E/N (Td)	Mean energy (eV)
0.100000	0.043
1.000000	0.2
""",
        encoding="utf-8",
    )
    table = parse_bolsig_output(output)
    supplement_missing_rate_coefficients(table, approximate_bolsig_table())
    assert "H2O_dissociation_s" in table.rate_coefficients
    assert "impact_ionization_s" in table.rate_coefficients
    assert "supplemented missing provisional rates" in table.warning


def test_bolsig_process_accepts_finished_output_before_process_exit(tmp_path):
    output = tmp_path / "fake_out.dat"
    log = tmp_path / "bolsiglog.txt"
    script = (
        "from pathlib import Path\n"
        "import time\n"
        "Path('fake_out.dat').write_text('E/N (Td)\\tMean energy (eV)\\n0.1\\t0.05\\n1.0\\t0.2\\n', encoding='utf-8')\n"
        "Path('bolsiglog.txt').write_text('FINISHED\\n', encoding='utf-8')\n"
        "time.sleep(10)\n"
    )
    result = _run_bolsig_process(
        [sys.executable, "-c", script],
        cwd=tmp_path,
        output_path=output,
        log_path=log,
        timeout_s=5.0,
    )
    table = parse_bolsig_output(output)
    assert table.reduced_field_td.tolist() == [0.1, 1.0]
    assert result.warning is not None
    assert "FINISHED" in result.warning
