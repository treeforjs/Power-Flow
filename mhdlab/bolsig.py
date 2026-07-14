"""BOLSIG+ deck generation, execution, and lightweight parsing."""

from __future__ import annotations

import hashlib
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


@dataclass
class BolsigTable:
    reduced_field_td: np.ndarray
    mean_energy_ev: np.ndarray
    rate_coefficients: dict[str, np.ndarray] = field(default_factory=dict)
    source_file: str | None = None
    log_file: str | None = None
    return_code: int | None = None
    warning: str | None = None

    def interp_mean_energy(self, en_td: np.ndarray | float) -> np.ndarray:
        return np.interp(en_td, self.reduced_field_td, self.mean_energy_ev)


@dataclass
class _BolsigProcessResult:
    return_code: int | None
    stdout: str
    stderr: str
    warning: str | None = None


class BolsigRunner:
    def __init__(self, bolsig_dir: str | Path, cache_dir: str | Path):
        self.bolsig_dir = Path(bolsig_dir).resolve()
        self.exe = self.bolsig_dir / "bolsigminus.exe"
        if not self.exe.exists():
            self.exe = self.bolsig_dir / "bolsigplus.exe"
        self.database = self.bolsig_dir / "SigloDataBase-LXCat-04Jun2013.txt"
        self.cache_dir = Path(cache_dir).resolve()
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def available(self) -> bool:
        return self.exe.exists() and self.database.exists()

    def run_table(
        self,
        species: list[str],
        fractions: list[float],
        en_min_td: float,
        en_max_td: float,
        count: int,
        gas_temperature_k: float = 300.0,
        gas_density_m3: float = 3.295e22,
        timeout_s: float = 20.0,
    ) -> BolsigTable:
        if not self.available():
            raise FileNotFoundError(f"BOLSIG+ bundle not found at {self.bolsig_dir}")
        key = self._hash_inputs(species, fractions, en_min_td, en_max_td, count, gas_temperature_k, gas_density_m3)
        input_path = self.cache_dir / f"bolsig_{key}.dat"
        output_path = self.cache_dir / f"bolsig_{key}_out.dat"
        stdout_path = self.cache_dir / f"bolsig_{key}_stdout.txt"
        stderr_path = self.cache_dir / f"bolsig_{key}_stderr.txt"
        log_path = self.cache_dir / "bolsiglog.txt"
        return_code: int | None = None
        warning: str | None = None
        if not output_path.exists():
            cache_database = self.cache_dir / self.database.name
            if not cache_database.exists():
                cache_database.write_bytes(self.database.read_bytes())
            if log_path.exists():
                log_path.unlink()
            input_path.write_text(
                self._input_deck(species, fractions, en_min_td, en_max_td, count, gas_temperature_k, gas_density_m3, output_path.name),
                encoding="utf-8",
            )
            completed = _run_bolsig_process(
                [str(self.exe), str(input_path.name)],
                cwd=self.cache_dir,
                output_path=output_path,
                log_path=log_path,
                timeout_s=timeout_s,
            )
            return_code = completed.return_code
            warning = completed.warning
            stdout_path.write_text(completed.stdout or "", encoding="utf-8")
            stderr_path.write_text(completed.stderr or "", encoding="utf-8")
            if not output_path.exists():
                raise RuntimeError(
                    f"BOLSIG+ exited with {return_code} and did not create {output_path.name}; "
                    f"see {stdout_path.name}, {stderr_path.name}, and bolsiglog.txt"
                )
            if return_code not in (0, None):
                warning = _join_warning(
                    warning,
                    f"BOLSIG+ exited with {return_code} after writing {output_path.name}; "
                    "accepting output because it parsed successfully",
                )
        table = parse_bolsig_output(output_path)
        table.log_file = str(log_path) if log_path.exists() else None
        table.return_code = return_code
        table.warning = warning
        return table

    @staticmethod
    def _hash_inputs(*items) -> str:
        h = hashlib.sha256()
        for item in items:
            h.update(repr(item).encode("utf-8"))
        return h.hexdigest()[:16]

    def _input_deck(
        self,
        species: list[str],
        fractions: list[float],
        en_min_td: float,
        en_max_td: float,
        count: int,
        gas_temperature_k: float,
        gas_density_m3: float,
        output_name: str,
    ) -> str:
        frac_text = " ".join(f"{f:.8g}" for f in fractions)
        species_text = " ".join(species)
        db_name = self.database.name
        return f"""/NOSCREEN
/NOLOGFILE

CLEARCOLLISIONS

READCOLLISIONS
{db_name}
{species_text}
1

CONDITIONS
10.
0.
0.
{gas_temperature_k:.8g}
{gas_temperature_k:.8g}
0.
0.
{gas_density_m3:.8g}
1.
1.
3
1
1
0.
200
0
200.
1e-10
1e-4
1000
{frac_text}
1

RUNSERIES
1
{en_min_td:.8g} {en_max_td:.8g}
{count}
3

SAVERESULTS
{output_name}
4
1
1
1
0
0
1
1
0

END
"""


def parse_bolsig_output(path: str | Path) -> BolsigTable:
    """Parse common E/N formatted BOLSIG output.

    The BOLSIG text format varies with output flags. This parser extracts numeric
    rows and uses the first two columns as E/N and mean energy when headers are
    unavailable; rate blocks are preserved later as this project accumulates
    curated decks.
    """
    path = Path(path)
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    en, mean = _parse_named_two_column_block(lines, "Mean energy (eV)")
    if en.size == 0:
        raise ValueError(f"could not find an E/N vs mean-energy block in {path}")
    rate_coefficients = _parse_rate_coefficient_blocks(lines)
    order = np.argsort(en)
    return BolsigTable(
        reduced_field_td=en[order],
        mean_energy_ev=mean[order],
        rate_coefficients=rate_coefficients,
        source_file=str(path),
    )


def _parse_named_two_column_block(lines: list[str], title: str) -> tuple[np.ndarray, np.ndarray]:
    rows: list[tuple[float, float]] = []
    in_block = False
    for raw in lines:
        line = raw.strip()
        if not in_block:
            if "E/N" in line and title in line:
                in_block = True
            continue
        if not line:
            if rows:
                break
            continue
        parts = line.replace(",", " ").split()
        if len(parts) < 2:
            if rows:
                break
            continue
        try:
            rows.append((float(parts[0]), float(parts[1])))
        except ValueError:
            if rows:
                break
    if not rows:
        return np.asarray([]), np.asarray([])
    data = np.asarray(rows, dtype=float)
    valid = np.isfinite(data[:, 0]) & np.isfinite(data[:, 1]) & (data[:, 0] >= 0.0)
    return data[valid, 0], data[valid, 1]


def _parse_rate_coefficient_blocks(lines: list[str]) -> dict[str, np.ndarray]:
    rates: dict[str, np.ndarray] = {}
    previous_nonempty = ""
    for index, raw in enumerate(lines):
        line = raw.strip()
        if "E/N" not in line or "Rate coefficient (m3/s)" not in line:
            if line:
                previous_nonempty = line
            continue
        if not previous_nonempty.startswith("C"):
            continue
        label = _sanitize_rate_label(previous_nonempty)
        _, values = _parse_numeric_rows_after(lines, index + 1)
        if values.size:
            rates[f"{label}_m3_s"] = values
        previous_nonempty = line
    return rates


def _parse_numeric_rows_after(lines: list[str], start: int) -> tuple[np.ndarray, np.ndarray]:
    rows: list[tuple[float, float]] = []
    for raw in lines[start:]:
        line = raw.strip()
        if not line:
            if rows:
                break
            continue
        parts = line.replace(",", " ").split()
        if len(parts) < 2:
            if rows:
                break
            continue
        try:
            rows.append((float(parts[0]), float(parts[1])))
        except ValueError:
            if rows:
                break
    if not rows:
        return np.asarray([]), np.asarray([])
    data = np.asarray(rows, dtype=float)
    return data[:, 0], data[:, 1]


def _sanitize_rate_label(label: str) -> str:
    import re

    label = label.strip().lower()
    label = label.replace("+", "plus")
    label = re.sub(r"[^a-z0-9]+", "_", label)
    return label.strip("_")


def _run_bolsig_process(
    command: list[str],
    cwd: Path,
    output_path: Path,
    log_path: Path,
    timeout_s: float,
) -> _BolsigProcessResult:
    proc = subprocess.Popen(
        command,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    deadline = time.monotonic() + float(timeout_s)
    while True:
        return_code = proc.poll()
        if return_code is not None:
            stdout, stderr = proc.communicate()
            return _BolsigProcessResult(int(return_code), stdout or "", stderr or "")
        if output_path.exists() and _bolsig_log_finished(log_path):
            stdout, stderr = _stop_process(proc)
            return _BolsigProcessResult(
                proc.returncode,
                stdout,
                stderr,
                f"BOLSIG+ wrote {output_path.name} and reported FINISHED before process exit; "
                "stopped the wrapper process and accepted parsed output",
            )
        if time.monotonic() >= deadline:
            stdout, stderr = _stop_process(proc)
            if not output_path.exists():
                raise TimeoutError(
                    f"BOLSIG+ timed out after {timeout_s:g} s and did not create {output_path.name}; "
                    "see stdout, stderr, and bolsiglog.txt"
                )
            finished_note = " after bolsiglog.txt reported FINISHED" if _bolsig_log_finished(log_path) else ""
            return _BolsigProcessResult(
                proc.returncode,
                stdout,
                stderr,
                f"BOLSIG+ timed out after {timeout_s:g} s{finished_note}; "
                f"accepting {output_path.name} because it parsed successfully",
            )
        time.sleep(0.25)


def _stop_process(proc: subprocess.Popen) -> tuple[str, str]:
    if proc.poll() is None:
        proc.terminate()
    try:
        stdout, stderr = proc.communicate(timeout=2.0)
    except subprocess.TimeoutExpired:
        proc.kill()
        stdout, stderr = proc.communicate()
    return _timeout_payload_to_text(stdout), _timeout_payload_to_text(stderr)


def _bolsig_log_finished(log_path: Path) -> bool:
    return log_path.exists() and "FINISHED" in log_path.read_text(encoding="utf-8", errors="ignore")


def _join_warning(first: str | None, second: str) -> str:
    return f"{first}; {second}" if first else second


def _timeout_payload_to_text(payload) -> str:
    if payload is None:
        return ""
    if isinstance(payload, bytes):
        return payload.decode("utf-8", errors="replace")
    return str(payload)


def approximate_bolsig_table() -> BolsigTable:
    en = np.geomspace(0.1, 1000.0, 80)
    mean = 0.25 + 0.085 * np.sqrt(en)
    rates = {
        "H2O_dissociation_s": 1.0e5 * np.exp(-20.0 / np.maximum(mean, 1e-6)),
        "impact_ionization_s": 1.0e4 * np.exp(-13.6 / np.maximum(mean, 1e-6)),
    }
    return BolsigTable(reduced_field_td=en, mean_energy_ev=mean, rate_coefficients=rates, source_file="approximate")


def supplement_missing_rate_coefficients(table: BolsigTable, fallback: BolsigTable) -> BolsigTable:
    added = []
    for name, values in fallback.rate_coefficients.items():
        if name in table.rate_coefficients:
            continue
        table.rate_coefficients[name] = np.interp(
            table.reduced_field_td,
            fallback.reduced_field_td,
            np.asarray(values, dtype=float),
        )
        added.append(name)
    if added:
        suffix = f"supplemented missing provisional rates: {', '.join(sorted(added))}"
        table.warning = f"{table.warning}; {suffix}" if table.warning else suffix
    return table
