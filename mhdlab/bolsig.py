"""BOLSIG+ deck generation, execution, and lightweight parsing."""

from __future__ import annotations

import hashlib
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


@dataclass
class BolsigTable:
    reduced_field_td: np.ndarray
    mean_energy_ev: np.ndarray
    rate_coefficients: dict[str, np.ndarray] = field(default_factory=dict)
    source_file: str | None = None

    def interp_mean_energy(self, en_td: np.ndarray | float) -> np.ndarray:
        return np.interp(en_td, self.reduced_field_td, self.mean_energy_ev)


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
        if not output_path.exists():
            cache_database = self.cache_dir / self.database.name
            if not cache_database.exists():
                cache_database.write_bytes(self.database.read_bytes())
            input_path.write_text(
                self._input_deck(species, fractions, en_min_td, en_max_td, count, gas_temperature_k, gas_density_m3, output_path.name),
                encoding="utf-8",
            )
            subprocess.run(
                [str(self.exe), str(input_path.name)],
                cwd=self.cache_dir,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=timeout_s,
            )
            if not output_path.exists():
                raise RuntimeError(f"BOLSIG+ completed but did not create {output_path.name}")
        return parse_bolsig_output(output_path)

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
    rows = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        parts = line.replace(",", " ").split()
        nums = []
        for part in parts:
            try:
                nums.append(float(part))
            except ValueError:
                pass
        if len(nums) >= 2:
            rows.append(nums)
    if not rows:
        raise ValueError(f"no numeric BOLSIG rows found in {path}")
    width = max(len(r) for r in rows)
    table = np.full((len(rows), width), np.nan)
    for i, row in enumerate(rows):
        table[i, : len(row)] = row
    en = table[:, 0]
    mean = table[:, 1]
    valid = np.isfinite(en) & np.isfinite(mean) & (en >= 0)
    order = np.argsort(en[valid])
    return BolsigTable(reduced_field_td=en[valid][order], mean_energy_ev=mean[valid][order], source_file=str(path))


def approximate_bolsig_table() -> BolsigTable:
    en = np.geomspace(0.1, 1000.0, 80)
    mean = 0.25 + 0.085 * np.sqrt(en)
    rates = {
        "H2O_dissociation_s": 1.0e5 * np.exp(-20.0 / np.maximum(mean, 1e-6)),
        "impact_ionization_s": 1.0e4 * np.exp(-13.6 / np.maximum(mean, 1e-6)),
    }
    return BolsigTable(reduced_field_td=en, mean_energy_ev=mean, rate_coefficients=rates, source_file="approximate")
