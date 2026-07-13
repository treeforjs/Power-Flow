"""Extract/save PRAB 26, 040401 cross-section provenance and starter tables.

This script saves two kinds of reusable assets:
1. The PRAB reaction/reference mapping for refs. 75-80.
2. Machine-readable starter electron-impact BEB tables from NIST where available.

The older ion-impact and charge-exchange references are preserved as manifest
entries with missing ``file`` values until digitized/measured tables are added.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen

import yaml
from pypdf import PdfReader


ROOT = Path(__file__).resolve().parents[1]
PDF_PATH = Path(r"C:\Users\tsmit\Documents\Papers\PhysRevAccelBeams.26.040401.pdf")
DATA_DIR = ROOT / "data" / "cross_sections"
NIST_DIR = DATA_DIR / "nist"


PRAB_REFERENCES = {
    75: "Y.-K. Kim and M. E. Rudd, Phys. Rev. A 50, 3954 (1994).",
    76: "G. H. Dunn and B. Van Zyl, Phys. Rev. 154, 40 (1967).",
    77: "M. B. Shah and H. B. Gilbody, J. Phys. B 14, 2361 (1981).",
    78: "W. Ott, E. Speth, and the W7-AS Team, Nucl. Fusion 42, 796 (2002).",
    79: "C. McGrath, M. B. Shah, P. C. E. McCartney, and J. W. McConkey, Phys. Rev. A 64, 062712 (2001).",
    80: "W. L. Fite, R. T. Brackmann, and W. R. Snow, Phys. Rev. 112, 1161 (1958).",
}


PRAB_REACTIONS = [
    {
        "name": "e_H_ionization_ref75",
        "reaction": "e + H -> 2e + H+",
        "incident": "e",
        "target": "H",
        "products": {"e": 2, "H+": 1},
        "process": "electron_impact_ionization",
        "reference_id": 75,
    },
    {
        "name": "e_H2_ionization_ref75",
        "reaction": "e + H2 -> 2e + H2+",
        "incident": "e",
        "target": "H2",
        "products": {"e": 2, "H2+": 1},
        "process": "electron_impact_ionization",
        "reference_id": 75,
        "nist_species": "H2",
    },
    {
        "name": "e_H2plus_dissociative_ionization_ref76",
        "reaction": "e + H2+ -> e + H+ + H",
        "incident": "e",
        "target": "H2+",
        "products": {"e": 1, "H+": 1, "H": 1},
        "process": "electron_impact_dissociative_ionization",
        "reference_id": 76,
    },
    {
        "name": "Hplus_H_ionization_ref77",
        "reaction": "H+ + H -> e + 2H+",
        "incident": "H+",
        "target": "H",
        "products": {"e": 1, "H+": 2},
        "process": "ion_impact_ionization",
        "reference_id": 77,
    },
    {
        "name": "H_H_ionization_ref78",
        "reaction": "H + H -> e + H+ + H",
        "incident": "H",
        "target": "H",
        "products": {"e": 1, "H+": 1, "H": 1},
        "process": "neutral_impact_ionization",
        "reference_id": 78,
    },
    {
        "name": "H_H2_target_ionization_ref78",
        "reaction": "H + H2 -> e + H+ + H2",
        "incident": "H",
        "target": "H2",
        "products": {"e": 1, "H+": 1, "H2": 1},
        "process": "neutral_impact_ionization",
        "reference_id": 78,
    },
    {
        "name": "H_H2_projectile_ionization_ref78",
        "reaction": "H + H2 -> e + H + H2+",
        "incident": "H",
        "target": "H2",
        "products": {"e": 1, "H": 1, "H2+": 1},
        "process": "neutral_impact_ionization",
        "reference_id": 78,
    },
    {
        "name": "Hplus_H2_ionization_ref78",
        "reaction": "H+ + H2 -> e + H+ + H2+",
        "incident": "H+",
        "target": "H2",
        "products": {"e": 1, "H+": 1, "H2+": 1},
        "process": "ion_impact_ionization",
        "reference_id": 78,
    },
    {
        "name": "H2plus_H_ionization_ref79",
        "reaction": "H2+ + H -> e + H+ + H2+",
        "incident": "H2+",
        "target": "H",
        "products": {"e": 1, "H+": 1, "H2+": 1},
        "process": "ion_impact_ionization",
        "reference_id": 79,
    },
    {
        "name": "Hplus_H_charge_exchange_ref80",
        "reaction": "H+ + H -> H + H+",
        "incident": "H+",
        "target": "H",
        "products": {"H": 1, "H+": 1},
        "process": "charge_exchange",
        "reference_id": 80,
    },
    {
        "name": "Hplus_H2_charge_exchange_ref80",
        "reaction": "H+ + H2 -> H + H2+",
        "incident": "H+",
        "target": "H2",
        "products": {"H": 1, "H2+": 1},
        "process": "charge_exchange",
        "reference_id": 80,
    },
    {
        "name": "H2plus_H_charge_exchange_ref80",
        "reaction": "H2+ + H -> H+ + H2",
        "incident": "H2+",
        "target": "H",
        "products": {"H+": 1, "H2": 1},
        "process": "charge_exchange",
        "reference_id": 80,
    },
]


NIST_TABLES = {
    "H2": "e_H2_ionization_ref75_nist_beb.csv",
    "H2O": "e_H2O_total_ionization_nist_beb.csv",
    "H2+": "e_H2plus_total_ionization_nist_beq.csv",
}


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    NIST_DIR.mkdir(parents=True, exist_ok=True)
    paper_context = extract_paper_context()
    nist_files = download_nist_tables()
    write_yaml(DATA_DIR / "prab_26_040401_reactions.yaml", {
        "source_pdf": str(PDF_PATH),
        "paper_context": paper_context,
        "references": PRAB_REFERENCES,
        "reactions": PRAB_REACTIONS,
    })
    write_yaml(DATA_DIR / "cross_section_manifest.yaml", build_manifest(nist_files))
    write_readme()
    print(f"Saved cross-section assets in {DATA_DIR}")


def extract_paper_context() -> str:
    if not PDF_PATH.exists():
        return "PRAB PDF not found; using hard-coded reaction/reference mapping."
    text = "\n".join(page.extract_text() or "" for page in PdfReader(str(PDF_PATH)).pages)
    start = text.find("Neutral ionization is mostly the result")
    if start < 0:
        return "Could not find neutral-ionization paragraph in PDF text."
    end = text.find("Dissociative recombination", start)
    if end < 0:
        end = start + 2200
    return " ".join(text[start:end].split())


def download_nist_tables() -> dict[str, str]:
    saved = {}
    for species, filename in NIST_TABLES.items():
        url = "https://physics.nist.gov/cgi-bin/Ionization/bebcsdwnload_ascii?" + quote(species, safe="")
        req = Request(url, headers={"User-Agent": "Mozilla/5.0"})
        text = urlopen(req, timeout=20).read().decode("utf-8", "replace")
        rows = parse_nist_ascii(text)
        out = NIST_DIR / filename
        with out.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(["energy_ev", "cross_section_m2"])
            writer.writerows(rows)
        saved[species] = f"nist/{filename}"
    return saved


def parse_nist_ascii(text: str) -> list[tuple[float, float]]:
    rows = []
    for line in text.splitlines():
        nums = re.findall(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][-+]?\d+)?", line)
        if len(nums) >= 2:
            energy_ev = float(nums[0])
            sigma_angstrom2 = float(nums[1])
            rows.append((energy_ev, sigma_angstrom2 * 1.0e-20))
    if not rows:
        raise ValueError("No numeric NIST cross-section rows parsed")
    return rows


def build_manifest(nist_files: dict[str, str]) -> dict:
    entries = []
    for reaction in PRAB_REACTIONS:
        ref_id = reaction["reference_id"]
        file_name = None
        source_note = "numeric table required from cited measured reference"
        if reaction.get("nist_species") in nist_files:
            file_name = nist_files[reaction["nist_species"]]
            source_note = "NIST BEB starter table; PRAB cites Kim and Rudd for this channel"
        entries.append({
            "name": reaction["name"],
            "reaction": reaction["reaction"],
            "incident": reaction["incident"],
            "target": reaction["target"],
            "products": reaction["products"],
            "process": reaction["process"],
            "reference": PRAB_REFERENCES[ref_id],
            "reference_id": ref_id,
            "file": file_name,
            "source_note": source_note,
            "units": {"energy": "eV", "cross_section": "m2"},
        })
    for species, file_name in nist_files.items():
        if species == "H2":
            continue
        entries.append({
            "name": f"nist_{species.replace('+', 'plus')}_total_ionization",
            "reaction": f"e + {species} total ionization",
            "incident": "e",
            "target": species,
            "products": {},
            "process": "electron_impact_total_ionization",
            "reference": "NIST Electron-Impact Cross Sections for Ionization and Excitation Database.",
            "reference_id": "NIST",
            "file": file_name,
            "source_note": "supporting BEB/BEQ table, not one of the PRAB refs. 75-80 plotted ion-impact channels",
            "units": {"energy": "eV", "cross_section": "m2"},
        })
    return {
        "schema": "mhdlab.cross_sections.v1",
        "source_paper": "N. Bennett et al., Phys. Rev. Accel. Beams 26, 040401 (2023).",
        "cross_sections": entries,
    }


def write_yaml(path: Path, data: dict) -> None:
    path.write_text(yaml.safe_dump(data, sort_keys=False, allow_unicode=False), encoding="utf-8")


def write_readme() -> None:
    (DATA_DIR / "README.md").write_text(
        """# Cross-Section Assets

This directory stores reusable collision cross-section provenance and tables.

- `prab_26_040401_reactions.yaml` is extracted from `PhysRevAccelBeams.26.040401.pdf` and records the refs. 75-80 reaction mapping.
- `cross_section_manifest.yaml` is the loader manifest used by `mhdlab.cross_sections.CrossSectionLibrary`.
- `nist/*.csv` contains NIST BEB/BEQ starter tables converted from Angstrom^2 to m^2.

Older measured ion-impact and charge-exchange tables from PRAB refs. 77-80 are listed in the manifest with `file: null` until digitized/tabulated data are added.
""",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
