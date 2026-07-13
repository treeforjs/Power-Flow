"""Collisional-radiative reference and data-manifest loaders."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class CRDataLibrary:
    reference_file: str
    reaction_file: str | None
    rovibrational_file: str | None
    references: dict[str, dict[str, Any]]
    reaction_overview: list[dict[str, Any]]
    rovibrational_sets: list[dict[str, Any]]

    @classmethod
    def from_config(cls, config: dict[str, Any], resolve_path) -> "CRDataLibrary | None":
        cr_cfg = config.get("cr_model", {})
        if not cr_cfg.get("enabled", True):
            return None
        references_path = resolve_path(config, cr_cfg.get("references"))
        if references_path is None:
            return None
        references_data = load_yaml(references_path)
        reaction_path = resolve_path(config, cr_cfg.get("reaction_overview"))
        rovib_path = resolve_path(config, cr_cfg.get("rovibrational_manifest"))
        reaction_data = load_yaml(reaction_path) if reaction_path else {}
        rovib_data = load_yaml(rovib_path) if rovib_path else {}
        return cls(
            reference_file=str(references_path),
            reaction_file=str(reaction_path) if reaction_path else None,
            rovibrational_file=str(rovib_path) if rovib_path else None,
            references=references_data.get("references", {}),
            reaction_overview=reaction_data.get("hydrogen_atom_reactions", []),
            rovibrational_sets=rovib_data.get("rovibrational_cross_sections", []),
        )

    def summary(self) -> dict[str, Any]:
        return {
            "reference_file": self.reference_file,
            "reaction_file": self.reaction_file,
            "rovibrational_file": self.rovibrational_file,
            "reference_count": len(self.references),
            "reaction_overview_count": len(self.reaction_overview),
            "rovibrational_set_count": len(self.rovibrational_sets),
            "bolsig_candidate_sets": [
                item["name"]
                for item in self.rovibrational_sets
                if item.get("can_feed_bolsig", False)
            ],
        }


def load_yaml(path: str | Path) -> dict[str, Any]:
    import yaml

    path = Path(path)
    if not path.exists():
        return {}
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    return data or {}
