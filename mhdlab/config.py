"""Configuration loading and path helpers."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def load_config(path: str | Path) -> dict[str, Any]:
    """Load a YAML or JSON configuration file."""
    cfg_path = Path(path).expanduser().resolve()
    text = cfg_path.read_text(encoding="utf-8")
    if cfg_path.suffix.lower() == ".json":
        data = json.loads(text)
    else:
        try:
            import yaml
        except ImportError as exc:
            raise RuntimeError("PyYAML is required for YAML configs; install requirements.txt") from exc
        data = yaml.safe_load(text)
    if not isinstance(data, dict):
        raise ValueError(f"configuration must be a mapping: {cfg_path}")
    data["_config_path"] = str(cfg_path)
    data["_config_dir"] = str(cfg_path.parent)
    return data


def resolve_path(config: dict[str, Any], value: str | Path | None) -> Path | None:
    """Resolve a path relative to the config file directory."""
    if value is None:
        return None
    path = Path(value).expanduser()
    if path.is_absolute():
        return path
    base = Path(config.get("_config_dir", "."))
    return (base / path).resolve()
