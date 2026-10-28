"""Validate a saved external receiver/population bundle without running its model."""

from __future__ import annotations

import argparse
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys


MODULE_PATH = Path(__file__).resolve().parents[1] / "mhdlab" / "raikou_handoff.py"
SPEC = importlib.util.spec_from_file_location("_current_raikou_handoff", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", type=Path, help="v1 manifest; evidence paths are relative to its directory")
    parser.add_argument("--instrument-config", type=Path, help="externally declared effective Gaussian profiles and bins")
    parser.add_argument("--allow-manufactured-control", action="store_true", help="explicitly allow synthetic projection with its blockers retained")
    parser.add_argument("--project-output", type=Path, help="save projection evidence to a new JSON file; otherwise include it in stdout")
    args = parser.parse_args()
    try:
        if (args.allow_manufactured_control or args.project_output) and not args.instrument_config:
            raise ValueError("projection options require --instrument-config")
        manifest_path = args.bundle.resolve(strict=True)
        manifest_bytes = manifest_path.read_bytes()
        document = MODULE.strict_json(manifest_bytes)
        artifacts = {}
        for name, reference in document["artifacts"].items():
            path = (manifest_path.parent / reference["relative_path"]).resolve(strict=True)
            if not path.is_relative_to(manifest_path.parent):
                raise ValueError("evidence paths must remain inside the bundle directory")
            artifacts[name] = path.read_bytes()
        report = MODULE.validate_raikou_handoff(document, artifacts)
        report["manifest_sha256"] = sha256(manifest_bytes).hexdigest()
        report["validator_source_sha256"] = sha256(MODULE_PATH.read_bytes()).hexdigest()
        report["external_model_evolution_calls"] = 0
        if args.instrument_config:
            instrument_bytes = args.instrument_config.read_bytes()
            prepared = MODULE.prepare_resolved_projection(document, artifacts, MODULE.strict_json(instrument_bytes),
                                                         allow_manufactured_control=args.allow_manufactured_control)
            projector_path = MODULE_PATH.parent / "resolved_emission_observation.py"
            spec = importlib.util.spec_from_file_location("_current_handoff_projector", projector_path)
            projector = importlib.util.module_from_spec(spec)
            sys.modules[spec.name] = projector
            spec.loader.exec_module(projector)
            kwargs = prepared.pop("projection_kwargs")
            kwargs["transitions"] = [projector.ResolvedTransition(**x) for x in kwargs["transitions"]]
            prepared["projection"] = projector.project_resolved_emission(**kwargs)
            prepared["instrument_config_sha256"] = sha256(instrument_bytes).hexdigest()
            prepared["projector_source_sha256"] = sha256(projector_path.read_bytes()).hexdigest()
            if args.project_output:
                raw_projection = (json.dumps(prepared, indent=2, allow_nan=False) + "\n").encode()
                args.project_output.parent.mkdir(parents=True, exist_ok=True)
                with args.project_output.open("xb") as stream:
                    stream.write(raw_projection)
                report["projection_output"] = {"path": str(args.project_output.resolve()),
                                               "sha256": sha256(raw_projection).hexdigest()}
            else:
                report["projection_evidence"] = prepared
    except (ValueError, KeyError, OSError, TypeError) as exc:
        parser.error(str(exc))
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
