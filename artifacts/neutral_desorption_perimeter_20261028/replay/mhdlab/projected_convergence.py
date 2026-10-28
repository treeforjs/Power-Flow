"""Qualification of independent resolution changes in projected observables.

Successive differences are sensitivity indicators, not discretization error
bounds. A physical contrast additionally needs independently justified numerical
uncertainty bounds supplied by the forward-model owner.
"""
from __future__ import annotations

import math
from numbers import Real
from copy import deepcopy

AXES = ("mesh", "time", "spectral", "state")


def _label(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} requires a nonempty label")
    return value


def _finite(value, name, *, nonnegative=False):
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} requires a finite number")
    value = float(value)
    if not math.isfinite(value) or (nonnegative and value < 0):
        raise ValueError(f"{name} requires a finite" + (" nonnegative" if nonnegative else "") + " number")
    return value


def qualify_projected_convergence(request: dict) -> dict:
    """Inspect a shared-feature, two-hypothesis resolution qualification request.

    Each axis changes exactly one positive integer resolution count relative to
    a common finest case, with all physics, gates and feature definitions held
    fixed. Larger counts mean finer resolution. This does not interpolate data.
    """
    if request.get("schema") != "mhdlab-projected-convergence-v1":
        raise ValueError("unsupported projected convergence schema")
    status = request.get("status")
    if status not in {"manufactured", "owner_supplied"}:
        raise ValueError("status requires manufactured or owner_supplied")
    def evidence(value, name):
        value = _label(value, name)
        if status == "owner_supplied" and value.strip().lower().replace("_", " ") in {
                "unknown", "unavailable", "missing", "pending", "none", "n/a", "tbd", "not provided"}:
            raise ValueError(f"{name} cannot be placeholder evidence")
        return value

    common = evidence(request.get("comparison_contract"), "comparison_contract")
    features = request.get("features")
    if not isinstance(features, list) or not features:
        raise ValueError("features must be a nonempty list")
    ids = []
    for feature in features:
        if not isinstance(feature, dict) or set(feature) != {"id", "units"}:
            raise ValueError("features require exactly id and units; no arbitrary normalization floor")
        name = _label(feature.get("id"), "feature id")
        if name in ids:
            raise ValueError("duplicate feature id")
        ids.append(name)
        _label(feature.get("units"), "feature units")
    claims = request.get("claim_feature_ids", [])
    if not isinstance(claims, list) or len(set(claims)) != len(claims) or set(claims) - set(ids):
        raise ValueError("claim_feature_ids must be unique declared feature IDs")
    hypotheses = request.get("hypotheses")
    if not isinstance(hypotheses, dict) or len(hypotheses) != 2:
        raise ValueError("two explicitly named hypotheses are required")

    def values(raw, name, nonnegative=False):
        if not isinstance(raw, dict) or set(raw) != set(ids):
            raise ValueError(f"{name} must specify every shared feature exactly")
        return {key: _finite(raw[key], name, nonnegative=nonnegative) for key in ids}

    def resolution(raw):
        if not isinstance(raw, dict) or set(raw) != set(AXES) or any(
                type(item) is not int or item <= 0 for item in raw.values()):
            raise ValueError("resolution requires positive mesh/time/spectral/state counts")
        return raw

    reports, finest, uncertainties, run_ids = {}, {}, {}, set()

    def run_id(value):
        label = _label(value, "run_id")
        if label in run_ids:
            raise ValueError("run IDs must identify distinct refinement evidence")
        run_ids.add(label)
        return label

    for name, case in hypotheses.items():
        _label(name, "hypothesis id")
        if case.get("comparison_contract") != common:
            raise ValueError("hypotheses must share the comparison contract")
        fine = case["finest"]
        run_id(fine.get("run_id"))
        evidence(fine.get("provenance"), "finest provenance")
        finest[name] = values(fine.get("values"), "finest values")
        fine_resolution = resolution(fine.get("resolution"))
        refinements = case.get("coarser_by_axis")
        if not isinstance(refinements, dict) or set(refinements) - set(AXES):
            raise ValueError("coarser_by_axis has unsupported axes")
        axes = {}
        changes = {key: 0.0 for key in ids}
        for axis, coarse in refinements.items():
            run_id(coarse.get("run_id"))
            evidence(coarse.get("provenance"), "coarse provenance")
            counts = resolution(coarse.get("resolution"))
            if counts[axis] >= fine_resolution[axis] or any(
                    counts[key] != fine_resolution[key] for key in AXES if key != axis):
                raise ValueError("refine exactly the declared axis while holding other counts fixed")
            coarse_values = values(coarse.get("values"), "coarse values")
            changes_axis, ratios = {}, {}
            for key in ids:
                difference = abs(coarse_values[key] - finest[name][key])
                scale = abs(finest[name][key])
                ratio = difference / scale if scale else (0.0 if difference == 0 else None)
                if not math.isfinite(difference) or ratio is not None and not math.isfinite(ratio):
                    raise ValueError("comparison exceeds finite numerical domain")
                changes_axis[key] = difference
                ratios[key] = ratio
                changes[key] += difference
                if not math.isfinite(changes[key]):
                    raise ValueError("resolution change budget exceeds finite numerical domain")
            axes[axis] = {"absolute_change": changes_axis, "relative_change": ratios,
                          "below_one_percent": all(r is not None and r < .01 for r in ratios.values())}
        bound = case.get("numerical_uncertainty")
        if bound is not None:
            evidence(bound.get("provenance"), "numerical uncertainty provenance")
            evidence(bound.get("method"), "numerical uncertainty method")
            uncertainties[name] = values(bound.get("bounds"), "numerical uncertainty bounds", nonnegative=True)
        reports[name] = {"axes": axes, "missing_axes": sorted(set(AXES) - axes.keys()),
            "resolution_criteria_met": len(axes) == len(AXES) and all(v["below_one_percent"] for v in axes.values()),
            "sum_absolute_resolution_changes": changes,
            "uncertainty_bound_supplied": bound is not None}
    first, second = hypotheses
    contrasts = {}
    for key in ids:
        contrast = abs(finest[first][key] - finest[second][key])
        if not math.isfinite(contrast):
            raise ValueError("hypothesis contrast exceeds finite numerical domain")
        change_proxy = sum(reports[name]["sum_absolute_resolution_changes"][key] for name in hypotheses)
        bound = sum(uncertainties[name][key] for name in hypotheses) if len(uncertainties) == 2 else None
        if not math.isfinite(change_proxy) or bound is not None and not math.isfinite(bound):
            raise ValueError("combined change or uncertainty exceeds finite numerical domain")
        contrasts[key] = {"absolute_contrast": contrast, "sum_resolution_change_proxy": change_proxy,
            "combined_owner_uncertainty_bound": bound,
            "uncertainty_below_ten_percent_contrast": contrast > 0 and bound is not None and bound < .1 * contrast}
    resolutions_pass = all(row["resolution_criteria_met"] for row in reports.values())
    return {"schema": "mhdlab-projected-convergence-report-v1", "status": status,
            "comparison_contract": common, "hypotheses": reports, "contrasts": contrasts,
            "input_evidence": deepcopy(request),
            "zero_reference_policy": "zero-to-zero change is zero; nonzero change to a zero finest feature has undefined relative change and fails admission",
            "resolution_criteria_met": resolutions_pass,
            "contrast_criteria_met": bool(claims) and resolutions_pass and all(
                contrasts[key]["uncertainty_below_ten_percent_contrast"] for key in claims),
            "physical_discrimination_established": False,
            "interpretation": "Successive differences are resolution sensitivity, not error bounds. Owner-supplied bounds require independent review; manufactured cases establish no physical result."}
