"""New manufactured polynomial/rule controls; no source kinetics are evaluated."""

from dataclasses import FrozenInstanceError
from fractions import Fraction
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys

import pytest

_CORE_PATH = Path(__file__).resolve().parents[1] / "mhdlab/boundary_quadrature.py"
_CORE_NAME = "_boundary_quadrature_test_module_0810"
_SPEC = importlib.util.spec_from_file_location(_CORE_NAME, _CORE_PATH)
assert _SPEC is not None and _SPEC.loader is not None
_CORE = importlib.util.module_from_spec(_SPEC)
sys.modules[_CORE_NAME] = _CORE
_SPEC.loader.exec_module(_CORE)
BoundaryPiece = _CORE.BoundaryPiece
BoundaryQuery = _CORE.BoundaryQuery
AdmittedBoundarySource = _CORE.AdmittedBoundarySource
admit_boundary_source = _CORE.admit_boundary_source
build_boundary_quadrature = _CORE.build_boundary_quadrature
combine_boundary_summaries = _CORE.combine_boundary_summaries
gauss_legendre5_rule = _CORE.gauss_legendre5_rule
make_piece_quadrature = _CORE.make_piece_quadrature
restrict_boundary_piece = _CORE.restrict_boundary_piece
summarize_piece_quadrature = _CORE.summarize_piece_quadrature


F = Fraction
KB = F("1.380649e-23")
QE = F("1.602176634e-19")
PROVENANCE = "Manufactured polynomial algebra only; not a source trajectory or molecular model."
QUANTITIES = (
    "released_number", "local_time_first_moment_s", "local_time_second_moment_s2",
    "temperature_first_moment_k", "temperature_second_moment_k2",
    "local_time_temperature_moment_k_s", "binding_energy_j", "internal_energy_j",
    "translational_energy_j", "gas_energy_j", "total_wall_energy_debit_j",
)
PHASES = ("target", "finite_rule", "represented", "rule_generation_defect",
          "export_projection_defect", "total_defect")


@pytest.fixture(scope="module", autouse=True)
def _record_import_boundary(record_testsuite_property):
    forbidden = sorted(name for name in sys.modules if name == "mhdlab" or name.startswith("mhdlab."))
    assert forbidden == []
    record_testsuite_property("direct_core_module_name", _CORE_NAME)
    record_testsuite_property("legacy_modules_after_direct_core_import", json.dumps(forbidden))


def _ratio(value):
    return F(int(value["numerator"]), int(value["denominator"]))


def _exact(record):
    return _ratio(record["exact"])


def _assert_record(record, expected):
    expected = F(expected)
    assert _exact(record) == expected
    assert record["value"] == float(expected)
    projection = record["projection"]
    assert projection["exact_zero"] is (expected == 0)
    assert projection["positive"] is (expected > 0)
    assert projection["negative"] is (expected < 0)
    assert projection["finite"] is True
    assert projection["underflowed"] is False
    assert projection["overflowed"] is False
    assert _ratio(projection["exact_rounding_error"]) == F(record["value"]) - expected


def _assert_stat(summary, phase, name, expected):
    stat = summary["statistics"][phase][name]
    assert stat["available"] is True and stat["reason"] is None
    _assert_record(stat["record"], expected)


def _piece(**changes):
    # u=(t-5)/2; M=2+u+2u^2+u^3; T=10+u+u^2+u^3.
    values = dict(
        start_s=F(5), stop_s=F(7),
        emitted_controls_m2=(2, F(7, 3), F(10, 3), 6),
        temperature_controls_k=(10, F(31, 3), 11, 13),
        temperature_min_k=8, temperature_max_k=16,
        initial_coverage_m2=10, state_id="synthetic_state",
        source_history_id="synthetic_polynomial", source_cell_index=0,
        parameter_status="manufactured", parameter_provenance=PROVENANCE,
    )
    values.update(changes)
    return BoundaryPiece(**values)


@pytest.fixture(scope="module")
def rule():
    return gauss_legendre5_rule()


@pytest.fixture(scope="module")
def quadrature(rule):
    return make_piece_quadrature(_piece(), desorption_energy_ev=F(1, 2),
                                 internal_energy_ev=F(1, 4), rule=rule)


@pytest.mark.parametrize("degree", range(10))
def test_positive_rule_integrates_independent_monomial_targets_through_degree_nine(rule, degree):
    # These exact integration targets do not use the production radical formula.
    assert len(rule.nodes) == len(rule.weights) == 5
    assert all(type(value) is F for value in (*rule.nodes, *rule.weights))
    assert all(0 < x < 1 and w > 0 for x, w in zip(rule.nodes, rule.weights))
    assert list(rule.nodes) == sorted(set(rule.nodes))
    finite = sum((w * x**degree for x, w in zip(rule.nodes, rule.weights)), F())
    assert abs(finite - F(1, degree + 1)) < F(1, 10**70)


def test_nodes_keep_the_full_time_temperature_number_and_energy_relation(quadrature, rule):
    assert len(quadrature.nodes) == 5
    for node, u, w in zip(quadrature.nodes, rule.nodes, rule.weights):
        expected_time = 5 + 2*u
        expected_temperature = 10 + u + u*u + u*u*u
        expected_weight = w * (1 + 4*u + 3*u*u)
        assert node.unprojected_number_weight_m2 == expected_weight
        assert node.unprojected_time_s == expected_time
        assert node.unprojected_temperature_k == expected_temperature
        assert node.unprojected_expected_translation_j == 2*KB*expected_temperature
        assert node.number_weight_m2 == float(expected_weight) > 0
        assert node.time_s == float(expected_time)
        assert node.temperature_k == float(expected_temperature)
        assert node.expected_translational_energy_j == float(2*KB*expected_temperature)
        assert node.desorption_energy_j == float(QE/2)
        assert node.internal_energy_j == float(QE/4)
    assert len({node.temperature_k for node in quadrature.nodes}) == 5
    assert len({node.time_s for node in quadrature.nodes}) == 5


def test_canonical_energy_and_actual_binary_temperature_rounding_are_distinct(quadrature):
    differences = []
    for node in quadrature.nodes:
        actual_t = F(node.temperature_k)
        canonical_energy = 2*KB*node.unprojected_temperature_k
        actual_energy = F(node.expected_translational_energy_j)
        relation_defect = actual_energy - 2*KB*actual_t
        assert relation_defect == (
            actual_energy - canonical_energy
            + 2*KB*(node.unprojected_temperature_k - actual_t)
        )
        differences.append(relation_defect)
    assert any(difference != 0 for difference in differences)


def test_degree_eight_polynomial_targets_and_joint_statistics_are_independent(quadrature):
    summary = summarize_piece_quadrature(quadrature, local_time_anchor_s=5, area_m2=1)
    expected = {
        "released_number": F(4), "local_time_first_moment_s": F(31, 6),
        "local_time_second_moment_s2": F(116, 15),
        "temperature_first_moment_k": F(691, 15),
        "temperature_second_moment_k2": F(224051, 420),
        "local_time_temperature_moment_k_s": F(12857, 210),
        "binding_energy_j": 2*QE, "internal_energy_j": QE,
        "translational_energy_j": F(1382, 15)*KB,
        "gas_energy_j": F(1382, 15)*KB + QE,
        "total_wall_energy_debit_j": F(1382, 15)*KB + 3*QE,
    }
    assert set(summary["target"]) == set(QUANTITIES)
    for name, value in expected.items():
        _assert_record(summary["target"][name], value)
        assert abs(_exact(summary["rule_generation_defect"][name])) < abs(value)*F(1, 10**70)
        assert abs(_exact(summary["total_defect"][name])) < abs(value)*F(1, 10**13)
    for name, value in {
        "mean_local_time_s": F(31, 24), "variance_local_time_s2": F(763, 2880),
        "mean_temperature_k": F(691, 60), "variance_temperature_k2": F(9199, 12600),
        "mean_local_time_temperature_k_s": F(12857, 840),
        "covariance_local_time_temperature_k_s": F(4337, 10080),
    }.items():
        _assert_stat(summary, "target", name, value)
    assert summary["statistics"]["represented"]["variance_temperature_k2"]["record"]["value"] > 0
    assert summary["statistics"]["represented"]["covariance_local_time_temperature_k_s"]["record"]["value"] > 0


def test_actual_binary_joint_sums_and_two_defect_layers_close_exactly(quadrature):
    area, anchor = F(3, 7), F(5)
    summary = summarize_piece_quadrature(quadrature, local_time_anchor_s=anchor, area_m2=area)
    expected = {name: F() for name in QUANTITIES}
    assert summary["node_count"] == len(summary["nodes"]) == 5
    for node in summary["nodes"]:
        assert all(type(node[name]) is float and math.isfinite(node[name]) for name in
                   ("number_weight", "time_s", "temperature_k", "expected_translational_energy_j",
                    "desorption_energy_j", "internal_energy_j"))
        w, tau, temperature = F(node["number_weight"]), F(node["time_s"]) - anchor, F(node["temperature_k"])
        translation, binding, internal = (F(node[key]) for key in
                                          ("expected_translational_energy_j", "desorption_energy_j", "internal_energy_j"))
        increments = (w, w*tau, w*tau*tau, w*temperature, w*temperature*temperature,
                      w*tau*temperature, w*binding, w*internal, w*translation,
                      w*(translation+internal), w*(translation+internal+binding))
        for name, value in zip(QUANTITIES, increments):
            expected[name] += value
    for name, value in expected.items():
        _assert_record(summary["represented"][name], value)
        target = _exact(summary["target"][name])
        finite = _exact(summary["finite_rule"][name])
        generation, projection, total = finite-target, value-finite, value-target
        _assert_record(summary["rule_generation_defect"][name], generation)
        _assert_record(summary["export_projection_defect"][name], projection)
        _assert_record(summary["total_defect"][name], total)
        assert total == generation + projection
    consistency = summary["energy_consistency"]
    temperature_energy = 2*KB*expected["temperature_first_moment_k"]
    _assert_record(consistency["represented_translation_from_temperature_j"], temperature_energy)
    _assert_record(consistency["represented_translation_field_minus_temperature_energy_j"],
                   expected["translational_energy_j"] - temperature_energy)
    _assert_record(consistency["represented_binding_field_minus_supplied_cost_j"],
                   expected["binding_energy_j"] - (QE/2)*expected["released_number"])
    _assert_record(consistency["represented_internal_field_minus_supplied_cost_j"],
                   expected["internal_energy_j"] - (QE/4)*expected["released_number"])
    assert _exact(consistency["represented_translation_field_minus_temperature_energy_j"]) != 0


def test_area_views_reuse_the_same_nodes_and_apply_exact_area_once(quadrature):
    original_nodes = quadrature.nodes
    unit = summarize_piece_quadrature(quadrature, local_time_anchor_s=5, area_m2=1)
    views = [summarize_piece_quadrature(quadrature, local_time_anchor_s=5, area_m2=area)
             for area in (F(2), F(1, 7))]
    for area, view in zip((F(2), F(1, 7)), views):
        assert _ratio(view["area_exact_m2"]) == area
        for source_node, exported in zip(original_nodes, view["nodes"]):
            assert exported["number_weight"] == float(area*source_node.unprojected_number_weight_m2)
            assert exported["node_id"] == source_node.node_id
            assert exported["time_s"] == source_node.time_s
            assert exported["temperature_k"] == source_node.temperature_k
        for phase in ("target", "finite_rule", "rule_generation_defect"):
            for name in QUANTITIES:
                _assert_record(view[phase][name], area*_exact(unit[phase][name]))
    center = original_nodes[2]
    direct = float(F(1, 7)*center.unprojected_number_weight_m2)
    after_areal_projection = float(F(1, 7)*F(center.number_weight_m2))
    assert direct.hex() == "0x1.3813813813814p-3"
    assert after_areal_projection.hex() == "0x1.3813813813813p-3"
    assert views[1]["nodes"][2]["number_weight"] == direct != after_areal_projection
    assert quadrature.nodes is original_nodes
    combined = combine_boundary_summaries(views, local_time_anchor_s=5)
    for phase in PHASES:
        for name in QUANTITIES:
            _assert_record(combined[phase][name], sum((_exact(view[phase][name]) for view in views), F()))


def test_mixture_reanchors_joint_moments_before_forming_covariance(quadrature, rule):
    second_piece = _piece(start_s=7, stop_s=9, temperature_controls_k=(14, F(43, 3), 15, 17),
                          temperature_max_k=20, source_history_id="synthetic_independent_piece", source_cell_index=1)
    second_q = make_piece_quadrature(second_piece, desorption_energy_ev=F(1, 2), internal_energy_ev=F(1, 4), rule=rule)
    first = summarize_piece_quadrature(quadrature, local_time_anchor_s=5, area_m2=1)
    second = summarize_piece_quadrature(second_q, local_time_anchor_s=7, area_m2=1)
    combined = combine_boundary_summaries((first, second), local_time_anchor_s=5)
    for phase in ("target", "finite_rule", "represented"):
        def q(document, name):
            return _exact(document[phase][name])
        _assert_record(combined[phase]["local_time_first_moment_s"],
                       q(first, "local_time_first_moment_s") + q(second, "local_time_first_moment_s") + 2*q(second, "released_number"))
        _assert_record(combined[phase]["local_time_second_moment_s2"],
                       q(first, "local_time_second_moment_s2") + q(second, "local_time_second_moment_s2")
                       + 4*q(second, "local_time_first_moment_s") + 4*q(second, "released_number"))
        _assert_record(combined[phase]["local_time_temperature_moment_k_s"],
                       q(first, "local_time_temperature_moment_k_s") + q(second, "local_time_temperature_moment_k_s")
                       + 2*q(second, "temperature_first_moment_k"))
    _assert_stat(combined, "target", "variance_local_time_s2", F(3643, 2880))
    _assert_stat(combined, "target", "variance_temperature_k2", F(59599, 12600))
    _assert_stat(combined, "target", "covariance_local_time_temperature_k_s", F(24497, 10080))


def test_restriction_changes_the_measure_once_and_preserves_source_identity(rule):
    restricted = restrict_boundary_piece(_piece(), start_s=F(11, 2), stop_s=F(13, 2))
    assert restricted.source_history_id == "synthetic_polynomial"
    assert restricted.source_cell_index == 0
    assert restricted.state_id == "synthetic_state"
    assert restricted.start_s == F(11, 2) and restricted.stop_s == F(13, 2)
    # M(3/4)-M(1/4) = 61/32; no further duration factor belongs in dM/dv.
    assert restricted.emitted_controls_m2[-1] - restricted.emitted_controls_m2[0] == F(61, 32)
    q = make_piece_quadrature(restricted, desorption_energy_ev=0, internal_energy_ev=0, rule=rule)
    assert len(q.nodes) == 5
    assert abs(sum((node.unprojected_number_weight_m2 for node in q.nodes), F()) - F(61, 32)) < F(1, 10**68)
    assert all(F(11, 2) < node.unprojected_time_s < F(13, 2) for node in q.nodes)


def test_exact_zero_piece_has_no_fake_record_and_retains_its_reason(rule):
    piece = _piece(emitted_controls_m2=(2, 2, 2, 2))
    q = make_piece_quadrature(piece, desorption_energy_ev=F(1, 2), internal_energy_ev=F(1, 4), rule=rule)
    assert q.nodes == ()
    assert q.skip_reason == "exact_zero_emitted_mass"
    assert q.piece.source_cell_index == 0
    assert q.piece.source_history_id == "synthetic_polynomial"


def test_zero_piece_summary_keeps_a_reasoned_skip_catalog_and_no_means(rule):
    q = make_piece_quadrature(_piece(emitted_controls_m2=(2, 2, 2, 2)),
                               desorption_energy_ev=0, internal_energy_ev=0, rule=rule)
    summary = summarize_piece_quadrature(q, local_time_anchor_s=5, area_m2=F(3, 7))
    assert summary["node_count"] == 0 and summary["nodes"] == []
    assert len(summary["skipped_pieces"]) == 1
    assert summary["skipped_pieces"][0]["reason"] == "exact_zero_emitted_mass"
    for phase in PHASES:
        for name in QUANTITIES:
            _assert_record(summary[phase][name], 0)
    for stats in summary["statistics"].values():
        for record in stats.values():
            assert record == {"available": False, "reason": "exact_zero_emitted_mass", "record": None}


def test_piece_and_nodes_are_immutable(quadrature):
    with pytest.raises((FrozenInstanceError, AttributeError)):
        quadrature.piece.state_id = "different"
    with pytest.raises((FrozenInstanceError, AttributeError)):
        quadrature.nodes[0].temperature_k = 1.0
    summary = summarize_piece_quadrature(quadrature, local_time_anchor_s=5, area_m2=1)
    json.dumps(summary, allow_nan=False)
    summary["target"]["released_number"]["exact"]["numerator"] = "999"
    fresh = summarize_piece_quadrature(quadrature, local_time_anchor_s=5, area_m2=1)
    _assert_record(fresh["target"]["released_number"], 4)


@pytest.mark.parametrize("changes", [
    {"emitted_controls_m2": (2, 4, 3, 6)},
    {"emitted_controls_m2": (-1, 2, 3, 6)},
    {"emitted_controls_m2": (2, 3, 5, 11)},
    {"temperature_controls_k": (10, 17, 11, 13)},
    {"temperature_controls_k": (10, True, 11, 13)},
    {"emitted_controls_m2": (2, 3, float("nan"), 6)},
    {"source_cell_index": True},
    {"source_history_id": " "},
    {"state_id": ""},
    {"parameter_status": "sourced"},
])
def test_piece_support_control_and_source_identity_guards(changes):
    with pytest.raises(ValueError):
        _piece(**changes)


@pytest.mark.parametrize("start,stop", [(4, 6), (6, 8), (6, 6), (7, 5)])
def test_restriction_rejects_extrapolation_or_empty_support(start, stop):
    with pytest.raises(ValueError):
        restrict_boundary_piece(_piece(), start_s=start, stop_s=stop)


def test_existing_piece_is_revalidated_after_control_forgery(rule):
    forged = _piece()
    object.__setattr__(forged, "temperature_controls_k", (F(10), F(100), F(11), F(13)))
    with pytest.raises(ValueError):
        make_piece_quadrature(forged, desorption_energy_ev=0, internal_energy_ev=0, rule=rule)


@pytest.mark.parametrize("kind", ["positive_weight_underflow", "positive_cost_underflow", "global_time_projection"])
def test_unsupported_positive_or_geometric_projections_reject_without_repair(rule, kind):
    with pytest.raises(ValueError):
        if kind == "positive_weight_underflow":
            tiny = F(1, 10**400)
            piece = _piece(emitted_controls_m2=(0, tiny, 2*tiny, 3*tiny))
            make_piece_quadrature(piece, desorption_energy_ev=0, internal_energy_ev=0, rule=rule)
        elif kind == "positive_cost_underflow":
            make_piece_quadrature(_piece(), desorption_energy_ev=F(1, 10**400), internal_energy_ev=0, rule=rule)
        else:
            origin = F(2**60)
            piece = _piece(start_s=origin, stop_s=origin + 2)
            make_piece_quadrature(piece, desorption_energy_ev=0, internal_energy_ev=0, rule=rule)


@pytest.mark.parametrize("changes", [
    {"area_m2": 0}, {"area_m2": F(-1, 10**1000)},
    {"area_m2": True}, {"local_time_anchor_s": float("inf")},
])
def test_summary_area_and_anchor_guards_do_not_normalize_or_clip(quadrature, changes):
    values = dict(local_time_anchor_s=5, area_m2=1)
    values.update(changes)
    with pytest.raises(ValueError):
        summarize_piece_quadrature(quadrature, **values)


def _invalid_admission_bundle():
    files = {alias: b"{}" for alias in ("report", "windows", "input_config", "window_config",
                                        "independent_completion", "independent_report", "scientific_review")}
    hashes = {alias: hashlib.sha256(raw).hexdigest() for alias, raw in files.items()}
    return files, hashes


@pytest.mark.parametrize("variant", ["missing_alias", "extra_alias", "wrong_hash", "duplicate_json", "nonfinite_json", "wrong_schema"])
def test_new_source_admission_rejects_invalid_bytes_without_old_module_calls(variant):
    files, hashes = _invalid_admission_bundle()
    if variant == "missing_alias":
        files.pop("windows")
    elif variant == "extra_alias":
        files["unexpected"] = b"{}"
    elif variant == "wrong_hash":
        files["report"] = b"not-json"
    else:
        files["report"] = {
            "duplicate_json": b'{"schema":"one","schema":"two"}',
            "nonfinite_json": b'{"schema":"one","probe":NaN}',
            "wrong_schema": b'{"schema":"manufactured-invalid-source"}',
        }[variant]
        hashes["report"] = hashlib.sha256(files["report"]).hexdigest()
    with pytest.raises(ValueError) as exc:
        admit_boundary_source(files, expected_sha256=hashes, parameter_provenance=PROVENANCE)
    if variant == "wrong_hash":
        assert "sha" in str(exc.value).lower() or "digest" in str(exc.value).lower()


@pytest.mark.parametrize("token,message", [
    (b"1e-400", "nonzero JSON float underflows its representation"),
    (b"-0", "negative-zero JSON integer is not canonical"),
])
def test_new_source_parser_rejects_token_projection_before_schema_admission(token, message):
    files, hashes = _invalid_admission_bundle()
    files["report"] = b'{"schema":"invalid-test-schema","probe":' + token + b'}'
    hashes["report"] = hashlib.sha256(files["report"]).hexdigest()
    with pytest.raises(ValueError, match=message):
        admit_boundary_source(files, expected_sha256=hashes, parameter_provenance=PROVENANCE)


def test_rule_ignores_hostile_caller_context_before_its_first_factory_call(rule, record_testsuite_property):
    # A fresh process prevents the production rule cache from masking ambient traps.
    script = r'''
import decimal
import importlib.util
import json
import sys

name = "_boundary_quadrature_hostile_test_0810"
spec = importlib.util.spec_from_file_location(name, sys.argv[1])
assert spec is not None and spec.loader is not None
module = importlib.util.module_from_spec(spec)
sys.modules[name] = module
spec.loader.exec_module(module)
forbidden = sorted(n for n in sys.modules if n == "mhdlab" or n.startswith("mhdlab."))
assert forbidden == []
def signature(context):
    return (context.prec, context.rounding, context.Emin, context.Emax,
            context.capitals, context.clamp, dict(context.traps), dict(context.flags))
before = signature(decimal.getcontext())
with decimal.localcontext() as caller:
    caller.prec = 2
    caller.rounding = decimal.ROUND_UP
    caller.Emin, caller.Emax = -2, 2
    caller.capitals, caller.clamp = 0, 1
    caller.clear_flags()
    for signal in caller.traps:
        caller.traps[signal] = True
    hostile_before = signature(caller)
    generated = module.gauss_legendre5_rule()
    assert signature(caller) == hostile_before
assert signature(decimal.getcontext()) == before
print(json.dumps({"nodes": [str(x) for x in generated.nodes],
                  "weights": [str(x) for x in generated.weights],
                  "forbidden_legacy_modules": forbidden,
                  "caller_context_restored": True}, sort_keys=True))
'''
    process = subprocess.run([sys.executable, "-W", "error", "-c", script, str(_CORE_PATH)],
                             cwd=_CORE_PATH.parents[1], text=True, capture_output=True, timeout=30)
    assert process.returncode == 0, process.stderr
    actual = json.loads(process.stdout)
    assert actual["nodes"] == [str(x) for x in rule.nodes]
    assert actual["weights"] == [str(x) for x in rule.weights]
    assert actual["forbidden_legacy_modules"] == []
    assert actual["caller_context_restored"] is True
    record_testsuite_property("hostile_context_subprocess_legacy_modules", json.dumps(actual["forbidden_legacy_modules"]))
    record_testsuite_property("hostile_context_caller_restored", "true")


def test_direct_source_outer_flags_and_provenance_reject_before_unusable_raw_data(record_testsuite_property):
    # These intentionally invalid declarations do not pretend to be accepted source data.
    base = {
        "schema": "mhdlab-boundary-saved-source-admission-v1",
        "status": "accepted_saved_joint_polynomial_source", "file_sha256": dict(_CORE.PINNED),
        "parameter_provenance": PROVENANCE,
        "inherited_qualification": "Accepted source interpolation and prior IVP evidence remain distinct from new polynomial quadrature. New query IVP error is not qualified.",
        "host_budget_accepted": False, "owner_implementation_verified": False,
        "documents": {alias: {} for alias in _CORE.PINNED},
    }
    query = BoundaryQuery("synthetic_guard", F(1, 4), F(3, 4), PROVENANCE)
    controls = (
        ({"host_budget_accepted": True}, "invalid source admission record"),
        ({"owner_implementation_verified": True}, "invalid source admission record"),
        ({"parameter_provenance": " "}, "parameter_provenance must be nonempty text"),
    )
    for changes, message in controls:
        with pytest.raises(ValueError, match=message):
            source = AdmittedBoundarySource(json.dumps({**base, **changes}), ())
            build_boundary_quadrature(source, (query,), parameter_status="manufactured", parameter_provenance=PROVENANCE)
    record_testsuite_property("forged_outer_declaration_rejections", str(len(controls)))
