"""Positive quadrature of saved joint emission-time/temperature polynomials.

This manufactured boundary reference makes no transport or new IVP claim.  The
finite Decimal rule, exact polynomial target, and actual binary export are
three separate numerical objects.  No source model is imported or evaluated.
"""

from dataclasses import dataclass, replace
from decimal import Context, Decimal, DivisionByZero, InvalidOperation, Overflow, localcontext, ROUND_HALF_EVEN
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
import json
import math
import re
import sys


KB = Fraction("1.380649e-23")
QE = Fraction("1.602176634e-19")
MAX_JSON_BYTES = 64 * 1024 * 1024
MAX_FRACTION_BITS = 65536
ALIASES = frozenset(("report", "windows", "input_config", "window_config",
                     "independent_completion", "independent_report", "scientific_review"))
QUANTITIES = ("released_number", "local_time_first_moment_s", "local_time_second_moment_s2",
              "temperature_first_moment_k", "temperature_second_moment_k2",
              "local_time_temperature_moment_k_s", "binding_energy_j", "internal_energy_j",
              "translational_energy_j", "gas_energy_j", "total_wall_energy_debit_j")
PHASES = ("target", "finite_rule", "represented", "rule_generation_defect",
          "export_projection_defect", "total_defect")
PINNED = {
    "report": "69907c8ba61545b31e67675a2c71bf8b9453006ea9f705b6f16d02e5b1a233ac",
    "windows": "64e0ce6b9cb271ad423109f3e52db02ab6017eb0d2591bbc8bdb30ca50aba731",
    "input_config": "22f21dae395f31f9a13186821c5364033b50f051aa7296f93839896c2fae4b56",
    "window_config": "2e8f1eb938a89957d2453f801321365ac6513fa8154b4dd3b925c1a67df1b18c",
    "independent_completion": "7c0715cb1a89b7126c391886ef9a1ea2100f7d9be6f50bb0cd4d07aa64ecfb2a",
    "independent_report": "80f6577f67dacb29ad5d1f53c13a85b89aa10ab160cd79daf000dd43ca48d733",
    "scientific_review": "29f5b4e5291eba8e42e46f8f64f41524771f875c6c720898feb1482d7b388de0",
}


def _text(value, name):
    if type(value) is not str or not value.strip():
        raise ValueError(name + " must be nonempty text")
    return value


def _fraction(value, name="quantity"):
    if type(value) not in (int, float, Fraction):
        raise ValueError(name + " must be an integer, Fraction or finite float")
    if type(value) is float and not math.isfinite(value):
        raise ValueError(name + " must be finite")
    result = Fraction(value)
    if max(result.numerator.bit_length(), result.denominator.bit_length()) > MAX_FRACTION_BITS:
        raise ValueError(name + " exceeds the exact arithmetic domain")
    return result


def _ratio(value):
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


def _read_ratio(value):
    if type(value) is not dict or set(value) != {"numerator", "denominator"}:
        raise ValueError("invalid canonical rational record")
    n, d = value["numerator"], value["denominator"]
    if (type(n) is not str or type(d) is not str or len(n) > 20000 or len(d) > 20000
            or not re.fullmatch(r"0|-?[1-9][0-9]*", n)
            or not re.fullmatch(r"[1-9][0-9]*", d)):
        raise ValueError("invalid canonical rational integers")
    result = _fraction(Fraction(int(n), int(d)))
    if _ratio(result) != value:
        raise ValueError("rational record is not reduced")
    return result


def _binary(value):
    try:
        result = float(value)
    except (OverflowError, ValueError) as error:
        raise ValueError("unsupported binary projection overflow") from error
    if not math.isfinite(result) or (value != 0 and result == 0):
        raise ValueError("unsupported nonzero binary projection")
    return result


def _record(value):
    value = _fraction(value)
    binary = _binary(value)
    represented = Fraction(binary)
    return {"exact": _ratio(value), "value": binary, "projection": {
        "exact_zero": value == 0, "positive": value > 0, "negative": value < 0,
        "finite": True, "underflowed": False, "overflowed": False,
        "rounded": represented != value,
        "subnormal": 0 < abs(binary) < sys.float_info.min,
        "exact_rounding_error": _ratio(represented - value)}}


def _exact(record):
    return _read_ratio(record["exact"])


def _canonical(document):
    return json.dumps(document, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _strict_json(raw):
    if type(raw) is not bytes or len(raw) > MAX_JSON_BYTES:
        raise ValueError("JSON byte input exceeds supported size")

    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError("duplicate JSON key")
            result[key] = value
        return result

    def integer(token):
        if token == "-0":
            raise ValueError("negative-zero JSON integer is not canonical")
        if len(token) > 4096:
            raise ValueError("JSON integer exceeds supported size")
        return int(token)

    def floating(token):
        if len(token) > 4096:
            raise ValueError("JSON float exceeds supported size")
        value = float(token)
        if not math.isfinite(value):
            raise ValueError("nonfinite JSON float")
        if value == 0 and Decimal(token) != 0:
            raise ValueError("nonzero JSON float underflows its representation")
        if value == 0 and token.startswith("-"):
            raise ValueError("negative-zero JSON float is not canonical")
        return value

    def constant(_):
        raise ValueError("nonfinite JSON constant")

    try:
        result = json.loads(raw.decode("utf-8"), object_pairs_hook=pairs,
                            parse_int=integer, parse_float=floating, parse_constant=constant)
    except (UnicodeError, RecursionError, json.JSONDecodeError) as error:
        raise ValueError("invalid or excessively nested JSON") from error
    stack = [(result, 0)]
    while stack:
        item, depth = stack.pop()
        if depth > 128:
            raise ValueError("JSON nesting exceeds supported domain")
        if isinstance(item, dict):
            stack.extend((x, depth + 1) for x in item.values())
        elif isinstance(item, list):
            stack.extend((x, depth + 1) for x in item)
    return result


def _controls(values, name):
    if type(values) not in (tuple, list) or len(values) != 4:
        raise ValueError(name + " requires four cubic Bernstein controls")
    return tuple(_fraction(x, name) for x in values)


@dataclass(frozen=True)
class BoundaryPiece:
    start_s: Fraction
    stop_s: Fraction
    emitted_controls_m2: tuple
    temperature_controls_k: tuple
    temperature_min_k: Fraction
    temperature_max_k: Fraction
    initial_coverage_m2: Fraction
    state_id: str
    source_history_id: str
    source_cell_index: int
    parameter_status: str
    parameter_provenance: str

    def __post_init__(self):
        for name in ("start_s", "stop_s", "temperature_min_k", "temperature_max_k", "initial_coverage_m2"):
            object.__setattr__(self, name, _fraction(getattr(self, name), name))
        for name in ("emitted_controls_m2", "temperature_controls_k"):
            object.__setattr__(self, name, _controls(getattr(self, name), name))
        if not self.start_s < self.stop_s or self.initial_coverage_m2 < 0:
            raise ValueError("invalid piece support or initial inventory")
        if not 0 < self.temperature_min_k <= self.temperature_max_k:
            raise ValueError("invalid temperature domain")
        m, t = self.emitted_controls_m2, self.temperature_controls_k
        if not 0 <= m[0] <= m[1] <= m[2] <= m[3] <= self.initial_coverage_m2:
            raise ValueError("emitted controls must be nonnegative monotone within inventory")
        if not all(self.temperature_min_k <= x <= self.temperature_max_k for x in t):
            raise ValueError("temperature controls leave the supported domain")
        if type(self.source_cell_index) is not int or self.source_cell_index < 0:
            raise ValueError("source_cell_index must be a nonnegative integer")
        for name in ("state_id", "source_history_id", "parameter_provenance"):
            _text(getattr(self, name), name)
        if self.parameter_status != "manufactured":
            raise ValueError("only explicit manufactured pieces are supported")


def _piece_copy(piece):
    if type(piece) is not BoundaryPiece:
        raise ValueError("expected BoundaryPiece")
    return replace(piece)


def _split(values, coordinate):
    rows = [tuple(values)]
    while len(rows[-1]) > 1:
        last = rows[-1]
        rows.append(tuple((1-coordinate)*a + coordinate*b for a, b in zip(last, last[1:])))
    return tuple(row[0] for row in rows), tuple(row[-1] for row in reversed(rows))


def restrict_boundary_piece(piece, *, start_s, stop_s):
    piece = _piece_copy(piece)
    start, stop = _fraction(start_s), _fraction(stop_s)
    if not piece.start_s <= start < stop <= piece.stop_s:
        raise ValueError("restriction must have positive finite contained support")
    a = (start-piece.start_s)/(piece.stop_s-piece.start_s)
    b = (stop-piece.start_s)/(piece.stop_s-piece.start_s)

    def restrict(values):
        left = _split(values, b)[0]
        return _split(left, a/b)[1] if a else left

    return replace(piece, start_s=start, stop_s=stop,
                   emitted_controls_m2=restrict(piece.emitted_controls_m2),
                   temperature_controls_k=restrict(piece.temperature_controls_k))


def _evaluate(values, coordinate):
    while len(values) > 1:
        values = tuple((1-coordinate)*a + coordinate*b for a, b in zip(values, values[1:]))
    return values[0]


def _product(a, b):
    m, n = len(a)-1, len(b)-1
    result = []
    for k in range(m+n+1):
        result.append(sum((a[i]*b[k-i]*Fraction(math.comb(m, i)*math.comb(n, k-i),
                                              math.comb(m+n, k))
                           for i in range(max(0, k-n), min(m, k)+1)), Fraction()))
    return tuple(result)


def _integral(values):
    return sum(values, Fraction())/len(values)


@dataclass(frozen=True)
class GL5Rule:
    nodes: tuple
    weights: tuple

    def __post_init__(self):
        if (type(self.nodes) is not tuple or type(self.weights) is not tuple
                or len(self.nodes) != 5 or len(self.weights) != 5):
            raise ValueError("GL5 requires five nodes and weights")
        if (any(type(x) is not Fraction for x in self.nodes + self.weights)
                or not all(0 < x < 1 for x in self.nodes)
                or tuple(sorted(set(self.nodes))) != self.nodes
                or not all(w > 0 for w in self.weights)):
            raise ValueError("invalid positive GL5 rule")

    def to_dict(self):
        return {"schema": "mhdlab-finite-decimal-gl5-rule-v1",
                "method": "closed_form_GL5_Decimal80_then_exact_rational",
                "precision": 80, "rounding": "ROUND_HALF_EVEN", "Emin": -10000, "Emax": 10000,
                "nodes": [_ratio(x) for x in self.nodes], "weights": [_ratio(w) for w in self.weights],
                "monomial_moment_defects": [_record(sum((w*x**k for x, w in zip(self.nodes, self.weights)), Fraction())
                                                       - Fraction(1, k+1)) for k in range(10)],
                "qualification": "Finite approximate rule; no renormalization or certified enclosure."}


@lru_cache(maxsize=1)
def _rule_parameters():
    deterministic = Context(prec=80, rounding=ROUND_HALF_EVEN, Emin=-10000, Emax=10000,
                            capitals=1, clamp=0, flags=[],
                            traps=[InvalidOperation, DivisionByZero, Overflow])
    with localcontext(deterministic):
        s = (Decimal(10)/7).sqrt()
        inner = (5-2*s).sqrt()/3
        outer = (5+2*s).sqrt()/3
        nodes = tuple(Fraction(x) for x in ((1-outer)/2, (1-inner)/2, Decimal("0.5"),
                                           (1+inner)/2, (1+outer)/2))
        w0 = (322-13*Decimal(70).sqrt())/1800
        w1 = (322+13*Decimal(70).sqrt())/1800
        weights = tuple(Fraction(x) for x in (w0, w1, Decimal(64)/225, w1, w0))
    return nodes, weights


def gauss_legendre5_rule():
    return GL5Rule(*_rule_parameters())


@dataclass(frozen=True)
class BoundaryNode:
    node_id: str
    node_index: int
    coordinate: Fraction
    unprojected_number_weight_m2: Fraction
    unprojected_time_s: Fraction
    unprojected_temperature_k: Fraction
    unprojected_expected_translation_j: Fraction
    number_weight_m2: float
    time_s: float
    temperature_k: float
    expected_translational_energy_j: float
    desorption_energy_j: float
    internal_energy_j: float


@dataclass(frozen=True)
class PieceQuadrature:
    piece: BoundaryPiece
    rule: GL5Rule
    nodes: tuple
    desorption_energy_j: Fraction
    internal_energy_j: Fraction
    skip_reason: str | None
    _fingerprint: str


def _quadrature_fingerprint(piece, rule, nodes, d, internal, skip):
    payload = [_piece_identity(piece), [_ratio(x) for x in rule.nodes],
               [_ratio(x) for x in rule.weights], _ratio(d), _ratio(internal), skip,
               [[str(getattr(node, key)) for key in node.__dataclass_fields__] for node in nodes]]
    return sha256(_canonical(payload).encode()).hexdigest()


def _quadrature_record(piece, rule, nodes, d, internal, skip):
    return PieceQuadrature(piece, rule, nodes, d, internal, skip,
                           _quadrature_fingerprint(piece, rule, nodes, d, internal, skip))


def _piece_identity(piece):
    return {"source_history_id": piece.source_history_id, "state_id": piece.state_id,
            "source_cell_index": piece.source_cell_index, "start_exact_s": _ratio(piece.start_s),
            "stop_exact_s": _ratio(piece.stop_s),
            "emitted_controls_m2": [_ratio(x) for x in piece.emitted_controls_m2],
            "temperature_controls_k": [_ratio(x) for x in piece.temperature_controls_k]}


def make_piece_quadrature(piece, *, desorption_energy_ev, internal_energy_ev, rule):
    piece = _piece_copy(piece)
    if type(rule) is not GL5Rule or rule != gauss_legendre5_rule():
        raise ValueError("only the deterministic finite Decimal80 GL5 rule is supported")
    d, internal = _fraction(desorption_energy_ev)*QE, _fraction(internal_energy_ev)*QE
    if d < 0 or internal < 0:
        raise ValueError("per-parent costs must be nonnegative")
    if piece.emitted_controls_m2[-1] == piece.emitted_controls_m2[0]:
        return _quadrature_record(piece, rule, (), d, internal, "exact_zero_emitted_mass")
    db, ib = _binary(d), _binary(internal)
    derivative = tuple(3*(b-a) for a, b in zip(piece.emitted_controls_m2, piece.emitted_controls_m2[1:]))
    identity = sha256(_canonical(_piece_identity(piece)).encode()).hexdigest()
    nodes = []
    for index, (coordinate, weight) in enumerate(zip(rule.nodes, rule.weights)):
        amount = weight*_evaluate(derivative, coordinate)
        t = piece.start_s+(piece.stop_s-piece.start_s)*coordinate
        temperature = _evaluate(piece.temperature_controls_k, coordinate)
        energy = 2*KB*temperature
        if amount <= 0:
            raise ValueError("positive piece requires positive interior rule amounts")
        node = BoundaryNode(identity+":"+str(index), index, coordinate, amount, t, temperature, energy,
                            _binary(amount), _binary(t), _binary(temperature), _binary(energy), db, ib)
        if not piece.start_s < Fraction(node.time_s) < piece.stop_s:
            raise ValueError("represented node time must remain strictly interior")
        if not piece.temperature_min_k <= Fraction(node.temperature_k) <= piece.temperature_max_k:
            raise ValueError("represented temperature leaves the declared domain")
        nodes.append(node)
    if len({x.time_s for x in nodes}) != 5:
        raise ValueError("represented node times are not distinct")
    return _quadrature_record(piece, rule, tuple(nodes), d, internal, None)


def _channels(mass, time1, time2, temp1, temp2, time_temp, binding, internal, translation):
    return dict(zip(QUANTITIES, (mass, time1, time2, temp1, temp2, time_temp,
                                binding, internal, translation, internal+translation,
                                binding+internal+translation)))


def _target(piece, anchor, area, d, internal):
    derivative = tuple(3*(b-a) for a, b in zip(piece.emitted_controls_m2, piece.emitted_controls_m2[1:]))
    tau = (piece.start_s-anchor, piece.stop_s-anchor)
    temperature = piece.temperature_controls_k
    integrands = (derivative, _product(derivative, tau), _product(_product(derivative, tau), tau),
                  _product(derivative, temperature), _product(_product(derivative, temperature), temperature),
                  _product(_product(derivative, tau), temperature))
    m, t1, t2, k1, k2, tk = (area*_integral(x) for x in integrands)
    return _channels(m, t1, t2, k1, k2, tk, m*d, m*internal, 2*KB*k1)


def _sums(nodes, anchor, area, d, internal, represented):
    values = {key: Fraction() for key in QUANTITIES}
    for node in nodes:
        if represented:
            w, tau, temperature = (Fraction(_binary(area*node.unprojected_number_weight_m2)),
                                    Fraction(node.time_s)-anchor, Fraction(node.temperature_k))
            bd, bi, bx = map(Fraction, (node.desorption_energy_j, node.internal_energy_j,
                                      node.expected_translational_energy_j))
        else:
            w, tau, temperature = (area*node.unprojected_number_weight_m2,
                                    node.unprojected_time_s-anchor, node.unprojected_temperature_k)
            bd, bi, bx = d, internal, node.unprojected_expected_translation_j
        increments = _channels(w, w*tau, w*tau*tau, w*temperature, w*temperature*temperature,
                               w*tau*temperature, w*bd, w*bi, w*bx)
        for key in values:
            values[key] += increments[key]
    return values


def _statistics(values):
    names = ("mean_local_time_s", "variance_local_time_s2", "mean_temperature_k",
             "variance_temperature_k2", "mean_local_time_temperature_k_s",
             "covariance_local_time_temperature_k_s")
    m = values["released_number"]
    if m == 0:
        return {key: {"available": False, "reason": "exact_zero_emitted_mass", "record": None} for key in names}
    if m < 0:
        raise ValueError("negative emitted mass in summary")
    t, t2, k, k2, tk = (values[key]/m for key in QUANTITIES[1:6])
    if t2-t*t < 0 or k2-k*k < 0:
        raise ValueError("negative exact variance in positive measure")
    return {key: {"available": True, "reason": None, "record": _record(value)}
            for key, value in zip(names, (t, t2-t*t, k, k2-k*k, tk, tk-t*k))}


def _finish_summary(anchor, area, nodes, target, finite, represented, consistency, skipped):
    phases = {"target": target, "finite_rule": finite, "represented": represented,
              "rule_generation_defect": {k: finite[k]-target[k] for k in QUANTITIES},
              "export_projection_defect": {k: represented[k]-finite[k] for k in QUANTITIES},
              "total_defect": {k: represented[k]-target[k] for k in QUANTITIES}}
    return {"local_time_anchor_exact_s": _ratio(anchor), "area_exact_m2": None if area is None else _ratio(area),
            "nodes": nodes, "node_count": len(nodes),
            **{phase: {k: _record(v) for k, v in values.items()} for phase, values in phases.items()},
            "statistics": {phase: _statistics(phases[phase]) for phase in ("target", "finite_rule", "represented")},
            "energy_consistency": {k: _record(v) for k, v in consistency.items()}, "skipped_pieces": skipped}


def summarize_piece_quadrature(quadrature, *, local_time_anchor_s, area_m2):
    if type(quadrature) is not PieceQuadrature:
        raise ValueError("expected PieceQuadrature")
    piece = _piece_copy(quadrature.piece)
    anchor, area = _fraction(local_time_anchor_s), _fraction(area_m2)
    if area <= 0:
        raise ValueError("area must be strictly positive")
    d, internal = quadrature.desorption_energy_j, quadrature.internal_energy_j
    # Revalidate an immutable record even if its fields were changed through a bypass.
    if (_quadrature_fingerprint(piece, quadrature.rule, quadrature.nodes, d, internal,
                                quadrature.skip_reason) != quadrature._fingerprint):
        raise ValueError("quadrature record does not match its source piece")
    nodes = []
    for node in quadrature.nodes:
        tau = Fraction(node.time_s)-anchor
        temperature_energy = 2*KB*Fraction(node.temperature_k)
        nodes.append({"node_id": node.node_id, "node_index": node.node_index,
                      "number_weight": _binary(area*node.unprojected_number_weight_m2),
                      "time_s": node.time_s, "temperature_k": node.temperature_k,
                      "expected_translational_energy_j": node.expected_translational_energy_j,
                      "desorption_energy_j": node.desorption_energy_j, "internal_energy_j": node.internal_energy_j,
                      "local_time_exact_s": _ratio(tau),
                      "translation_from_represented_temperature_j": _record(temperature_energy),
                      "translation_field_minus_temperature_energy_j": _record(Fraction(node.expected_translational_energy_j)-temperature_energy)})
    target = _target(piece, anchor, area, d, internal)
    finite = _sums(quadrature.nodes, anchor, area, d, internal, False)
    represented = _sums(quadrature.nodes, anchor, area, d, internal, True)
    temp_energy = 2*KB*represented["temperature_first_moment_k"]
    consistency = {"represented_translation_from_temperature_j": temp_energy,
                   "represented_translation_field_minus_temperature_energy_j": represented["translational_energy_j"]-temp_energy,
                   "represented_binding_field_minus_supplied_cost_j": represented["binding_energy_j"]-d*represented["released_number"],
                   "represented_internal_field_minus_supplied_cost_j": represented["internal_energy_j"]-internal*represented["released_number"]}
    skipped = [] if quadrature.skip_reason is None else [{**_piece_identity(piece), "reason": quadrature.skip_reason}]
    return _finish_summary(anchor, area, nodes, target, finite, represented, consistency, skipped)


def _shift(values, delta):
    result = dict(values)
    m, t, k = values[QUANTITIES[0]], values[QUANTITIES[1]], values[QUANTITIES[3]]
    result[QUANTITIES[1]] += delta*m
    result[QUANTITIES[2]] += 2*delta*t+delta*delta*m
    result[QUANTITIES[5]] += delta*k
    return result


def combine_boundary_summaries(summaries, *, local_time_anchor_s):
    anchor = _fraction(local_time_anchor_s)
    summaries = tuple(summaries)
    if not summaries:
        raise ValueError("at least one boundary summary is required")
    totals = {phase: {k: Fraction() for k in QUANTITIES} for phase in ("target", "finite_rule", "represented")}
    nodes, skipped, consistency, areas = [], [], {}, []
    for summary in summaries:
        old_anchor = _read_ratio(summary["local_time_anchor_exact_s"])
        for phase in totals:
            values = {k: _exact(summary[phase][k]) for k in QUANTITIES}
            for k, value in _shift(values, old_anchor-anchor).items():
                totals[phase][k] += value
        for node in summary["nodes"]:
            copied = dict(node)
            copied["local_time_exact_s"] = _ratio(Fraction(node["time_s"])-anchor)
            nodes.append(copied)
        skipped.extend(summary["skipped_pieces"])
        for k, record in summary["energy_consistency"].items():
            consistency[k] = consistency.get(k, Fraction())+_exact(record)
        areas.append(summary.get("area_exact_m2"))
    result = _finish_summary(anchor, None, nodes, totals["target"], totals["finite_rule"],
                             totals["represented"], consistency, skipped)
    result["component_areas_exact_m2"] = areas
    return result


@dataclass(frozen=True)
class BoundaryQuery:
    window_id: str
    start_fraction: Fraction
    stop_fraction: Fraction
    parameter_provenance: str

    def __post_init__(self):
        _text(self.window_id, "window_id")
        _text(self.parameter_provenance, "parameter_provenance")
        for name in ("start_fraction", "stop_fraction"):
            object.__setattr__(self, name, _fraction(getattr(self, name), name))
        if not 0 <= self.start_fraction < self.stop_fraction <= 1:
            raise ValueError("query fractions must have positive contained support")


@dataclass(frozen=True)
class AdmittedBoundarySource:
    _serialized: str
    _raw_files: tuple

    def to_dict(self):
        return json.loads(self._serialized)


@dataclass(frozen=True)
class BoundaryQuadrature:
    _serialized: str

    def to_dict(self):
        return json.loads(self._serialized)


def _require(document, key, expected):
    actual = document.get(key)
    if type(actual) is not type(expected) or actual != expected:
        raise ValueError("source association mismatch: " + key)


def _admission_metadata(parameter_provenance):
    return {"schema": "mhdlab-boundary-saved-source-admission-v1",
            "status": "accepted_saved_joint_polynomial_source", "file_sha256": dict(PINNED),
            "parameter_provenance": _text(parameter_provenance, "parameter_provenance"),
            "inherited_qualification": "Accepted source interpolation and prior IVP evidence remain distinct from new polynomial quadrature. New query IVP error is not qualified.",
            "host_budget_accepted": False, "owner_implementation_verified": False}


def admit_boundary_source(files, *, expected_sha256, parameter_provenance):
    _text(parameter_provenance, "parameter_provenance")
    if set(files) != ALIASES or set(expected_sha256) != ALIASES:
        raise ValueError("source aliases must match the seven required files")
    if sum(len(raw) for raw in files.values() if type(raw) is bytes) > 256*1024*1024:
        raise ValueError("source byte bundle exceeds supported size")
    for alias in sorted(ALIASES):
        pin, raw = expected_sha256[alias], files[alias]
        if type(pin) is not str or not re.fullmatch("[0-9a-f]{64}", pin) or type(raw) is not bytes:
            raise ValueError("invalid source sha256 or bytes")
        if sha256(raw).hexdigest() != pin:
            raise ValueError("source sha256 mismatch: " + alias)
    documents = {alias: _strict_json(files[alias]) for alias in sorted(ALIASES)}
    if any(type(d) is not dict for d in documents.values()):
        raise ValueError("source documents must be JSON objects")
    schemas = {
        "report": "mhdlab-coupled-emission-windows-export-run-v1",
        "windows": "mhdlab-coupled-emission-windows-v1",
        "input_config": "mhdlab-thermal-width-power-manufactured-config-v1",
        "window_config": "mhdlab-coupled-emission-windows-manufactured-config-v1",
        "independent_completion": "mitl-independent-coupled-window-completion-v1",
        "independent_report": "mitl-independent-coupled-window-query-report-v1",
        "scientific_review": "mitl-coupled-window-static-oracle-review-v1"}
    for alias, schema in schemas.items():
        _require(documents[alias], "schema", schema)
    # This reference admits one explicitly accepted immutable source bundle.
    # Arbitrary rebound hashes are insufficient evidence of source acceptance.
    if dict(expected_sha256) != PINNED:
        raise ValueError("source byte pins do not identify the accepted source bundle")
    report, windows = documents["report"], documents["windows"]
    completion, reference, review = (documents[x] for x in ("independent_completion", "independent_report", "scientific_review"))
    _require(report, "result_sha256", PINNED["windows"])
    _require(report, "config_sha256", PINNED["window_config"])
    _require(windows, "parameter_status", "manufactured")
    _require(windows, "host_budget_accepted", False)
    _require(windows, "owner_implementation_verified", False)
    for document in (completion, reference):
        _require(document, "status", "passed")
        _require(document, "all_checks_passed", True)
        _require(document, "producer_report_sha256", PINNED["report"])
    _require(completion, "report_sha256", PINNED["independent_report"])
    _require(completion, "final_static_review_sha256", PINNED["scientific_review"])
    _require(review, "status", "accepted")
    for key in ("actual_report_sha256", "owner_report_sha256"):
        _require(review, key, PINNED["independent_report"])
    _require(review, "window_document_sha256", PINNED["windows"])
    _require(review, "further_changes_requested", False)
    if any(item.get("status") != "repaired_before_first_query" for item in review["findings"]):
        raise ValueError("unresolved source acceptance finding")
    config = documents["input_config"]
    complete = [x for x in windows["source_windows"] if x["query_id"] == "complete_source_support"]
    if len(complete) != 3 or len({x["source_history_id"] for x in complete}) != 3:
        raise ValueError("three complete source histories are required")
    state = config["common_scalar_problem"]["reservoir"]["state_id"]
    duration = _fraction(config["duration_s"])
    for row in complete:
        _require(row, "state_id", state)
        if (_read_ratio(row["window_start_exact_s"]) != 0
                or _read_ratio(row["window_stop_exact_s"]) != duration or len(row["pieces"]) != 64):
            raise ValueError("incomplete constant-power source support")
        previous = None
        for index, item in enumerate(row["pieces"]):
            piece = _saved_piece(item, row, config)
            if piece.source_cell_index != index or piece.start_s != duration*index/64 or piece.stop_s != duration*(index+1)/64:
                raise ValueError("source cell support association mismatch")
            if item["emitted_controls_m2"] != item["source_emitted_controls_m2"] or item["temperature_controls_k"] != item["source_temperature_controls_k"]:
                raise ValueError("source complete controls are restricted")
            if previous is not None and (previous.emitted_controls_m2[-1] != piece.emitted_controls_m2[0]
                                         or previous.temperature_controls_k[-1] != piece.temperature_controls_k[0]):
                raise ValueError("source piece endpoints are discontinuous")
            previous = piece
    admitted = {**_admission_metadata(parameter_provenance), "documents": documents}
    return AdmittedBoundarySource(_canonical(admitted), tuple((alias, files[alias]) for alias in sorted(ALIASES)))


def _saved_piece(item, row, config):
    problem = config["common_scalar_problem"]
    return BoundaryPiece(_read_ratio(item["start_exact_s"]), _read_ratio(item["stop_exact_s"]),
                         tuple(map(_read_ratio, item["emitted_controls_m2"])),
                         tuple(map(_read_ratio, item["temperature_controls_k"])),
                         problem["temperature_min_k"], problem["temperature_max_k"], problem["initial_coverage_m2"],
                         row["state_id"], row["source_history_id"], item["source_cell_index"], "manufactured",
                         row["parameter_provenance"])


def _geometry_ratio(record):
    ratio = _read_ratio({k: record[k] for k in ("numerator", "denominator")})
    if _fraction(record["value"]) != ratio:
        raise ValueError("geometry value differs from its exact rational")
    return ratio


def _node_cache(node, piece_id):
    return {"node_id": node.node_id, "piece_id": piece_id, "node_index": node.node_index,
            "coordinate_exact": _ratio(node.coordinate),
            "unprojected": {"number_weight_m2": _ratio(node.unprojected_number_weight_m2),
                            "time_s": _ratio(node.unprojected_time_s),
                            "temperature_k": _ratio(node.unprojected_temperature_k),
                            "expected_translational_energy_j": _ratio(node.unprojected_expected_translation_j)},
            "binary": {key: getattr(node, key) for key in ("number_weight_m2", "time_s", "temperature_k",
                                                          "expected_translational_energy_j", "desorption_energy_j", "internal_energy_j")},
            "energy_consistency": {"translation_field_minus_temperature_energy_j":
                                   _record(Fraction(node.expected_translational_energy_j)-2*KB*Fraction(node.temperature_k))}}


def build_boundary_quadrature(source, queries, *, parameter_status, parameter_provenance):
    if type(source) is not AdmittedBoundarySource or parameter_status != "manufactured":
        raise ValueError("accepted source and explicit manufactured status are required")
    _text(parameter_provenance, "parameter_provenance")
    admitted = source.to_dict()
    metadata = _admission_metadata(admitted.get("parameter_provenance"))
    if (set(admitted) != set(metadata) | {"documents"}
            or any(_canonical(admitted[k]) != _canonical(v) for k, v in metadata.items())):
        raise ValueError("invalid source admission record")
    documents = admitted["documents"]
    if type(documents) is not dict or set(documents) != ALIASES:
        raise ValueError("invalid admitted document aliases")
    # Immutable byte association is replayed without source evaluation.
    if (type(source._raw_files) is not tuple or len(source._raw_files) != len(ALIASES)
            or {a for a, _ in source._raw_files} != ALIASES):
        raise ValueError("invalid admitted byte association")
    for alias, raw in source._raw_files:
        if type(raw) is not bytes or sha256(raw).hexdigest() != PINNED[alias]:
            raise ValueError("admitted source bytes changed")
        if _canonical(_strict_json(raw)) != _canonical(documents[alias]):
            raise ValueError("admitted document differs from its original bytes")
    queries = tuple(replace(query) if type(query) is BoundaryQuery else None for query in queries)
    if not queries or any(q is None for q in queries) or len({q.window_id for q in queries}) != len(queries):
        raise ValueError("queries must have unique identities")
    if len(queries) > 16:
        raise ValueError("too many queries for the bounded reference")
    config, windows = documents["input_config"], documents["windows"]
    problem, reservoir = config["common_scalar_problem"], config["common_scalar_problem"]["reservoir"]
    d, internal = _fraction(reservoir["desorption_energy_ev"]), _fraction(reservoir["internal_energy_ev"])
    duration = _fraction(config["duration_s"])
    edges = [_geometry_ratio(x) for x in config["geometry"]["width_edges_m"]]
    length = _geometry_ratio(config["geometry"]["common_length_m"])
    areas = [(b-a)*length for a, b in zip(edges, edges[1:])]
    if len(areas) != 3 or min(areas) <= 0:
        raise ValueError("three positive exact patch areas are required")
    rule = gauss_legendre5_rule()
    complete = sorted((x for x in windows["source_windows"] if x["query_id"] == "complete_source_support"), key=lambda x: x["source_history_id"])
    quadratures, cache_entries, node_entries = {}, [], []
    source_queries, patch_queries, aggregate_queries, contrasts = [], [], [], []
    source_map = {}
    for row in complete:
        for query in queries:
            start, stop = duration*query.start_fraction, duration*query.stop_fraction
            pieces, summaries, piece_ids = [], [], []
            for saved in row["pieces"]:
                piece = _saved_piece(saved, row, config)
                a, b = max(start, piece.start_s), min(stop, piece.stop_s)
                if a >= b:
                    continue
                restricted = restrict_boundary_piece(piece, start_s=a, stop_s=b)
                key = (row["source_history_id"], piece.source_cell_index, a, b)
                if key not in quadratures:
                    quadrature = make_piece_quadrature(restricted, desorption_energy_ev=d,
                                                       internal_energy_ev=internal, rule=rule)
                    quadratures[key] = quadrature
                    piece_id = sha256(_canonical(_piece_identity(restricted)).encode()).hexdigest()
                    cache_entries.append({"piece_id": piece_id, **_piece_identity(restricted),
                                          "source_start_exact_s": _ratio(piece.start_s), "source_stop_exact_s": _ratio(piece.stop_s),
                                          "source_emitted_controls_m2": [_ratio(x) for x in piece.emitted_controls_m2],
                                          "source_temperature_controls_k": [_ratio(x) for x in piece.temperature_controls_k],
                                          "node_ids": [n.node_id for n in quadrature.nodes], "skip_reason": quadrature.skip_reason})
                    node_entries.extend(_node_cache(node, piece_id) for node in quadrature.nodes)
                quadrature = quadratures[key]
                piece_id = sha256(_canonical(_piece_identity(restricted)).encode()).hexdigest()
                pieces.append(quadrature)
                piece_ids.append(piece_id)
                summaries.append(summarize_piece_quadrature(quadrature, local_time_anchor_s=start, area_m2=1))
            summary = combine_boundary_summaries(summaries, local_time_anchor_s=start)
            summary["area_exact_m2"] = _ratio(Fraction(1))
            key = row["source_history_id"]+"/"+query.window_id
            for node in summary["nodes"]:
                node.update({"view_id": key+"/"+node["node_id"], "source_history_id": row["source_history_id"], "query_id": query.window_id})
            result = {"source_history_id": row["source_history_id"], "query_id": query.window_id,
                      "state_id": row["state_id"], "window_start_exact_s": _ratio(start), "window_stop_exact_s": _ratio(stop),
                      "local_time_anchor_exact_s": _ratio(start), "piece_ids": piece_ids,
                      "node_ids": [n["node_id"] for n in summary["nodes"]], "summary": summary,
                      "canonical_history_sha256": row["canonical_history_sha256"],
                      "raw_history_file_sha256": row["raw_history_file_sha256"],
                      "parent_complete_source_window_key": row["source_window_key"],
                      "parent_source_RK_comparison": row["source_RK_comparison"],
                      "parent_independent_window_reference": row["independent_window_reference"],
                      "new_query_IVP_error_qualified": False}
            source_queries.append(result)
            source_map[key] = (result, pieces)
    for query in queries:
        start = duration*query.start_fraction
        aggregates = {}
        for hypothesis, assignments in config["hypotheses"].items():
            indices, summaries = [], []
            for index, case_id in enumerate(assignments):
                source_key = case_id+"/"+query.window_id
                source_row, pieces = source_map[source_key]
                summary = combine_boundary_summaries([summarize_piece_quadrature(x, local_time_anchor_s=start, area_m2=areas[index])
                                                       for x in pieces], local_time_anchor_s=start)
                summary["area_exact_m2"] = _ratio(areas[index])
                for node in summary["nodes"]:
                    node.update({"view_id": hypothesis+"/"+str(index)+"/"+query.window_id+"/"+node["node_id"],
                                 "hypothesis": hypothesis, "patch_index": index,
                                 "source_history_id": case_id, "query_id": query.window_id})
                indices.append(len(patch_queries))
                patch_queries.append({"hypothesis": hypothesis, "patch_index": index, "source_history_id": case_id,
                                      "query_id": query.window_id, "state_id": source_row["state_id"],
                                      "geometry": {"width_start_exact_m": _ratio(edges[index]), "width_stop_exact_m": _ratio(edges[index+1]),
                                                   "length_exact_m": _ratio(length), "area_exact_m2": _ratio(areas[index])},
                                      "source_query_key": source_key, "summary": summary})
                summaries.append(summary)
            combined = combine_boundary_summaries(summaries, local_time_anchor_s=start)
            combined["area_exact_m2"] = _ratio(sum(areas, Fraction()))
            aggregate_queries.append({"hypothesis": hypothesis, "query_id": query.window_id, "state_id": reservoir["state_id"],
                                      "patch_query_indices": indices, "summary": combined})
            aggregates[hypothesis] = combined
        uniform, nonuniform = aggregates["uniform"], aggregates["nonuniform"]
        stats = {}
        for phase in ("target", "finite_rule", "represented"):
            stats[phase] = {}
            for key in uniform["statistics"][phase]:
                u, n = uniform["statistics"][phase][key], nonuniform["statistics"][phase][key]
                stats[phase][key] = ({"available": True, "reason": None, "record": _record(_exact(n["record"])-_exact(u["record"]))}
                                     if u["available"] and n["available"] else
                                     {"available": False, "reason": "one_hypothesis_has_zero_emitted_mass", "record": None})
        contrasts.append({"query_id": query.window_id,
                          "phase_difference_nonuniform_minus_uniform": {phase: {k: _record(_exact(nonuniform[phase][k])-_exact(uniform[phase][k])) for k in QUANTITIES} for phase in PHASES},
                          "statistics_difference_nonuniform_minus_uniform": stats,
                          "energy_consistency_difference": {k: _record(_exact(nonuniform["energy_consistency"][k])-_exact(uniform["energy_consistency"][k])) for k in uniform["energy_consistency"]}})
    document = {"schema": "mhdlab-boundary-quadrature-v1", "model_identity": "manufactured-positive-joint-boundary-quadrature-v1",
                "parameter_status": parameter_status, "parameter_provenance": parameter_provenance,
                "source_admission": {k: v for k, v in admitted.items() if k != "documents"}, "rule": rule.to_dict(),
                "queries": [{"window_id": q.window_id, "start_fraction": _ratio(q.start_fraction), "stop_fraction": _ratio(q.stop_fraction),
                             "parameter_provenance": q.parameter_provenance} for q in queries],
                "piece_cache": cache_entries, "node_cache": node_entries, "source_queries": source_queries,
                "patch_queries": patch_queries, "aggregate_queries": aggregate_queries, "contrasts": contrasts,
                "cache": {"unique_piece_supports": len(cache_entries), "unique_nodes": len(node_entries),
                          "source_queries": len(source_queries), "source_query_node_references": sum(x["summary"]["node_count"] for x in source_queries),
                          "patch_queries": len(patch_queries), "patch_node_references": sum(x["summary"]["node_count"] for x in patch_queries),
                          "aggregates": len(aggregate_queries), "contrasts": len(contrasts),
                          "source_evolution_calls": 0, "RHS_calls": 0, "assembly_calls": 0, "old_window_export_calls": 0},
                "quantity_units": dict(zip(QUANTITIES, ("parents", "parent s", "parent s2", "parent K", "parent K2", "parent K s", "J", "J", "J", "J", "J"))),
                "limits": ["Polynomial source only; inherited interpolation, source RK and IVP qualifications remain separate.",
                           "Finite Decimal80 rule and binary export errors are explicit, without renormalization.",
                           "No velocity, position or species-mass sampling; no transport or owner acceptance.",
                           "Local time is represented global time minus the exact query anchor."],
                "host_budget_accepted": False, "owner_implementation_verified": False, "new_query_IVP_error_qualified": False}
    return BoundaryQuadrature(_canonical(document))
