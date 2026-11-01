"""Pinned H2 ground-state dissociation rate reference; no event admission.

Raw MCCC arrays stay outside redistribution artifacts. This module reads a
specific primary table. Its interpolation is a declared numerical operation,
not a measured differential collision model or full H2 encounter network.
"""
from dataclasses import dataclass
from decimal import (Decimal, Context, localcontext, ROUND_HALF_EVEN,
                     InvalidOperation, DivisionByZero, Overflow)
import hashlib
import json
import math
from pathlib import Path
import re

SOURCE_SHA256 = '9bf47765e4e381e272b1ff4fa03108c06b95b655584912a426893cb5f0f0a732'
BLOB_SHA256 = '4ebf093eaf5a7fed50504cf4410b04b6183463faa59d67736da16777f348cafc'
MANIFEST_SHA256 = 'decc4f78e3f6ebd8c51eda9e662bb72d36597124c0b7026bea46dbc31c9ad9bd'
SPEC_SHA256 = '805a75e99fa4fb3f23a1c12c10c98cf8fd2f39f8ea2c552387230f7ed8867846'
POLICY_SHA256 = 'd9fb6769f7a1fcc9bd6023bcbee2122ccc2f9eefc658e2da38ccd7495954bf09'
REACTION = '# This file: e + H2(X1Sg,vi=0) -> e + H2(b3Su) (dissociative excitation)'
UNITS = '# Energy (eV)   CS (a_0^2)'
NUMBER = re.compile(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?\Z')
_CAPABILITY = object()


def _context(precision):
    # No caller rounding, exponent, trap, or sticky-flag inheritance.
    return Context(prec=precision, rounding=ROUND_HALF_EVEN,
                   Emin=-999999, Emax=999999, capitals=1, clamp=0,
                   flags=[], traps=[InvalidOperation, DivisionByZero, Overflow])


class Refusal(ValueError):
    pass


def _need(condition, reason):
    if not condition:
        raise Refusal(reason)


def _read_pinned(path, expected, maximum):
    path = Path(path)
    with path.open('rb') as stream:
        raw = stream.read(maximum + 1)
    _need(len(raw) <= maximum, 'source exceeds bounded reader size')
    _need(hashlib.sha256(raw).hexdigest() == expected, 'pinned source identity differs')
    return raw


def _number(text):
    _need(isinstance(text, str) and len(text) <= 80 and NUMBER.fullmatch(text), 'invalid finite numeric token')
    value = Decimal(text)
    _need(value.is_finite(), 'nonfinite numeric token')
    binary = float(value)
    _need(math.isfinite(binary) and (value == 0 or binary != 0), 'numeric token outside finite binary64 range')
    return value


@dataclass(frozen=True)
class _Parsed:
    threshold: Decimal
    rows: tuple


def _parse_source(raw):
    """Bounded syntax checks, used on manufactured parser fixtures in tests."""
    _need(isinstance(raw, bytes) and len(raw) <= 65536, 'bounded source bytes required')
    try:
        lines = raw.decode('ascii').splitlines()
    except UnicodeDecodeError as exc:
        raise Refusal('source must be ASCII') from exc
    _need(lines.count(REACTION) == 1, 'wrong or duplicate reaction/state metadata')
    _need(lines.count(UNITS) == 1, 'wrong or duplicate energy/area units')
    threshold_lines = [s for s in lines if s.startswith('# Threshold:')]
    _need(len(threshold_lines) == 1, 'missing or duplicate threshold')
    fields = threshold_lines[0].split()
    _need(len(fields) == 4 and fields[-1] == 'eV', 'threshold units/columns')
    threshold = _number(fields[2])
    _need(threshold > 0, 'nonpositive threshold')
    rows = []
    for line in lines:
        if not line.strip() or line.startswith('#'):
            continue
        fields = line.split()
        _need(len(fields) == 2, 'cross-section row must have two columns')
        energy, sigma = map(_number, fields)
        _need(energy > 0 and sigma >= 0, 'invalid energy or negative cross section')
        _need(not rows or energy > rows[-1][0], 'energy grid is not strictly increasing')
        rows.append((energy, sigma))
        _need(len(rows) <= 4096, 'too many table rows')
    _need(len(rows) >= 2, 'insufficient table rows')
    _need(rows[0] == (threshold, Decimal(0)), 'threshold must match first zero anchor')
    _need(any(s > 0 for _, s in rows), 'no positive table data')
    return _Parsed(threshold, tuple(rows))


@dataclass(frozen=True)
class RegistryOperands:
    electron_mass_eV: float
    target_mass_eV: float
    hydrogen_mass_eV: float
    target_reference_eV: float
    hydrogen_reference_eV: float
    bohr_m: float

    def __post_init__(self):
        vals = (self.electron_mass_eV, self.target_mass_eV, self.hydrogen_mass_eV,
                self.target_reference_eV, self.hydrogen_reference_eV, self.bohr_m)
        _need(all(type(v) is float and math.isfinite(v) for v in vals), 'finite binary64 registry operands required')
        _need(all(v > 0 for v in (self.electron_mass_eV, self.target_mass_eV,
                                 self.hydrogen_mass_eV, self.bohr_m)), 'nonpositive mass or area conversion')
        _need(self.electron_mass_eV < self.target_mass_eV, 'projectile/target mass ordering')
        _need(self.reference_cost_eV > 0 and math.isfinite(self.reference_cost_eV), 'invalid dissociation reference cost')

    @property
    def reference_cost_eV(self):
        # Same addition/subtraction of the actual binary64 reference operands
        # as the current chemistry table builder. Never subtract GeV masses.
        return (self.hydrogen_reference_eV + self.hydrogen_reference_eV) - self.target_reference_eV


def _load_operands(root):
    root = Path(root) / 'data/tables'
    m = json.loads(_read_pinned(root/'chemistry-v3.2-manifest.json', MANIFEST_SHA256, 2**22))
    _read_pinned(root/'chemistry-v3.2.blob', BLOB_SHA256, 2**24)
    _read_pinned(root/'chemistry-v3.2-spec.json', SPEC_SHA256, 2**22)
    policy = json.loads(_read_pinned(root/'chemistry-v3.2-coverage-policy.json', POLICY_SHA256, 2**22))
    _need(m['blob_sha256'] == BLOB_SHA256 and policy['blob_sha256'] == BLOB_SHA256, 'registry/blob binding')
    _need(policy['admitted_unbounded_regions'] == [], 'strict unknown-region policy required')
    all_species = m['species']
    s = {r['name']: r for r in all_species}
    _need(len(s) == len(all_species), 'duplicate registry species')
    for name, atoms, charge, kind in [('e', {}, -1, 'electron'), ('H', {'H':1}, 0, 'neutral'), ('H2', {'H':2}, 0, 'neutral')]:
        _need(s[name]['atoms'] == atoms and s[name]['charge'] == charge and s[name]['kind'] == kind, 'registry target/product identity')
    _need('e+H2' not in {g['name'] for g in m['groups']}, 'base changed: encounter group is no longer absent')
    return RegistryOperands(s['e']['mass_eV'], s['H2']['mass_eV'], s['H']['mass_eV'],
                            s['H2']['ref_eV'], s['H']['ref_eV'], m['constants']['bohr_m'])


@dataclass(frozen=True)
class SourceEnergy:
    value: Decimal
    unit: str = 'eV'
    coordinate: str = 'published-source-abscissa'
    target: str = 'H2(X1Sg,vi=0)'

    def __post_init__(self):
        _need(type(self.value) is Decimal and self.value.is_finite() and self.value >= 0, 'finite nonnegative Decimal energy required')
        _need(self.unit == 'eV', 'energy units must be eV')
        _need(self.coordinate == 'published-source-abscissa', 'source/target-frame energy mapping is unqualified')
        _need(self.target == 'H2(X1Sg,vi=0)', 'wrong target state or isotope')


@dataclass(frozen=True)
class DiagnosticRate:
    sigma_m2: Decimal
    source_sha256: str
    mode: str = 'diagnostic-source-interpolation'
    interpolation: str = 'piecewise-linear in source eV and sigma'
    event_admitted: bool = False
    network_complete: bool = False


class PrimaryRateComponent:
    __slots__ = ('_data', '_operands')

    def __init__(self, capability, data, operands):
        _need(capability is _CAPABILITY, 'construct only through pinned loader')
        object.__setattr__(self, '_data', data)
        object.__setattr__(self, '_operands', operands)

    def __setattr__(self, name, value):
        raise Refusal('loaded source component is immutable')

    @property
    def operands(self):
        return self._operands

    def diagnostic_rate(self, energy):
        _need(type(energy) is SourceEnergy, 'explicit source-energy operand required')
        # Recheck operand semantics even if a caller defeats frozen dataclass
        # construction externally; this object is not a capability to admit events.
        energy.__post_init__()
        e = energy.value
        rows = self._data.rows
        _need(rows[0][0] <= e <= rows[-1][0], 'outside published table; no endpoint hold or unknown zero')
        with localcontext(_context(50)):
            for i, (x, y) in enumerate(rows):
                if e == x:
                    sigma = y
                    break
                if e < x:
                    a, b = rows[i-1]
                    sigma = b + (y-b) * (e-a) / (x-a)
                    break
            else:
                raise Refusal('internal table search failure')
            a0 = Decimal.from_float(self._operands.bohr_m)
            return DiagnosticRate(sigma * a0 * a0, SOURCE_SHA256)

    def threshold_assessment(self):
        o = self._operands
        with localcontext(_context(70)):
            m, M, Q = (Decimal.from_float(x) for x in (o.electron_mass_eV, o.target_mass_eV, o.reference_cost_eV))
            # From s=(m+M)^2+2MT and final invariant mass m+M+Q.
            threshold = Q * (2*(m+M)+Q) / (2*M)
            source = self._data.threshold
            printed_unit = Decimal('0.00001')
            return {
                'schema': 'h2-threshold-compatibility-v1',
                'source_threshold_eV': str(source),
                'registry_reference_cost_binary64_eV': o.reference_cost_eV,
                'target_rest_incident_threshold_eV': str(threshold),
                'source_threshold_minus_registry_cost_eV': str(source-Q),
                'source_threshold_minus_target_rest_threshold_eV': str(source-threshold),
                'source_printed_threshold_last_digit_eV': str(printed_unit),
                'gap_to_registry_exceeds_one_printed_unit': abs(source-Q) > printed_unit,
                'gap_to_target_rest_threshold_exceeds_one_printed_unit': abs(source-threshold) > printed_unit,
                'printed_resolution_is_not_a_source_uncertainty_bound': True,
                'source_first_positive_table_energy_eV': str(next(e for e,s in self._data.rows if s>0)),
                'source_energy_recoil_convention_verified': False,
                'threshold_alignment_qualified': False,
                'event_kinematics_qualified': False,
                'production_qualified': False,
                'reason': 'unresolved theoretical source onset/reference and energy-frame compatibility; no shift, clamp, or event admission'}

    def event_rate(self, energy):
        _need(type(energy) is SourceEnergy, 'explicit energy operand required')
        energy.__post_init__()
        raise Refusal('event admission refused: threshold/frame, product differential law, and remaining encounter channels unqualified')

    def identity(self):
        return {'schema':'h2-primary-dissociation-rate-v1.1', 'source_sha256':SOURCE_SHA256,
                'blob_sha256':BLOB_SHA256, 'manifest_sha256':MANIFEST_SHA256,
                'policy_sha256':POLICY_SHA256, 'spec_sha256':SPEC_SHA256,
                'target':'H2(X1Sg,vi=0)', 'channel':'b3Su dissociative excitation',
                'interpolation':'piecewise-linear in source eV and sigma',
                'bohr_m_binary64':self._operands.bohr_m,
                'decimal_policy':{'rounding':'ROUND_HALF_EVEN', 'Emin':-999999, 'Emax':999999,
                                  'interpolation_precision':50, 'threshold_precision':70,
                                  'traps':['InvalidOperation','DivisionByZero','Overflow']},
                'network_complete':False, 'production_qualified':False,
                'physical_model_validated':False, 'native_adopted':False,
                'row_count':len(self._data.rows),
                'source_energy_range_eV':[str(self._data.rows[0][0]),str(self._data.rows[-1][0])]}


def load_primary(source_path, chemistry_root):
    data = _parse_source(_read_pinned(source_path, SOURCE_SHA256, 65536))
    _need(len(data.rows) == 54 and data.threshold == Decimal('4.47713') and
          data.rows[-1][0] == Decimal('150'), 'pinned table shape differs')
    return PrimaryRateComponent(_CAPABILITY, data, _load_operands(chemistry_root))
