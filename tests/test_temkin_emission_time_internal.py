"""Representation contracts for frozen Temkin conditional event times.

Parameters are manufactured mathematical cases, with no species/material fit.
"""

import math

import numpy as np
import pytest

from mhdlab.temkin_emission_time import (
    FrozenTemkinEmission, TemkinEmissionTimeRangeError,
)
from mhdlab.temkin_release import TemkinReservoir


def pool(**overrides):
    values = dict(name="manufactured", state_id="M:ground", attempt_frequency_s=1.,
                  activation_energy_ev=0., desorption_energy_ev=0., internal_energy_ev=0.,
                  coverage_coefficient=.5, reference_coverage_m2=1., max_coverage_m2=1.,
                  kinetic_provenance="manufactured mathematical rate",
                  energy_provenance="manufactured independent zero binding/internal cost")
    values.update(overrides)
    return TemkinReservoir(**values)


def test_conditional_fraction_does_not_divide_rounded_subnormal_inventory():
    smallest = float(np.nextafter(0., 1.))
    source = FrozenTemkinEmission(pool(reference_coverage_m2=smallest, max_coverage_m2=smallest),
                                 smallest, 1., 1.)
    # Both extensive arrays round as if all inventory left, despite an exact
    # first-order released fraction 1-exp(-1). Timing must retain that fraction.
    assert source.released_m2 == smallest and source.remaining_coverage_m2 == 0
    expected = -math.log1p(-.5 * (-math.expm1(-1.)))
    assert source.log_released_fraction == pytest.approx(math.log(-math.expm1(-1.)), rel=2e-13)
    assert source.quantile(.5) == pytest.approx(expected, rel=2e-13)
    assert source.cdf(expected) == pytest.approx(.5, rel=2e-13)


def test_positive_but_coarsely_quantized_offset_requires_log_interface():
    source = FrozenTemkinEmission(pool(reference_coverage_m2=1e300, max_coverage_m2=1e300),
                                 1e300, 1., 1e-320)
    probability = .123456789
    expected_log_time = math.fsum((math.log(source.dt_s), math.log(probability)))
    assert source.log_quantile(probability) == pytest.approx(expected_log_time, abs=1e-12)
    with pytest.raises(TemkinEmissionTimeRangeError, match="coarsely resolved"):
        source.quantile(probability)
    # A represented subnormal offset still has a well-defined law CDF.
    offset = float(np.nextafter(0., 1.))
    assert source.cdf(offset) == pytest.approx(offset / source.dt_s, rel=2e-13)


def test_survivor_complement_is_kept_for_probabilities_next_to_one():
    source = FrozenTemkinEmission(pool(), 1., 1., 40.)
    draw = float(np.nextafter(1., 0.))
    expected = -math.log(math.fsum((1 - draw, draw * math.exp(-40.))))
    # 1-exp(-40) rounds to one.  Retaining the survivor avoids incorrectly
    # using -log(1-draw), whose difference is meaningful at this endpoint.
    assert -math.expm1(-40.) == 1.
    assert source.quantile(draw) == pytest.approx(expected, rel=2e-13)
    assert abs(source.quantile(draw) + math.log1p(-draw)) > .01
