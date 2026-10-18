"""Independent finite-window references for the frozen Temkin source law."""

from dataclasses import FrozenInstanceError
from decimal import Decimal, localcontext
import json
import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from mhdlab.temkin_emission_window import ReducedTemkinEmissionWindow, TemkinWindowRangeError


# Rounded excerpts of the saved 160-digit manufactured, exactly float-projected
# reference. The original rate is a float input there; this API instead takes
# its rounded logarithm, so these tests allow that input projection and the
# declared runtime log-mass qualification. No reference solver runs here.
NONLINEAR_REFERENCE_SHA256 = "96a11fc7457afb93d408196f157e463f7d877350218862d318e5216c09d3b2ed"
CONSTANT_REFERENCE_SHA256 = "ff3748687ef487ea76d820c98fa095d145b8343b564ce224cb832df3a5d569c8"
LOG_MASS_TOLERANCE = 1e-7
REFERENCE_LOG_TOLERANCE = 2e-7
PROBABILITIES = np.array([.1, .5, .9])
REFERENCES = [
    dict(label="moderate_nonlinear", coupling=5., rate=3., start=.1, stop=.7, parent=1.,
         start_depletion=.18738994632687447, local_depletion=.3658456039300304,
         log_local_coupling=1.4220479661072258, log_fraction=-1.370286045087885,
         log_probability=-.6290862116353607,
         log_quantiles=[-3.6496665547007723, -1.6952441489530166, -.7200999631159033]),
    dict(label="maximum_qualified_coupling", coupling=512., rate=1e6, start=.001, stop=.01, parent=.1,
         start_depletion=.025963282838229948, local_depletion=.004616930488210417,
         log_local_coupling=6.212361342201278, log_fraction=-5.406296051102761,
         log_probability=-2.042555242341232,
         log_quantiles=[-8.258974066357244, -6.136596355361112, -4.969981606959333]),
    dict(label="late_endpoint_cdfs_round_one", coupling=4., rate=100., start=100., stop=101., parent=200.,
         start_depletion=185.12367826577307, local_depletion=1.831563888873418,
         log_local_coupling=-183.73738390465317, log_fraction=-185.2982255941575,
         log_probability=-185.2982255941575,
         log_quantiles=[-3.038762646084969, -1.212898198815656, -.2615901323669729]),
    dict(label="narrow_local_window", coupling=3., rate=4., start=1., stop=1.0000000001, parent=2.,
         start_depletion=1.0788805725407125, local_depletion=5.5223705192503486e-11,
         log_local_coupling=.019731716127397304, log_fraction=-24.698509385504423,
         log_probability=-24.45422123982178,
         log_quantiles=[-25.328435940244333, -23.718998027787922, -23.131211362863493]),
    dict(label="late_number_below_float", coupling=4., rate=100., start=1000., stop=1001., parent=2000.,
         start_depletion=1833.5311782518493, local_depletion=1.831563888873418,
         log_local_coupling=-1832.1448838907295, log_fraction=-1833.7057255802338,
         log_probability=-1833.7057255802338,
         log_quantiles=[-3.038762646084969, -1.212898198815656, -.2615901323669729]),
    dict(label="late_fraction_underflow_rescued_by_population", coupling=4., rate=1.,
         start=45000., stop=45001., parent=46000., start_depletion=826.1710393714694,
         local_depletion=.01831563888873418, log_local_coupling=-824.7847450103495,
         log_fraction=-830.1801832133433, log_probability=-830.1801832133433,
         log_quantiles=[-2.310820801557156, -.6977365413054781, -.10628120591272655]),
    dict(label="late_origin_loses_delay_spread", coupling=2., rate=1e-12,
         start=1e12, stop=1000000000000.0003, parent=1000000000001.,
         start_depletion=.5925489801598065, local_depletion=9.984077038474061e-17,
         log_local_coupling=.10059820040013878, log_fraction=-37.43550403326797,
         log_probability=-36.63049499400069,
         log_quantiles=[-10.62035125971339, -9.010913347279288, -8.42312668237717]),
]

REFERENCE_LOG_MASS_DECIMALS = {
    "moderate_nonlinear": "-1.370286045087885164308181415870614352492180316684052015989784673884712106459757122667861813804337342204589066990396804002155575413254545519441743225689803208765",
    "late_endpoint_cdfs_round_one": "-185.2982255941574984626801101680681761680907354208632438749675710765415936883931342264409565772983216330895839284349950630118134552949461916886818509660701710118",
    "late_fraction_underflow_rescued_by_population": "-830.1801832133433312607647036725538083988905292255257441424230886171808572063488102762376253111888405326272673841116707693619550313181382112954683652954199236504",
}


def window(**changes):
    arguments = dict(coupling=5., log_initial_rate_s=math.log(3.), parent_duration_s=1.,
                     start_s=.1, stop_s=.7, log_mass_tolerance=LOG_MASS_TOLERANCE)
    arguments.update(changes)
    return ReducedTemkinEmissionWindow(**arguments)


def from_reference(reference):
    return window(coupling=reference["coupling"], log_initial_rate_s=math.log(reference["rate"]),
                  parent_duration_s=reference["parent"], start_s=reference["start"], stop_s=reference["stop"])


def log_rate_projection_allowance(reference):
    """Bound the effect of supplying log(r) instead of the reference's r.

    For positive original emission flux j, dlog(j)/dlog(r) is
    1-(1+c*exp(-s))*ds/dlog(r), and 0<=ds/dlog(r)<=r*t.
    A positive integral of flux over the window therefore has log derivative
    bounded in magnitude by1+(1+c)*r*b. Extend that bound over the short
    input-projection interval using its largest rate. Decimal precision here
    is independent of the runtime inversion and its reported estimate.
    """
    with localcontext() as context:
        context.prec = 100
        rate = Decimal.from_float(reference["rate"])
        projected_log = Decimal.from_float(math.log(reference["rate"]))
        log_error = abs(projected_log-rate.ln())
        if rate != 1:
            log_error += Decimal("1e-95")
        maximum_rate = rate * log_error.exp()
        sensitivity = 1 + (1 + Decimal.from_float(reference["coupling"])) * maximum_rate * Decimal.from_float(reference["stop"])
        return +(log_error * sensitivity)


@pytest.mark.parametrize("reference", [REFERENCES[0], REFERENCES[2], REFERENCES[5]], ids=lambda item: item["label"])
def test_reference_declared_log_mass_envelope_contains_independent_error_and_input_projection(reference):
    sample = from_reference(reference)
    with localcontext() as context:
        context.prec = 100
        measured_error = abs(Decimal.from_float(sample.log_source_fraction) -
                             Decimal(REFERENCE_LOG_MASS_DECIMALS[reference["label"]]))
        estimated_error = Decimal.from_float(sample.numerical_diagnostics["estimated_log_source_fraction_error"])
        projection = log_rate_projection_allowance(reference)
        # The copied reference has over145 decimal places. This allowance also
        # dominates the finite Decimal arithmetic used in this comparison.
        allowed = estimated_error + projection + Decimal("1e-90")
        assert measured_error <= allowed, (measured_error, estimated_error, projection)
        if reference["rate"] == 1.:
            assert projection == 0


def decimal_ein(coupling):
    """Independent convergent Ein series for the selected0<c<=5 inputs."""
    c = Decimal.from_float(coupling)
    power_over_factorial, result = c, c
    for index in range(2, 400):
        power_over_factorial *= c / Decimal(index)
        contribution = power_over_factorial / Decimal(index)
        result += contribution if index % 2 else -contribution
        if index > coupling and contribution < Decimal("1e-90"):
            return result
    raise AssertionError("independent alternating Ein series did not converge")


@pytest.mark.parametrize("coupling,start", [(1., 271829.), (4., 5459816.), (5., 14841316.)])
def test_large_tail_window_rejects_or_matches_an_independent_ein_log_mass_bound(coupling, start):
    try:
        sample = window(coupling=coupling, log_initial_rate_s=0., parent_duration_s=start+2.,
                        start_s=start, stop_s=start+1., log_mass_tolerance=1e-5)
    except TemkinWindowRangeError:
        return  # Explicit numerical refusal is permitted in this new domain.
    assert sample.start_depletion > 1e5
    with localcontext() as context:
        context.prec = 100
        rate_tail = (-Decimal.from_float(coupling)).exp()
        # I(s,c)=exp(c)[s-Ein(c)+Ein(c*exp(-s))]. At these tails,
        # the final positive correction is smaller than1e-1000.
        expected_start = Decimal.from_float(start)*rate_tail + decimal_ein(coupling)
        expected_log_mass = -expected_start + (1-(-rate_tail).exp()).ln()
        measured_error = abs(Decimal.from_float(sample.log_source_fraction)-expected_log_mass)
        estimated_error = Decimal.from_float(sample.numerical_diagnostics["estimated_log_source_fraction_error"])
        assert measured_error <= estimated_error + Decimal("1e-80"), (measured_error, estimated_error)
        assert estimated_error <= Decimal.from_float(1e-5)


@pytest.mark.parametrize("reference", REFERENCES, ids=lambda item: item["label"])
def test_saved_nonlinear_log_amounts_and_local_quantiles(reference):
    sample = from_reference(reference)
    assert sample.start_depletion == pytest.approx(reference["start_depletion"], rel=0., abs=REFERENCE_LOG_TOLERANCE)
    assert sample.local_depletion == pytest.approx(reference["local_depletion"], rel=3e-7, abs=0.)
    assert sample.log_local_coupling == pytest.approx(reference["log_local_coupling"], rel=0., abs=REFERENCE_LOG_TOLERANCE)
    assert sample.log_source_fraction == pytest.approx(reference["log_fraction"], rel=0., abs=REFERENCE_LOG_TOLERANCE)
    assert sample.log_parent_conditional_probability == pytest.approx(reference["log_probability"], rel=0., abs=REFERENCE_LOG_TOLERANCE)
    np.testing.assert_allclose(sample.log_quantile(PROBABILITIES), reference["log_quantiles"],
                               rtol=0., atol=REFERENCE_LOG_TOLERANCE)
    delays = sample.quantile(PROBABILITIES)
    assert np.all(delays > 0) and np.all(delays < sample.local_duration_s)
    np.testing.assert_allclose(sample.cdf(delays), PROBABILITIES, rtol=3e-7, atol=0.)
    assert sample.cdf(0.) == 0.
    assert sample.cdf(sample.local_duration_s) == 1.


def test_late_survival_is_uncapped_and_its_positive_log_amount_is_retained():
    sample = from_reference(REFERENCES[4])
    assert sample.start_depletion > 1800.
    assert math.isfinite(sample.log_source_fraction)
    assert sample.log_source_fraction < -1800.
    assert math.exp(sample.log_source_fraction) == 0.
    # The saved reference's original endpoint CDFs both round to one.
    # Their zero difference cannot replace this finite logarithmic source mass.
    assert sample.log_parent_conditional_probability < -1800.
    assert np.all(np.isfinite(sample.log_quantile(PROBABILITIES)))


def test_underflowed_fraction_can_still_rescue_a_representable_initial_population_count():
    sample = from_reference(REFERENCES[5])
    assert math.exp(sample.log_source_fraction) == 0.
    log_count = math.log(1e300) + sample.log_source_fraction
    assert log_count == pytest.approx(-139.40465531512962, rel=0., abs=REFERENCE_LOG_TOLERANCE)
    assert 0. < math.exp(log_count) < 1.
    assert math.isfinite(sample.log_local_coupling) and sample.log_local_coupling < -800.


def test_delay_coordinates_preserve_spread_lost_in_log_absolute_times():
    sample = from_reference(REFERENCES[6])
    delays = sample.quantile(PROBABILITIES)
    assert np.all(np.diff(delays) > 0)
    assert np.all(np.diff(sample.log_quantile(PROBABILITIES)) > 0)
    assert len(set(np.log(sample.start_s + delays).tolist())) == 1
    with pytest.raises(ValueError):
        sample.cdf(sample.start_s + delays)


def test_zero_coupling_matches_independent_truncated_exponential_formulas():
    rate, start, stop, parent = 1.7, 2., 2.4, 3.
    sample = window(coupling=0., log_initial_rate_s=math.log(rate), parent_duration_s=parent,
                    start_s=start, stop_s=stop)
    duration = stop - start
    local_release = -math.expm1(-rate * duration)
    expected_log_mass = -rate * start + math.log(local_release)
    expected_log_probability = expected_log_mass - math.log(-math.expm1(-rate * parent))
    assert sample.log_source_fraction == pytest.approx(expected_log_mass, rel=0., abs=2e-12)
    assert sample.log_parent_conditional_probability == pytest.approx(expected_log_probability, rel=0., abs=2e-12)
    times = np.array([0., .03, .1, .2, duration])
    expected_cdf = -np.expm1(-rate * times) / local_release
    np.testing.assert_allclose(sample.cdf(times), expected_cdf, rtol=3e-12, atol=0.)
    expected_quantiles = -np.log1p(-PROBABILITIES * local_release) / rate
    np.testing.assert_allclose(sample.quantile(PROBABILITIES), expected_quantiles, rtol=3e-12, atol=0.)


def test_local_cdf_at_new_delay_points_matches_a_direct_original_time_ode():
    # New moderate manufactured inputs, not a rerun of the saved oracle cases.
    coupling, rate, start, stop = 2.3, 7.2, .04, .8
    sample = window(coupling=coupling, log_initial_rate_s=math.log(rate), parent_duration_s=1.,
                    start_s=start, stop_s=stop)

    def rhs(_time, state):
        return rate * np.exp(coupling * np.expm1(-state))

    solution = solve_ivp(rhs, (0., stop), [0.], method="DOP853", dense_output=True,
                         rtol=2e-12, atol=1e-14)
    assert solution.success
    start_depletion, stop_depletion = (float(solution.sol(time)[0]) for time in (start, stop))
    local_release = -math.expm1(-(stop_depletion-start_depletion))
    delays = (stop-start) * np.array([.01, .1, .3, .7, .95])
    local_depletion = solution.sol(start+delays)[0] - start_depletion
    expected_cdf = -np.expm1(-local_depletion) / local_release
    np.testing.assert_allclose(sample.cdf(delays), expected_cdf, rtol=3e-10, atol=0.)
    expected_log_mass = -start_depletion + math.log(local_release)
    assert sample.log_source_fraction == pytest.approx(expected_log_mass, rel=0., abs=3e-10)


def test_tiny_hazard_has_uniform_local_time_without_erasing_log_amount_or_parent_normalizer():
    sample = window(coupling=0., log_initial_rate_s=math.log(1e-200), parent_duration_s=1e-100,
                    start_s=0., stop_s=1e-200)
    assert sample.log_source_fraction == pytest.approx(math.log(1e-200) * 2, rel=0., abs=2e-12)
    assert sample.log_parent_conditional_probability == pytest.approx(math.log(1e-100), rel=0., abs=2e-12)
    np.testing.assert_allclose(sample.quantile(PROBABILITIES), PROBABILITIES * 1e-200,
                               rtol=3e-12, atol=0.)
    np.testing.assert_allclose(sample.cdf(PROBABILITIES * 1e-200), PROBABILITIES, rtol=3e-12, atol=0.)


def test_log_cdf_retains_a_positive_probability_that_underflows_ordinary_cdf():
    log_rate, duration, delay = math.log(1e-200), 1e200, 1e-300
    sample = window(coupling=0., log_initial_rate_s=log_rate, parent_duration_s=duration,
                    start_s=0., stop_s=duration)
    local_depletion = math.exp(math.fsum((log_rate, math.log(duration))))
    expected = math.fsum((log_rate, math.log(delay), -math.log(-math.expm1(-local_depletion))))
    assert sample.log_cdf(delay) == pytest.approx(expected, rel=0., abs=3e-10)
    assert math.isfinite(sample.log_cdf(delay))
    assert sample.cdf(delay) == 0.
    assert sample.log_cdf(0.) == -math.inf
    assert sample.log_cdf(duration) == 0.


@pytest.mark.parametrize("reference", [REFERENCES[0], REFERENCES[2], REFERENCES[5]], ids=lambda item: item["label"])
def test_partition_conservation_is_checked_inside_the_selected_window_in_log_space(reference):
    whole = from_reference(reference)
    midpoint = reference["start"] + .5 * (reference["stop"] - reference["start"])
    common = dict(coupling=reference["coupling"], log_initial_rate_s=math.log(reference["rate"]),
                  parent_duration_s=reference["parent"])
    first = window(**common, start_s=reference["start"], stop_s=midpoint)
    second = window(**common, start_s=midpoint, stop_s=reference["stop"])
    log_sum = np.logaddexp(first.log_source_fraction, second.log_source_fraction)
    assert log_sum == pytest.approx(whole.log_source_fraction, rel=0., abs=3e-7)


def test_cdf_quantile_replay_at_noncentral_probabilities():
    sample = window()
    probabilities = np.array([1e-8, 1e-4, .1, .5, .9, 1.-1e-4, 1.-1e-8])
    replay = sample.cdf(sample.quantile(probabilities))
    np.testing.assert_allclose(replay, probabilities, rtol=5e-7, atol=0.)
    np.testing.assert_allclose(1.-replay, 1.-probabilities, rtol=5e-4, atol=0.)


def test_log_delays_survive_when_requested_linear_projection_underflows():
    sample = window(coupling=0., log_initial_rate_s=math.log(1e300), parent_duration_s=1e-299,
                    start_s=0., stop_s=1e-300)
    logarithm = sample.log_quantile(1e-30)
    assert math.isfinite(logarithm) and logarithm < math.log(np.nextafter(0., 1.))
    with pytest.raises(TemkinWindowRangeError):
        sample.quantile(1e-30)


def test_unresolved_subnormal_delay_projection_cannot_claim_high_relative_accuracy():
    sample = window(coupling=0., log_initial_rate_s=math.log(1e308), parent_duration_s=1e-319,
                    start_s=0., stop_s=1e-320)
    # The median happens to land on the subnormal grid; probability.1 does not.
    assert math.isfinite(sample.log_quantile(.1))
    with pytest.raises(TemkinWindowRangeError):
        sample.quantile(.1)


@pytest.mark.parametrize("changes", [
    {"coupling": -.1}, {"coupling": 513.}, {"start_s": -.1}, {"start_s": .7},
    {"stop_s": .1}, {"stop_s": 1.1}, {"parent_duration_s": 0.},
    {"log_mass_tolerance": 0.}, {"stop_s": float("1.0000000000000001"), "start_s": 1., "parent_duration_s": 2.},
])
def test_invalid_or_projected_empty_windows_reject(changes):
    with pytest.raises(ValueError):
        window(**changes)


@pytest.mark.parametrize("name", ["coupling", "log_initial_rate_s", "parent_duration_s", "start_s", "stop_s", "log_mass_tolerance"])
@pytest.mark.parametrize("bad", [True, "1", np.nan, np.inf, 1.+0j, np.ma.array(1., mask=False)])
def test_constructor_rejects_coercions_nonfinite_and_masked_inputs(name, bad):
    with pytest.raises(ValueError):
        window(**{name: bad})


@pytest.mark.parametrize("method,bad", [
    ("cdf", -.1), ("cdf", .7), ("cdf", [0., True]), ("cdf", [".1"]),
    ("cdf", np.array([[.1]])), ("cdf", np.ma.array([.1], mask=False)),
    ("quantile", 0.), ("quantile", 1.), ("quantile", -.1), ("quantile", [True, .5]),
    ("quantile", [np.nan]), ("quantile", [np.inf]), ("quantile", [.5+0j]),
    ("quantile", np.array([[.5]])), ("quantile", np.ma.array([.5], mask=False)),
    ("log_quantile", 0.), ("log_quantile", 1.), ("log_quantile", "0.5"),
])
def test_local_time_and_probability_inputs_have_strict_supported_domains(method, bad):
    with pytest.raises(ValueError):
        getattr(window(), method)(bad)


def test_requested_mass_precision_must_be_numerically_qualified():
    with pytest.raises(TemkinWindowRangeError):
        window(log_mass_tolerance=1e-30)
    with pytest.raises(TemkinWindowRangeError):
        window(coupling=0., log_initial_rate_s=0., parent_duration_s=1e12+2.,
               start_s=1e12, stop_s=1e12+1., log_mass_tolerance=1e-8)


def test_qualification_diagnostics_do_not_substitute_time_residual_for_log_mass_error():
    sample = from_reference(REFERENCES[4])
    diagnostics = sample.numerical_diagnostics
    assert diagnostics["estimated_log_source_fraction_error"] <= LOG_MASS_TOLERANCE
    assert diagnostics["estimated_log_parent_probability_error"] <= LOG_MASS_TOLERANCE
    assert diagnostics["start_depletion_error"] > 0.
    assert diagnostics["local_log_depletion_error"] > 0.
    assert diagnostics["local_coupling_approximated_as_zero"]
    assert math.isfinite(diagnostics["log_local_coupling_approximation_error_bound"])
    assert diagnostics["quadrature_assumption"]


@pytest.mark.parametrize("changes", [
    {"coupling": 0.},
    {"coupling": 4., "log_initial_rate_s": math.log(100.), "parent_duration_s": 2000.,
     "start_s": 1000., "stop_s": 1001.},
    {"coupling": 0., "log_initial_rate_s": math.log(1e-200), "parent_duration_s": 1e-100,
     "start_s": 0., "stop_s": 1e-200},
])
def test_manifest_is_strict_json_even_for_exact_zero_or_underflowed_derived_quantities(changes):
    sample = window(**changes)
    manifest = sample.manifest()
    assert math.isfinite(sample.log_source_fraction)
    json.dumps(manifest, allow_nan=False)


def test_window_and_returned_arrays_are_immutable_and_inputs_are_not_changed():
    sample = window()
    probabilities = np.array([.1, .5, .9])
    before = probabilities.copy()
    delays, logarithms = sample.quantile(probabilities), sample.log_quantile(probabilities)
    cdf = sample.cdf(delays)
    np.testing.assert_array_equal(probabilities, before)
    for values in (delays, logarithms, cdf):
        with pytest.raises(ValueError):
            values.setflags(write=True)
    with pytest.raises((FrozenInstanceError, AttributeError)):
        sample.start_s = .2
    manifest = sample.manifest()
    json.dumps(manifest, allow_nan=False)
    manifest.clear()
    assert sample.manifest()
