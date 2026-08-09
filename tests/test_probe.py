"""Verification of the band probe and the analytic predictor.

Three things carry the controller and are checked here rather than assumed:
the two-body predictor really is a two-body flow, the band stack really is
shared -- one evaluation serving every candidate, which is the whole cost
argument -- and the sign of the estimate is right, since getting it backwards
leaves every magnitude untouched and reverses every cancellation decision.

The diagnostics are checked for what they are *not*, too: ``kappa`` and the
retained fraction are pinned to disagree on a case where the transport rotates
the defect, because reading one for the other is the mistake D74 records.
"""

from __future__ import annotations

import numpy as np
import pytest

from conftest import MU_MOON, R_MOON
from tda.analytic import PointMassField
from tda.kepler import KeplerError, propagate_two_body, stumpff_c, stumpff_s
from tda.probe import (
    band_direction,
    candidate_window,
    direction_accuracy,
    plan_probe_points,
    probe_overhead,
    required_degrees,
    retained_fraction,
)
from tda.spectrum import DifferencingBandStack

CANDIDATES = (20, 40, 60, 80, 100)


# ---------------------------------------------------------------------------
# Stumpff functions
# ---------------------------------------------------------------------------


def test_stumpff_series_and_closed_forms_agree_at_the_cutoff() -> None:
    """The series exists because the closed forms cancel near zero.

    Both must be right where they meet, or the predictor has a seam.
    """
    for z in (-0.1000001, -0.0999999, 0.0999999, 0.1000001):
        assert stumpff_c(z) == pytest.approx(stumpff_c(z * 1.0000001),
                                             rel=1e-9)
        assert stumpff_s(z) == pytest.approx(stumpff_s(z * 1.0000001),
                                             rel=1e-9)


def test_stumpff_values_at_zero() -> None:
    assert stumpff_c(0.0) == pytest.approx(0.5)
    assert stumpff_s(0.0) == pytest.approx(1.0 / 6.0)


def test_stumpff_matches_the_closed_forms_away_from_zero() -> None:
    z = 2.0
    assert stumpff_c(z) == pytest.approx((1 - np.cos(np.sqrt(z))) / z)
    assert stumpff_s(z) == pytest.approx(
        (np.sqrt(z) - np.sin(np.sqrt(z))) / z**1.5)
    z = -2.0
    assert stumpff_c(z) == pytest.approx((np.cosh(np.sqrt(-z)) - 1) / (-z))


# ---------------------------------------------------------------------------
# The predictor
# ---------------------------------------------------------------------------


@pytest.fixture
def circular():
    r = R_MOON + 1.0e5
    return np.array([r, 0.0, 0.0, 0.0, np.sqrt(MU_MOON / r), 0.0]), r


def test_a_circular_orbit_closes(circular) -> None:
    state, r = circular
    period = 2.0 * np.pi * np.sqrt(r**3 / MU_MOON)
    assert np.allclose(propagate_two_body(state, MU_MOON, period), state,
                       rtol=0.0, atol=1.0e-6)


def test_the_predictor_agrees_with_a_numerical_two_body_flow(circular) -> None:
    """The claim is analytic exactness, so it is checked against integration.

    A predictor that was merely close would put the probe points off the
    trajectory, and the argument for forward probing is that they are on it.
    """
    from scipy.integrate import solve_ivp

    state, _ = circular
    span = 600.0

    def rhs(_, y):
        r = y[0:3]
        return np.concatenate([y[3:6],
                               -MU_MOON * r / np.linalg.norm(r) ** 3])

    reference = solve_ivp(rhs, (0.0, span), state, rtol=1e-13, atol=1e-9,
                          method="DOP853").y[:, -1]
    assert np.allclose(propagate_two_body(state, MU_MOON, span), reference,
                       rtol=1e-9, atol=1e-6)


def test_energy_is_conserved_over_a_long_step(circular) -> None:
    state, _ = circular
    out = propagate_two_body(state, MU_MOON, 1.0e5)

    def energy(s):
        return 0.5 * s[3:6] @ s[3:6] - MU_MOON / np.linalg.norm(s[0:3])

    assert energy(out) == pytest.approx(energy(state), rel=1e-12)


def test_a_hyperbolic_state_propagates(circular) -> None:
    """Universal variables exist so the conic does not need a case split."""
    state, r = circular
    fast = state.copy()
    fast[4] *= 1.6                              # comfortably unbound
    out = propagate_two_body(fast, MU_MOON, 3600.0)
    assert np.linalg.norm(out[0:3]) > r


def test_forward_then_backward_returns_the_state(circular) -> None:
    state, _ = circular
    there = propagate_two_body(state, MU_MOON, 431.0)
    back = propagate_two_body(there, MU_MOON, -431.0)
    assert np.allclose(back, state, rtol=1e-10, atol=1e-6)


def test_a_zero_step_is_the_identity(circular) -> None:
    state, _ = circular
    assert np.array_equal(propagate_two_body(state, MU_MOON, 0.0), state)


@pytest.mark.parametrize(("state", "mu", "message"), [
    (np.zeros(3), 1.0, r"shape \(6,\)"),
    (np.ones(6), 0.0, "mu must be positive"),
    (np.zeros(6), 1.0, "at the origin"),
])
def test_the_predictor_rejects_degenerate_input(state, mu, message) -> None:
    with pytest.raises(ValueError, match=message):
        propagate_two_body(state, mu, 1.0)


def test_a_radial_state_is_a_degenerate_conic_not_an_error() -> None:
    """Free fall from rest is a valid orbit, and the formulation knows it.

    A zero-angular-momentum ellipse is a genuine conic; universal variables
    continue analytically through the centre rather than failing.  Documented
    as behaviour rather than guarded against, because the campaign's states
    are orbits and a guard here would be dead code that reads as a safeguard.
    """
    radial = np.array([R_MOON, 0.0, 0.0, 0.0, 0.0, 0.0])
    out = propagate_two_body(radial, MU_MOON, 600.0)

    def energy(s):
        return 0.5 * s[3:6] @ s[3:6] - MU_MOON / np.linalg.norm(s[0:3])

    assert np.all(np.isfinite(out))
    assert np.linalg.norm(out[0:3]) < R_MOON          # it has fallen inward
    assert energy(out) == pytest.approx(energy(radial), rel=1e-9)


def test_non_convergence_is_raised_rather_than_returned(monkeypatch) -> None:
    """The guard exists and is reachable.

    Triggered by starving the iteration rather than by finding a pathological
    state, because the formulation is robust enough that a natural one is hard
    to construct -- which is a property of the method, not a reason to leave
    the failure path untested.
    """
    import tda.kepler as kepler

    monkeypatch.setattr(kepler, "_NEWTON_STEPS", 1)
    monkeypatch.setattr(kepler, "_NEWTON_TOLERANCE", 1.0e-300)
    state = np.array([R_MOON + 1.0e5, 0.0, 0.0, 0.0, 1633.0, 0.0])
    with pytest.raises(KeplerError, match="did not converge"):
        kepler.propagate_two_body(state, MU_MOON, 5.0e4)


# ---------------------------------------------------------------------------
# Probe geometry
# ---------------------------------------------------------------------------


def test_probe_points_are_interior_midpoints(circular) -> None:
    """Endpoints would sample a boundary twice and leave the middle unprobed."""
    state, _ = circular
    plan = plan_probe_points(state, MU_MOON, epoch=100.0, span=120.0,
                             n_points=4)
    assert len(plan) == 4
    assert np.allclose(plan.times, 100.0 + np.array([15.0, 45.0, 75.0, 105.0]))
    assert plan.positions.shape == (4, 3)


def test_probe_points_lie_on_the_two_body_arc(circular) -> None:
    state, _ = circular
    plan = plan_probe_points(state, MU_MOON, 0.0, 120.0, 3)
    for t, position in zip(plan.times, plan.positions, strict=True):
        assert np.allclose(propagate_two_body(state, MU_MOON, t)[0:3],
                           position, rtol=0.0, atol=1e-9)


def test_the_predictor_error_bound_is_carried(circular) -> None:
    """Read against ``pi r / N``; a few metres against kilometres."""
    state, _ = circular
    plan = plan_probe_points(state, MU_MOON, 0.0, 120.0, 4,
                             perturbing_acceleration=4.0e-4)
    assert plan.predictor_error_bound == pytest.approx(
        0.5 * 4.0e-4 * 105.0**2)
    assert plan.predictor_error_bound < 0.01 * np.pi * 1.8e6 / 300


@pytest.mark.parametrize(("span", "n_points", "message"), [
    (0.0, 3, "span must be positive"),
    (120.0, 0, "n_points must be positive"),
])
def test_probe_planning_rejects_degenerate_input(circular, span, n_points,
                                                 message) -> None:
    state, _ = circular
    with pytest.raises(ValueError, match=message):
        plan_probe_points(state, MU_MOON, 0.0, span, n_points)


# ---------------------------------------------------------------------------
# The candidate window and the shared stack
# ---------------------------------------------------------------------------


def test_the_window_is_in_indices_not_degrees() -> None:
    """The half-width is an index count, so refining the grid narrows the span.

    That is what keeps the controller's complexity invariant when the
    candidate grid changes.
    """
    coarse = candidate_window((10, 30, 50, 70), 50, 1)
    fine = candidate_window((10, 20, 30, 40, 50, 60, 70), 50, 1)
    assert len(coarse) == len(fine) == 3
    assert coarse == (30, 50, 70)
    assert fine == (40, 50, 60)


def test_the_window_clips_rather_than_extrapolates() -> None:
    assert candidate_window(CANDIDATES, 20, 2) == (20, 40, 60)
    assert candidate_window(CANDIDATES, 100, 2) == (60, 80, 100)


@pytest.mark.parametrize(("planned", "half", "message"), [
    (30, 1, "not a tabulated candidate"),
    (40, -1, "non-negative"),
])
def test_the_window_rejects_bad_input(planned, half, message) -> None:
    with pytest.raises(ValueError, match=message):
        candidate_window(CANDIDATES, planned, half)


def test_required_degrees_is_the_union_of_both_ends() -> None:
    assert list(required_degrees((20, 40), 3)) == [20, 23, 40, 43]


def test_the_stack_is_evaluated_once_for_the_whole_window() -> None:
    """Candidates are partial sums of one stack, not independent probes.

    Whether that stack costs one Legendre pass or one per degree is the
    kernel's business (D120); what this pins is that the *number of degrees
    requested* is the union and not the per-candidate sum.
    """
    stack = DifferencingBandStack(PointMassField(mu=MU_MOON))
    window = (20, 40, 60)
    depth = 3
    sigma = np.ones(200)
    band_direction(stack, np.array([1.8e6, 0.0, 0.0]), 0.0, window, depth,
                   sigma)
    assert stack.total_syntheses == len(required_degrees(window, depth)) == 6
    assert stack.total_syntheses < 2 * len(window) * depth


# ---------------------------------------------------------------------------
# The direction estimate
# ---------------------------------------------------------------------------


class _LinearStack:
    """A stack whose cumulative acceleration is a known function of degree.

    ``a_{<=n} = n * e``, so the band sum over ``k`` degrees is exactly
    ``k * e`` and the estimate can be predicted in closed form.
    """

    def __init__(self, direction):
        self.direction = np.asarray(direction, dtype=float)
        self.total_syntheses = 0

    def cumulative(self, r, t, degrees):
        degrees = np.asarray(degrees)
        self.total_syntheses += degrees.size
        return degrees[:, None] * self.direction[None, :]

    def syntheses_for(self, n_degrees):
        return int(n_degrees)


def test_the_estimate_points_against_the_omitted_tail() -> None:
    """The sign, which no magnitude check would catch.

    ``v_hat`` estimates the *defect* ``a_N - a_ref``, which is minus the
    omitted tail, so the probed bands enter negated.  Getting this backwards
    reverses every cancellation decision while leaving every norm intact.
    """
    stack = _LinearStack([1.0, 0.0, 0.0])
    sigma = np.zeros(200)
    sigma[1:] = 1.0
    out = band_direction(stack, np.zeros(3), 0.0, (20,), 3, sigma)
    assert out[20][0] < 0.0


def test_the_estimate_scales_with_the_completion_factor() -> None:
    """Direction from the bands, amplitude from the spectrum."""
    from tda.spectrum import tail_completion_factor

    stack = _LinearStack([1.0, 0.0, 0.0])
    sigma = np.zeros(200)
    sigma[1:] = 1.0
    depth = 3
    out = band_direction(stack, np.zeros(3), 0.0, (20,), depth, sigma)
    gamma = tail_completion_factor(sigma, 20, depth)
    assert out[20] == pytest.approx(-gamma * depth * np.array([1.0, 0.0, 0.0]))


def test_a_spectrum_too_short_for_the_window_is_refused() -> None:
    stack = _LinearStack([1.0, 0.0, 0.0])
    with pytest.raises(ValueError, match="empty"):
        band_direction(stack, np.zeros(3), 0.0, (20,), 3, np.ones(21))


# ---------------------------------------------------------------------------
# Diagnostics, and what they are not
# ---------------------------------------------------------------------------


def test_direction_accuracy_is_a_cosine() -> None:
    assert direction_accuracy(np.array([1.0, 0, 0]),
                              np.array([2.0, 0, 0])) == pytest.approx(1.0)
    assert direction_accuracy(np.array([1.0, 0, 0]),
                              np.array([0.0, 1, 0])) == pytest.approx(0.0)
    assert direction_accuracy(np.array([1.0, 0, 0]),
                              np.array([-1.0, 0, 0])) == pytest.approx(-1.0)


def test_direction_accuracy_is_undefined_for_a_null_vector() -> None:
    """``nan``, not zero: zero would read as orthogonal."""
    assert np.isnan(direction_accuracy(np.zeros(3), np.array([1.0, 0, 0])))


def test_kappa_and_the_retained_fraction_disagree() -> None:
    """The mistake D74 records, made concrete.

    A perfect ``kappa`` of one forces the retained fraction to one, but the
    converse fails: here the estimate is at forty-five degrees to the truth
    and yet retains the term exactly, because what matters is the angle each
    makes with ``z``, not the angle between them.
    """
    truth = np.array([1.0, 0.0, 0.0])
    estimate = np.array([1.0, 1.0, 0.0])
    z = np.array([1.0, 0.0, 0.0])

    assert direction_accuracy(estimate, truth) == pytest.approx(
        1.0 / np.sqrt(2.0))
    assert retained_fraction(estimate, truth, z) == pytest.approx(1.0)


def test_the_retained_fraction_is_undefined_where_the_term_vanishes() -> None:
    """There the true contribution is zero and any probe error is spurious."""
    assert np.isnan(retained_fraction(np.array([1.0, 0, 0]),
                                      np.array([0.0, 1, 0]),
                                      np.array([1.0, 0, 0])))


# ---------------------------------------------------------------------------
# Cost
# ---------------------------------------------------------------------------


def test_the_overhead_is_proportional_to_the_degree() -> None:
    """The property that makes a single quoted percentage meaningless (D141)."""
    common = {"speed": 1600.0, "radius": 1.8e6, "duration": 6.048e5,
              "rhs_calls": 121_000}
    low = probe_overhead(degree=120, **common)
    high = probe_overhead(degree=600, **common)

    # The ratio is the property; the levels depend on the arc, and this one is
    # a tight circular orbit where the probe is dearer than on an eccentric
    # arc whose long apolune stretches need almost none.
    assert high / low == pytest.approx(5.0, rel=1e-12)
    assert low == pytest.approx(0.1697, rel=1e-3)
    assert high == pytest.approx(0.8485, rel=1e-3)


def test_the_overhead_accepts_a_sampled_arc() -> None:
    """The budget takes the arc integral, not a local rate at perilune.

    Where the correlation time is shortest the vehicle spends the least time,
    so reading the local rate there overstates the cost.
    """
    n = 500
    phase = np.linspace(0.0, 2.0 * np.pi, n)
    radius = 1.8e6 + 1.0e6 * (1.0 - np.cos(phase))
    speed = np.sqrt(MU_MOON * (2.0 / radius - 1.0 / 2.8e6))
    integrated = probe_overhead(300, speed, radius, 6.048e5, 121_000)
    at_perilune = probe_overhead(300, speed[0], radius[0], 6.048e5, 121_000)
    assert integrated < at_perilune


@pytest.mark.parametrize(("kwargs", "message"), [
    ({"rhs_calls": 0}, "rhs_calls must be positive"),
    ({"duration": 0.0}, "duration must be positive"),
])
def test_the_overhead_rejects_degenerate_input(kwargs, message) -> None:
    args = {"degree": 300, "speed": 1600.0, "radius": 1.8e6,
            "duration": 6.048e5, "rhs_calls": 121_000}
    args.update(kwargs)
    with pytest.raises(ValueError, match=message):
        probe_overhead(**args)
