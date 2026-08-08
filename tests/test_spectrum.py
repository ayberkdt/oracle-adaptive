"""The omitted tail: degree variances, amplitude completion, band stacks."""

from __future__ import annotations

import numpy as np
import pytest

from tda.analytic import PointMassField
from tda.spectrum import (
    DifferencingBandStack,
    acceleration_degree_rms,
    coefficient_degree_variance,
    tail_completion_factor,
)

# ---------------------------------------------------------------------------
# Amplitude completion
# ---------------------------------------------------------------------------


def test_gamma_is_one_when_the_probe_covers_the_whole_tail() -> None:
    """If the probed bands are the entire tail, no completion is needed."""
    sigma = np.array([0.0, 1.0, 0.5, 0.25])          # degrees 0..3
    assert tail_completion_factor(sigma, degree=1, depth=2) == pytest.approx(1.0)


def test_gamma_exceeds_one_when_the_tail_continues() -> None:
    sigma = np.array([0.0, 1.0, 0.5, 0.25, 0.125])
    gamma = tail_completion_factor(sigma, degree=1, depth=1)
    expected = np.sqrt((0.5**2 + 0.25**2 + 0.125**2) / 0.5**2)
    assert gamma == pytest.approx(expected)
    assert gamma > 1.0


def test_gamma_uses_the_measured_spectrum_not_a_power_law() -> None:
    """A kinked spectrum must give a different answer from a fitted slope.

    The lunar spectrum is not a single power law -- the previous campaign
    measured a spectral slope of 2.13 against an effective tail exponent of
    1.76 -- so this pins that the tabulated values are what is used.
    """
    n = np.arange(0, 40, dtype=float)
    power_law = np.where(n > 0, np.power(np.maximum(n, 1.0), -1.76), 0.0)
    kinked = power_law.copy()
    kinked[20:] *= 3.0                                # a bump the fit misses

    g_fit = tail_completion_factor(power_law, degree=10, depth=3)
    g_measured = tail_completion_factor(kinked, degree=10, depth=3)
    assert g_measured > g_fit


@pytest.mark.parametrize(("degree", "depth", "message"), [
    (5, 0, "depth must be positive"),
    (5, -1, "depth must be positive"),
    (-1, 3, "degree must be non-negative"),
    (9, 3, "empty"),
])
def test_gamma_rejects_a_degenerate_window(degree, depth, message) -> None:
    """An undefined ratio is an error, never a large number."""
    sigma = np.ones(10)
    with pytest.raises(ValueError, match=message):
        tail_completion_factor(sigma, degree=degree, depth=depth)


def test_gamma_rejects_a_silent_zero_denominator() -> None:
    sigma = np.zeros(10)
    sigma[0] = 1.0
    with pytest.raises(ValueError, match="no power"):
        tail_completion_factor(sigma, degree=2, depth=3)


# ---------------------------------------------------------------------------
# Band stack cost reporting
# ---------------------------------------------------------------------------


def test_band_stack_reports_its_own_cost() -> None:
    """The overhead in ``B_plus`` must be the cost of the code that ran.

    The manuscript's §7.2 model assumes one Legendre pass for the whole
    stack; the kernel has no such entry point yet (``DECISIONS.md`` D120), and
    the differencing implementation must not let that assumption pass
    unnoticed.
    """
    stack = DifferencingBandStack(PointMassField(mu=4.9028e12))
    assert stack.total_syntheses == 0
    assert stack.syntheses_for(5) == 5

    degrees = np.array([10, 11, 12])
    out = stack.cumulative(np.array([1.8e6, 0.0, 0.0]), 0.0, degrees)
    assert out.shape == (3, 3)
    assert stack.total_syntheses == 3

    stack.cumulative(np.array([1.8e6, 0.0, 0.0]), 60.0, degrees)
    assert stack.total_syntheses == 6, "the counter must accumulate"


def test_band_stack_rejects_an_empty_request() -> None:
    """An empty stack sums to zero, which reads as "no omitted acceleration"."""
    stack = DifferencingBandStack(PointMassField(mu=4.9028e12))
    with pytest.raises(ValueError, match="empty"):
        stack.cumulative(np.array([1.8e6, 0.0, 0.0]), 0.0, np.array([], int))


# ---------------------------------------------------------------------------
# Degree variance and per-degree acceleration
# ---------------------------------------------------------------------------


def test_degree_variance_sums_orders_up_to_the_degree() -> None:
    c = np.zeros((3, 3))
    s = np.zeros((3, 3))
    c[2, 0], c[2, 1], s[2, 2] = 3.0, 4.0, 12.0
    assert coefficient_degree_variance(c, s) == pytest.approx(
        [0.0, 0.0, 9.0 + 16.0 + 144.0])


def test_degree_variance_ignores_the_upper_triangle() -> None:
    """A model read from a file is not obliged to zero its padding."""
    c = np.zeros((3, 3))
    c[1, 1] = 2.0
    c[0, 2] = 99.0                      # order 2 at degree 0: not a coefficient
    c[1, 2] = 99.0
    assert coefficient_degree_variance(c, np.zeros((3, 3))) == pytest.approx(
        [0.0, 4.0, 0.0])


@pytest.mark.parametrize(("c", "s"), [
    (np.zeros((2, 3)), np.zeros((2, 3))),
    (np.zeros((3, 3)), np.zeros((2, 2))),
])
def test_degree_variance_rejects_a_bad_shape(c, s) -> None:
    with pytest.raises(ValueError, match="equal square arrays"):
        coefficient_degree_variance(c, s)


def test_acceleration_rms_matches_the_closed_form() -> None:
    """Radial and tangential parts give ``(n+1)(2n+1)`` under the root."""
    mu, r_ref, r = 4.9028e12, 1.7374e6, 1.8374e6
    sigma2 = np.array([0.0, 0.0, 4.0e-8])
    got = acceleration_degree_rms(sigma2, r, r_ref, mu)
    n = 2.0
    expected = ((mu / r**2) * (r_ref / r) ** n
                * np.sqrt((n + 1.0) * (2.0 * n + 1.0)) * np.sqrt(4.0e-8))
    assert got[2] == pytest.approx(expected, rel=1e-12)


def test_acceleration_rms_attenuates_with_altitude() -> None:
    """Upward continuation suppresses high degrees faster than low ones.

    This factor is why the completion ratio is not scale-free and why the
    bare coefficient variances may not be substituted for it.
    """
    sigma2 = np.ones(40)
    low = acceleration_degree_rms(sigma2, 1.8e6, 1.7374e6, 4.9028e12)
    high = acceleration_degree_rms(sigma2, 2.5e6, 1.7374e6, 4.9028e12)
    assert (high[30] / high[2]) < (low[30] / low[2])


def test_acceleration_rms_passes_zero_power_through_exactly() -> None:
    sigma2 = np.array([0.0, 1.0e-8, 0.0, 1.0e-8])
    got = acceleration_degree_rms(sigma2, 1.8e6, 1.7374e6, 4.9028e12)
    assert got[0] == 0.0
    assert got[2] == 0.0
    assert got[1] > 0.0


@pytest.mark.parametrize(("radius", "r_ref", "mu", "message"), [
    (0.0, 1.7e6, 4.9e12, "must be positive"),
    (1.8e6, 0.0, 4.9e12, "must be positive"),
    (1.8e6, 1.7e6, 0.0, "must be positive"),
])
def test_acceleration_rms_rejects_bad_scales(radius, r_ref, mu,
                                             message) -> None:
    with pytest.raises(ValueError, match=message):
        acceleration_degree_rms(np.ones(4), radius, r_ref, mu)


def test_acceleration_rms_rejects_negative_power() -> None:
    with pytest.raises(ValueError, match="non-negative"):
        acceleration_degree_rms(np.array([1.0, -1.0]), 1.8e6, 1.7e6, 4.9e12)


def test_completion_factor_uses_attenuated_power_not_raw_coefficients() -> None:
    """The two give different answers, and the attenuated one is correct.

    A guard against the shortcut of feeding ``coefficient_degree_variance``
    straight into the completion factor: upward continuation weights the far
    tail down, so the raw form overstates gamma.
    """
    sigma2 = np.ones(60)
    sigma_a = acceleration_degree_rms(sigma2, 1.8374e6, 1.7374e6, 4.9028e12)
    attenuated = tail_completion_factor(sigma_a, degree=20, depth=3)
    raw = tail_completion_factor(np.sqrt(sigma2), degree=20, depth=3)
    assert attenuated < raw
