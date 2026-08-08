"""Configuration provenance and the field-side arithmetic.

Nothing here touches the spherical-harmonic kernel.  What is tested is the
part of :mod:`tda.field` that is arithmetic rather than synthesis -- the
amplitude completion factor and the input guards -- together with the
immutability and hashing that the manifest chain depends on.
"""

from __future__ import annotations

import dataclasses
import math

import numpy as np
import pytest

from tda.analytic import PointMassField
from tda.config import (
    OMEGA_MOON_RAD_S,
    GradientConfig,
    IntegratorConfig,
    RunConfig,
)
from tda.field import (
    DifferencingBandStack,
    _RotationZ,
    degree_variance_gamma,
)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------


def test_config_is_frozen() -> None:
    """A setting must not change after a run has started."""
    cfg = RunConfig(label="x")
    with pytest.raises(dataclasses.FrozenInstanceError):
        cfg.integrator.rtol = 1.0e-9          # type: ignore[misc]
    with pytest.raises(dataclasses.FrozenInstanceError):
        cfg.label = "y"                        # type: ignore[misc]


def test_digest_is_stable_and_order_independent() -> None:
    """Two equal configurations must hash alike, whatever the build order."""
    a = RunConfig(label="run", gradient=GradientConfig(degree=300),
                  integrator=IntegratorConfig(rtol=1e-11))
    b = RunConfig(label="run", integrator=IntegratorConfig(rtol=1e-11),
                  gradient=GradientConfig(degree=300))
    assert a.digest() == b.digest()
    assert len(a.digest()) == 16


def test_digest_separates_settings_that_change_a_number() -> None:
    """The gradient degree is the open question Q13; it must be in the hash."""
    base = RunConfig(label="run")
    changed = RunConfig(label="run", gradient=GradientConfig(degree=120))
    assert base.digest() != changed.digest()


def test_absolute_tolerance_is_a_vector_not_a_scalar() -> None:
    """Position and velocity differ by orders of magnitude in SI.

    The previous campaign traced a noise-floor artefact to a scalar ``atol``
    that was far too loose on the velocity components; the vector form is the
    fix, and this pins it.
    """
    cfg = IntegratorConfig()
    atol = cfg.atol_vector(n_extra=36)
    assert len(atol) == 42
    assert atol[0] == cfg.atol_position_m
    assert atol[3] == cfg.atol_velocity_m_s
    assert atol[3] < atol[0]


# ---------------------------------------------------------------------------
# Amplitude completion
# ---------------------------------------------------------------------------


def test_gamma_is_one_when_the_probe_covers_the_whole_tail() -> None:
    """If the probed bands are the entire tail, no completion is needed."""
    sigma = np.array([0.0, 1.0, 0.5, 0.25])          # degrees 0..3
    assert degree_variance_gamma(sigma, degree=1, depth=2) == pytest.approx(1.0)


def test_gamma_exceeds_one_when_the_tail_continues() -> None:
    sigma = np.array([0.0, 1.0, 0.5, 0.25, 0.125])
    gamma = degree_variance_gamma(sigma, degree=1, depth=1)
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

    g_fit = degree_variance_gamma(power_law, degree=10, depth=3)
    g_measured = degree_variance_gamma(kinked, degree=10, depth=3)
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
        degree_variance_gamma(sigma, degree=degree, depth=depth)


def test_gamma_rejects_a_silent_zero_denominator() -> None:
    sigma = np.zeros(10)
    sigma[0] = 1.0
    with pytest.raises(ValueError, match="no power"):
        degree_variance_gamma(sigma, degree=2, depth=3)


# ---------------------------------------------------------------------------
# Body rotation -- must match the archive exactly
# ---------------------------------------------------------------------------


def test_rotation_rate_is_the_archive_expression() -> None:
    """Pins that the rate is computed, not copied as a rounded literal.

    ``rev3_common.OMEGA_MOON = 2*pi/(27.321661*86400)``.  A literal truncated
    to seven digits would drift the body-fixed longitude by metres over a
    seven-day arc -- invisible to a magnitude statistic and fatal to a signed
    one.
    """
    archive_expression = 2.0 * math.pi / (27.321661 * 86400.0)
    assert archive_expression == OMEGA_MOON_RAD_S


def test_body_transform_matches_the_archive_formula() -> None:
    """The inertial/body convention is reproduced verbatim, sign for sign.

    Compared exactly, not approximately: the implementation and the archive
    perform the same IEEE operations on the same values, so any difference
    would be a changed formula rather than accumulated roundoff.  ``math`` is
    used rather than ``numpy`` for the trigonometry because the two need not
    agree in the last unit in the last place.
    """
    rot = _RotationZ(OMEGA_MOON_RAD_S)
    r = np.array([1.8e6, -4.0e5, 9.0e4])
    t = 12345.6

    theta = OMEGA_MOON_RAD_S * t
    c, s = math.cos(theta), math.sin(theta)
    expected_body = (c * r[0] + s * r[1], -s * r[0] + c * r[1], r[2])
    assert rot.to_body(r, t) == pytest.approx(expected_body, rel=0, abs=0)

    a_body = (1.0, 2.0, 3.0)
    expected_inertial = np.array([c * 1.0 - s * 2.0, s * 1.0 + c * 2.0, 3.0])
    assert rot.to_inertial(a_body, t) == pytest.approx(expected_inertial,
                                                       rel=0, abs=0)


def test_body_transform_round_trips() -> None:
    """A vector taken to the body frame and back is unchanged."""
    rot = _RotationZ(OMEGA_MOON_RAD_S)
    r = np.array([1.8e6, -4.0e5, 9.0e4])
    t = 98765.4
    assert rot.to_inertial(rot.to_body(r, t), t) == pytest.approx(r, rel=1e-15)


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
