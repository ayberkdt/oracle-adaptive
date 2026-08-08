"""Configuration immutability and provenance.

The digest is what ties a reported number to the settings that produced it,
so the tests here are about the hash covering everything that can move a
result -- not about the values of the defaults, which are decisions recorded
in ``DECISIONS.md`` and not properties of the code.
"""

from __future__ import annotations

import dataclasses
import math

import pytest

from tda.config import (
    OMEGA_MOON_RAD_S,
    GradientConfig,
    IntegratorConfig,
    RunConfig,
)


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
# Body rotation rate
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
