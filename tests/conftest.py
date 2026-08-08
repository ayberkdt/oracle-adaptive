"""Fixtures and constants shared across the suite.

Kept here rather than duplicated per module so that a change to the test
orbit changes it everywhere at once; a suite in which two files disagree
about what "the reference orbit" means is a suite whose failures are hard to
localise.

The constants are lunar values used to make the fixtures physically sensible.
They are **not** the campaign's adopted values -- those come from the archive
through :mod:`tda.config` and the orbit configurations -- and nothing here
should be imported by production code.
"""

from __future__ import annotations

import numpy as np
import pytest

from tda.config import IntegratorConfig

MU_MOON = 4.902800118e12
"""Lunar gravitational parameter, m^3 s^-2."""

R_MOON = 1.7374e6
"""Lunar mean radius, m."""

J2_MOON = 2.033e-4
"""Unnormalised second zonal coefficient of the Moon."""


@pytest.fixture(scope="session")
def circular_state() -> np.ndarray:
    """A 100 km circular orbit, inclined so that :math:`J_2` acts.

    Inclined deliberately: an equatorial orbit makes the out-of-plane blocks
    of the state-transition matrix decouple, and a decoupled block is one a
    mis-indexed variational term can survive.
    """
    r = R_MOON + 1.0e5
    speed = np.sqrt(MU_MOON / r)
    inc = np.deg2rad(45.0)
    return np.array([r, 0.0, 0.0,
                     0.0, speed * np.cos(inc), speed * np.sin(inc)])


@pytest.fixture(scope="session")
def orbital_period() -> float:
    """Period of :func:`circular_state`, seconds."""
    r = R_MOON + 1.0e5
    return float(2.0 * np.pi * np.sqrt(r**3 / MU_MOON))


@pytest.fixture(scope="session")
def tight() -> IntegratorConfig:
    """Tolerances tighter than the campaign's, so tests measure the method.

    A test that failed because the integrator was loose would be a test of
    the tolerance, not of the code under it.
    """
    return IntegratorConfig(rtol=1.0e-12, atol_position_m=1.0e-9,
                            atol_velocity_m_s=1.0e-12, max_step_s=60.0)
