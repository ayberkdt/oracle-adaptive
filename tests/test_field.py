"""Field-side arithmetic that does not need the synthesis kernel.

Only the inertial/body transform is exercised here.  The kernel adapter
itself is covered by the admissibility check of WP0, which compares against
archived numbers and needs the archive on the path; that is an integration
test and does not belong in this file.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from tda.config import OMEGA_MOON_RAD_S
from tda.field import _RotationZ

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
