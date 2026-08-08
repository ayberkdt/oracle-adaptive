"""Closed-form gravity fields.

Two roles, and they are the same code:

*Production.*  The ablation ``abl-stm`` replaces the controller's low-degree
state-transition matrix with a Keplerian one, to price whether a closed form
is enough.  That needs a point-mass field with an exact gradient, which is
:class:`PointMassField`.

*Verification.*  :mod:`tda.dynamics` must be testable without the numerical
kernel, both so that the tests run anywhere and so that the state-transition
matrix can be checked against a field whose gradient is known exactly rather
than differenced.  An error in the variational block would otherwise hide
behind an error in the gradient.

Both classes satisfy :class:`tda.field.GravityField`, so they drop into
:func:`tda.dynamics.propagate` unchanged.

References
----------
.. [Battin1999] R. H. Battin, *An Introduction to the Mathematics and Methods
   of Astrodynamics*, rev. ed., AIAA, 1999, §9.3 -- two-body gravity gradient.
.. [Vallado2013] D. A. Vallado, *Fundamentals of Astrodynamics and
   Applications*, 4th ed., Microcosm, 2013, §8.6 -- the :math:`J_2`
   acceleration in Cartesian components.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

__all__ = ["J2Field", "PointMassField"]

Vec3 = NDArray[np.float64]
Mat3 = NDArray[np.float64]


@dataclass(frozen=True, slots=True)
class PointMassField:
    """Central gravity, :math:`\\mathbf a=-\\mu\\mathbf r/r^3`.

    ``degree`` is accepted and ignored, so that the class is substitutable
    wherever a truncated field is expected.  The defect against itself is
    identically zero, which makes it a useful null: an allocation run on this
    field must produce no signal at all, and any that appears is a bug in the
    machinery rather than a property of the Moon.

    Parameters
    ----------
    mu:
        Gravitational parameter, m^3 s^-2.
    reference_radius:
        Carried only to satisfy the protocol; unused.
    """

    mu: float
    reference_radius: float = 1.0

    @property
    def max_degree(self) -> int:
        """Zero: a point mass has no expansion."""
        return 0

    def acceleration(self, r: Vec3, t: float, degree: int = 0) -> Vec3:
        """Central acceleration; ``t`` and ``degree`` are ignored."""
        r = np.asarray(r, dtype=float)
        rn = float(np.linalg.norm(r))
        return -self.mu * r / rn**3

    def defect(self, r: Vec3, t: float, degree: int,
               reference_degree: int) -> Vec3:
        """Identically zero: a point mass has no truncation."""
        return np.zeros(3)

    def exact_gradient(self, r: Vec3) -> Mat3:
        """Closed-form :math:`\\partial\\mathbf a/\\partial\\mathbf r`.

        .. math::
            \\mathbf G = -\\frac{\\mu}{r^{3}}
                         \\bigl(\\mathbf I - 3\\,\\hat{\\mathbf r}
                                \\hat{\\mathbf r}^{\\top}\\bigr),

        the two-body gravity gradient [Battin1999]_.  Symmetric and traceless
        away from the origin, which is what the tests assert.
        """
        r = np.asarray(r, dtype=float)
        rn = float(np.linalg.norm(r))
        rhat = r / rn
        return -(self.mu / rn**3) * (np.eye(3) - 3.0 * np.outer(rhat, rhat))


@dataclass(frozen=True, slots=True)
class J2Field:
    """Central gravity plus the oblateness term.

    Present because a point mass is degenerate for a state-transition matrix
    test: its flow is integrable and its STM has closed form, so an integrator
    can look correct on it for the wrong reason.  Adding :math:`J_2` breaks the
    degeneracy while keeping the gradient available by symbolic
    differentiation of a short expression -- here obtained by differencing the
    *analytic* acceleration, which is exact to roundoff in a way that
    differencing a spherical-harmonic synthesis is not.

    Parameters
    ----------
    mu:
        Gravitational parameter, m^3 s^-2.
    reference_radius:
        Equatorial reference radius, m.
    j2:
        Unnormalised second zonal coefficient.

    Notes
    -----
    Acceleration in the body-equatorial frame follows [Vallado2013]_ §8.6;
    the field is treated as non-rotating, which is what a verification fixture
    needs and is not a model of the Moon.
    """

    mu: float
    reference_radius: float
    j2: float

    @property
    def max_degree(self) -> int:
        """Two: central term plus the second zonal."""
        return 2

    def acceleration(self, r: Vec3, t: float, degree: int = 2) -> Vec3:
        """Central plus :math:`J_2`; ``degree < 2`` drops the oblateness."""
        r = np.asarray(r, dtype=float)
        rn = float(np.linalg.norm(r))
        a = -self.mu * r / rn**3
        if degree < 2:
            return a
        z_r = r[2] / rn
        k = 1.5 * self.j2 * self.mu * self.reference_radius**2 / rn**5
        factor = 5.0 * z_r * z_r
        return a + k * np.array([
            r[0] * (factor - 1.0),
            r[1] * (factor - 1.0),
            r[2] * (factor - 3.0),
        ])

    def defect(self, r: Vec3, t: float, degree: int,
               reference_degree: int) -> Vec3:
        """Difference of two truncations of this field."""
        return (self.acceleration(r, t, degree)
                - self.acceleration(r, t, reference_degree))
