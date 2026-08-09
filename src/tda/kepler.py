"""Analytic two-body propagation, in universal variables.

The controller needs to know *where* it will be over the interval it is about
to commit to, so that it can probe there.  It cannot integrate the field to
find out --- that would cost the evaluations the probe is trying to price ---
and it does not need to: what the prediction must achieve is a position well
inside :math:`\\pi r/N`, kilometres, not an accurate trajectory.

A two-body propagation delivers that analytically.  Its error is dominated by
the perturbation it omits, :math:`\\tfrac12\\lVert\\mathbf a_{\\mathrm{pert}}
\\rVert\\Delta t_{\\mathrm{dec}}^2`, which at a hundred kilometres over a
two-minute interval is a few metres --- three orders inside the tolerance.

Universal variables rather than a case split on the conic: one formulation
covers ellipse, parabola and hyperbola, and the campaign's populations include
eccentricities where a near-parabolic special case would otherwise have to be
detected and handled [Battin1999]_.  The Stumpff functions are evaluated by
series near zero, where their closed forms cancel catastrophically.

References
----------
.. [Battin1999] R. H. Battin, *An Introduction to the Mathematics and Methods
   of Astrodynamics*, rev. ed., AIAA, 1999, §4.5 -- universal variables and
   the Stumpff functions.
.. [Bate1971] R. R. Bate, D. D. Mueller, J. E. White, *Fundamentals of
   Astrodynamics*, Dover, 1971, §4.4 -- the universal Kepler equation and its
   Newton iteration.
.. [Vallado2013] D. A. Vallado, *Fundamentals of Astrodynamics and
   Applications*, 4th ed., Microcosm, 2013, §2.2 -- ``KEPLER`` and the f and g
   expressions used here.
"""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import NDArray

__all__ = ["KeplerError", "propagate_two_body", "stumpff_c", "stumpff_s"]

Arr = NDArray[np.float64]

_SERIES_CUTOFF = 0.1
"""Below this ``|z|`` the closed forms cancel and the series is used instead."""

_NEWTON_STEPS = 60
_NEWTON_TOLERANCE = 1.0e-12


class KeplerError(RuntimeError):
    """Raised when the universal Kepler equation does not converge.

    Never absorbed into a fallback.  A predictor that silently returned a
    slightly wrong position would move the probe points off the trajectory,
    and the whole argument for forward probing is that the points are on it.
    """


def stumpff_c(z: float) -> float:
    """Stumpff :math:`C(z)=(1-\\cos\\sqrt z)/z`, series near zero.

    Examples
    --------
    >>> round(stumpff_c(0.0), 12)
    0.5
    """
    if z > _SERIES_CUTOFF:
        root = math.sqrt(z)
        return (1.0 - math.cos(root)) / z
    if z < -_SERIES_CUTOFF:
        root = math.sqrt(-z)
        return (math.cosh(root) - 1.0) / (-z)
    # 1/2! - z/4! + z^2/6! - z^3/8! + ...
    return 0.5 - z / 24.0 + z * z / 720.0 - z**3 / 40320.0


def stumpff_s(z: float) -> float:
    """Stumpff :math:`S(z)=(\\sqrt z-\\sin\\sqrt z)/z^{3/2}`, series near zero.

    Examples
    --------
    >>> round(stumpff_s(0.0), 12)
    0.166666666667
    """
    if z > _SERIES_CUTOFF:
        root = math.sqrt(z)
        return (root - math.sin(root)) / root**3
    if z < -_SERIES_CUTOFF:
        root = math.sqrt(-z)
        return (math.sinh(root) - root) / root**3
    # 1/3! - z/5! + z^2/7! - z^3/9! + ...
    return 1.0 / 6.0 - z / 120.0 + z * z / 5040.0 - z**3 / 362880.0


def propagate_two_body(state: Arr, mu: float, dt: float) -> Arr:
    """Propagate a Cartesian state under a point mass, analytically.

    Parameters
    ----------
    state:
        ``[r, v]``, shape ``(6,)``, SI.
    mu:
        Gravitational parameter, m^3 s^-2.
    dt:
        Time step, seconds; may be negative.

    Returns
    -------
    ndarray, shape (6,)

    Raises
    ------
    ValueError
        If the state is misshaped, at the origin, or ``mu`` is not positive.
    KeplerError
        If the Newton iteration does not converge.

    Notes
    -----
    Solves the universal Kepler equation

    .. math::
        \\sqrt\\mu\\,\\Delta t=\\frac{\\mathbf r_0\\cdot\\mathbf v_0}
        {\\sqrt\\mu}\\chi^2C(z)+(1-\\alpha r_0)\\chi^3S(z)+r_0\\chi ,
        \\qquad z=\\alpha\\chi^2 ,

    by Newton's method [Bate1971]_ and then forms the state through the
    :math:`f` and :math:`g` expressions [Vallado2013]_.

    Examples
    --------
    A circular orbit returns to its start after one period.

    >>> import numpy as np
    >>> mu, r = 4.9028e12, 1.8374e6
    >>> s0 = np.array([r, 0, 0, 0, np.sqrt(mu / r), 0])
    >>> period = 2 * np.pi * np.sqrt(r**3 / mu)
    >>> bool(np.allclose(propagate_two_body(s0, mu, period), s0, atol=1e-6))
    True
    """
    state = np.asarray(state, dtype=float)
    if state.shape != (6,):
        raise ValueError(f"state must have shape (6,), got {state.shape}")
    if mu <= 0.0:
        raise ValueError(f"mu must be positive, got {mu}")
    if dt == 0.0:
        return state.copy()

    r0_vec, v0_vec = state[0:3], state[3:6]
    r0 = float(np.linalg.norm(r0_vec))
    if r0 == 0.0:
        raise ValueError("state is at the origin; the two-body flow is "
                         "undefined there")
    root_mu = math.sqrt(mu)
    rv = float(r0_vec @ v0_vec) / root_mu
    alpha = 2.0 / r0 - float(v0_vec @ v0_vec) / mu      # 1 / semi-major axis

    if alpha > 1.0e-12:                                  # ellipse
        chi = root_mu * dt * alpha
    elif alpha < -1.0e-12:                               # hyperbola
        sign = math.copysign(1.0, dt)
        chi = sign * math.sqrt(-1.0 / alpha) * math.log(
            (-2.0 * mu * alpha * dt)
            / (rv * root_mu + sign * math.sqrt(-mu / alpha)
               * (1.0 - r0 * alpha)))
    else:                                                # near-parabolic
        chi = root_mu * dt / r0

    for _ in range(_NEWTON_STEPS):
        z = alpha * chi * chi
        c, s = stumpff_c(z), stumpff_s(z)
        radius = chi * chi * c + rv * chi * (1.0 - z * s) + r0 * (1.0 - z * c)
        residual = (rv * chi * chi * c + (1.0 - alpha * r0) * chi**3 * s
                    + r0 * chi - root_mu * dt)
        if abs(residual) < _NEWTON_TOLERANCE * max(1.0, abs(root_mu * dt)):
            break
        if radius == 0.0:
            raise KeplerError(
                "the Newton iteration hit a zero radius; the state is on a "
                "collision trajectory")
        chi -= residual / radius
    else:
        raise KeplerError(
            f"the universal Kepler equation did not converge in "
            f"{_NEWTON_STEPS} steps for dt={dt:.6g} s at r0={r0:.6g} m")

    z = alpha * chi * chi
    c, s = stumpff_c(z), stumpff_s(z)
    f = 1.0 - chi * chi * c / r0
    g = dt - chi**3 * s / root_mu
    r_vec = f * r0_vec + g * v0_vec
    r_new = float(np.linalg.norm(r_vec))
    if r_new == 0.0:
        raise KeplerError("the propagated state is at the origin")
    f_dot = root_mu * chi * (z * s - 1.0) / (r_new * r0)
    g_dot = 1.0 - chi * chi * c / r_new
    return np.concatenate([r_vec, f_dot * r0_vec + g_dot * v0_vec])
