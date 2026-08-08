"""Gravity-field adapter.

Evaluating the field that is *retained*: the acceleration itself, the
*truncation defect* that is the allocation's input, and the gradient that
generates the state-transition matrix.

What a truncation leaves *out* -- the omitted degree bands and the spectrum
that restores their amplitude -- lives in :mod:`tda.spectrum`, which depends
on this module and not the other way round.

Nothing here reimplements spherical-harmonic synthesis.  The kernel is the
archive's, a Pines-type singularity-free formulation [Pines1973]_ with the
fully normalised recursions of Holmes and Featherstone [Holmes2002]_, and it
is called verbatim so that the admissibility check of WP0 can compare against
archived numbers.  What this module adds is (i) the inertial/body transform in
the exact form the archive uses, (ii) a one-pass defect that avoids paying for
two independent syntheses, and (iii) an explicit, recorded choice of gradient
degree, which the plan flags as an open question rather than a detail.

Conventions
-----------
* ``r`` is an inertial Cartesian position in metres, shape ``(3,)``.
* ``degree`` is the truncation degree :math:`N`; the model's own maximum is
  :attr:`GravityField.max_degree`.
* The body-fixed frame is reached by a uniform rotation about the inertial
  z-axis, reproducing the archive's simplification rather than improving on
  it (see :data:`tda.config.OMEGA_MOON_RAD_S`).

References
----------
.. [Pines1973] S. Pines, "Uniform representation of the gravitational
   potential and its derivatives", *AIAA Journal* 11(11), 1973.
.. [Holmes2002] S. A. Holmes and W. E. Featherstone, "A unified approach to
   the Clenshaw summation and the recursive computation of very high degree
   and order normalised associated Legendre functions", *Journal of Geodesy*
   76, 2002.
.. [Press2007] W. H. Press et al., *Numerical Recipes*, 3rd ed., CUP, 2007,
   §5.7 -- central differences and step-size choice.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import numpy as np
from numpy.typing import NDArray

from tda.config import GradientConfig

__all__ = [
    "GravityField",
    "LunarisField",
    "gravity_gradient",
]

Vec3 = NDArray[np.float64]
Mat3 = NDArray[np.float64]


# ---------------------------------------------------------------------------
# Interface
# ---------------------------------------------------------------------------


@runtime_checkable
class GravityField(Protocol):
    """Minimal contract the rest of the package depends on.

    Keeping this a :class:`~typing.Protocol` rather than a base class means
    the tests can substitute an analytic field -- a point mass, or a
    :math:`J_2` field with a known closed-form gradient -- without importing
    the numerical kernel at all.  That is what makes
    :mod:`tda.dynamics` testable.
    """

    @property
    def max_degree(self) -> int:
        """Highest degree the loaded coefficient set supports."""

    @property
    def mu(self) -> float:
        """Gravitational parameter, m^3 s^-2."""

    @property
    def reference_radius(self) -> float:
        """Reference radius of the expansion, m."""

    def acceleration(self, r: Vec3, t: float, degree: int) -> Vec3:
        """Inertial acceleration, truncated at ``degree``.

        Parameters
        ----------
        r:
            Inertial position, shape ``(3,)``, m.
        t:
            Epoch, seconds from arc start; enters through the body rotation.
        degree:
            Truncation degree :math:`N`.

        Returns
        -------
        ndarray, shape (3,)
            Acceleration, m s^-2.
        """

    def defect(self, r: Vec3, t: float, degree: int,
               reference_degree: int) -> Vec3:
        """Truncation defect at one point.

        :math:`\\Delta\\mathbf a=\\mathbf a_N-\\mathbf a_{\\mathrm{ref}}`,
        the input to the whole allocation.

        Parameters
        ----------
        r:
            Inertial position, shape ``(3,)``, m.
        t:
            Epoch, seconds from arc start.
        degree:
            Truncated (policy) degree :math:`N`.
        reference_degree:
            The orbit's adopted reference degree.

        Returns
        -------
        ndarray, shape (3,)
            Defect, m s^-2.
        """


# ---------------------------------------------------------------------------
# Archive-backed implementation
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class _RotationZ:
    """Uniform rotation about the inertial z-axis.

    Written out rather than delegated so that the inertial/body convention is
    visible at the point of use and can be diffed against the archive's
    ``rev13_variational_check.accel_inertial``.  The two must agree exactly;
    a sign flip here would be invisible in magnitudes and fatal in a signed
    method.
    """

    omega: float

    def to_body(self, r: Vec3, t: float) -> tuple[float, float, float]:
        """Rotate an inertial position into the body-fixed frame."""
        c, s = math.cos(self.omega * t), math.sin(self.omega * t)
        x, y, z = float(r[0]), float(r[1]), float(r[2])
        return c * x + s * y, -s * x + c * y, z

    def to_inertial(self, a_body: tuple[float, float, float],
                    t: float) -> Vec3:
        """Rotate a body-fixed vector back into the inertial frame."""
        c, s = math.cos(self.omega * t), math.sin(self.omega * t)
        ax, ay, az = a_body
        return np.array([c * ax - s * ay, s * ax + c * ay, az])


class LunarisField:
    """Adapter over the archive's spherical-harmonic kernel.

    Parameters
    ----------
    model:
        The loaded gravity model object of the archive.
    kernel_args:
        The pre-unpacked argument tuple the numba kernel expects.  Held as a
        tuple rather than re-derived per call because unpacking dominates the
        cost of a low-degree evaluation.
    omega:
        Body rotation rate, s^-1.  Read from the archive, not assumed.

    Notes
    -----
    Construct through :meth:`from_archive` in production code; the explicit
    constructor exists so that a test can inject a stub kernel.

    The instance is **not** thread-safe: the underlying kernel writes into
    pre-allocated Legendre workspaces held by ``kernel_args``.  Parallelism in
    this campaign is by process, matching the archive.
    """

    __slots__ = ("_args", "_dual", "_model", "_rot")

    def __init__(self, model, kernel_args: tuple, omega: float) -> None:
        """Wrap a loaded model; see the class docstring for the parameters."""
        self._model = model
        self._args = kernel_args
        self._rot = _RotationZ(omega)
        self._dual = _load_dual_kernel()

    # -- construction -------------------------------------------------------

    @classmethod
    def from_archive(cls, degree: int, omega: float) -> LunarisField:
        """Load the archive's model at ``degree`` and warm the JIT.

        The loaders are imported from ``rev3_common``, which is where they are
        defined; ``rev10_sobol_confirmatory`` re-exports the same function
        objects and is what the archive's own variational scripts import, but
        going to the definition avoids pulling in a large module for three
        names.

        The warm-up call is not optional.  Numba compiles on first call, and
        an uncompiled first evaluation inside a timed section would corrupt
        the work measurements the budget accounting rests on.

        Requires the archive's ``python_codes`` directory on ``sys.path``;
        :mod:`tda.archive` is what puts it there in production code.
        """
        from rev3_common import kernel_args, load_model, warmup

        model = load_model(degree)
        args = kernel_args(model)
        warmup(model, args)
        return cls(model, args, omega)

    # -- properties ---------------------------------------------------------

    @property
    def max_degree(self) -> int:
        """Highest degree the loaded coefficient set supports."""
        return int(self._model.max_degree)

    @property
    def mu(self) -> float:
        """Gravitational parameter of the model, m^3 s^-2."""
        return float(self._model.mu)

    @property
    def reference_radius(self) -> float:
        """Reference radius of the expansion, m."""
        return float(self._model.r_ref)

    # -- evaluations --------------------------------------------------------

    def acceleration(self, r: Vec3, t: float, degree: int) -> Vec3:
        """Inertial acceleration truncated at ``degree``.

        One full synthesis.  Cost scales as :math:`N^2`, which is the quantity
        the campaign's budget counts.
        """
        from lunaris.physics.spherical_harmonics import sh_accel_fixed_numba

        xb, yb, zb = self._rot.to_body(r, t)
        a_body = sh_accel_fixed_numba(xb, yb, zb, degree, *self._args)
        return self._rot.to_inertial(a_body, t)

    def defect(self, r: Vec3, t: float, degree: int,
               reference_degree: int) -> Vec3:
        """Truncation defect in **one** Legendre pass.

        :math:`\\Delta\\mathbf a(t,N)=\\mathbf a_N-\\mathbf a_{\\mathrm{ref}}`
        is the input to the whole allocation, and the naive implementation --
        two calls to :meth:`acceleration` -- pays for the associated Legendre
        recursion twice even though the higher truncation already computes
        every term the lower one needs.  The kernel exposes a dual-degree
        entry point that accumulates both truncations in a single pass, and
        this uses it.

        The saving is not cosmetic: the defect table of WP1 is built at every
        accumulation epoch for every candidate degree, and it is the largest
        single computation of stage M1.

        Parameters
        ----------
        r:
            Inertial position, shape ``(3,)``, m.
        t:
            Epoch, seconds from arc start; enters through the body rotation.
        degree:
            The truncated (policy) degree :math:`N`.
        reference_degree:
            The orbit's adopted reference degree.

        Returns
        -------
        ndarray, shape (3,)
            Inertial defect, m s^-2.  Zero to machine precision when
            ``degree == reference_degree``.

        Raises
        ------
        ValueError
            If ``degree`` exceeds ``reference_degree``; the defect is defined
            as a *truncation* error and the sign convention of the manuscript
            assumes the policy degree is the lower one.
        """
        if degree > reference_degree:
            raise ValueError(
                f"degree {degree} exceeds reference degree {reference_degree}; "
                "the defect is defined for truncations of the reference field"
            )
        xb, yb, zb = self._rot.to_body(r, t)
        # The guard above fixes the ordering, so `degree` is the low
        # truncation and `reference_degree` the high one.  Passing them
        # through min/max here would look defensive but would silently flip
        # the sign of the defect if the guard were ever removed.
        ax_lo, ay_lo, az_lo, ax_hi, ay_hi, az_hi = self._dual(
            xb, yb, zb, degree, reference_degree, *self._args
        )
        d_body = (ax_lo - ax_hi, ay_lo - ay_hi, az_lo - az_hi)
        return self._rot.to_inertial(d_body, t)


def _load_dual_kernel():
    """Return the kernel's dual-degree entry point.

    Imported lazily and by name so that the dependency is visible in a
    traceback rather than hidden behind a wildcard, and so that importing
    :mod:`tda.field` does not require the kernel to be installed -- the
    analytic test fields do not need it.
    """
    from lunaris.physics.spherical_harmonics import (
        _compute_sh_acceleration_dual_numba as dual,
    )

    return dual


# ---------------------------------------------------------------------------
# Gravity gradient
# ---------------------------------------------------------------------------


def gravity_gradient(field: GravityField, r: Vec3, t: float,
                     degree: int, cfg: GradientConfig) -> Mat3:
    """Gravity gradient :math:`\\partial\\mathbf a/\\partial\\mathbf r`.

    Estimated by central differences, six syntheses per call [Press2007]_.
    This is the dominant cost of building a state-transition matrix and the
    reason a variational integration is not free: at gradient degree equal to
    the reference degree, a right-hand side costs about seven syntheses rather
    than one.

    Symmetry
    --------
    The exact gradient is the negated Hessian of a potential and is therefore
    symmetric.  A central-difference estimate is not, to roundoff.  With
    ``cfg.symmetrise`` the estimate is projected onto the symmetric part,
    which (i) halves the differencing noise, since the antisymmetric part is
    pure error, and (ii) makes the resulting linear system Hamiltonian, so
    that the state-transition matrix it generates is symplectic and
    :func:`tda.dynamics.symplectic_defect` measures integration error rather
    than gradient noise.

    Degree
    ------
    ``cfg.degree`` is deliberately explicit.  The previous campaign measured
    that a degree-120 gradient is inadequate on 31 km-perilune orbits, and
    this campaign carries a low-perilune population; see ``DECISIONS.md`` Q13.
    ``None`` means "match the trajectory's reference degree", the conservative
    choice and the default.

    Parameters
    ----------
    field:
        Gravity field adapter; only :meth:`GravityField.acceleration` is used,
        which is what lets an analytic field be substituted for verification.
    r:
        Inertial position, shape ``(3,)``, m.
    t:
        Epoch, seconds from arc start; enters through the body rotation.
    degree:
        Reference degree of the trajectory, used when ``cfg.degree is None``.
    cfg:
        Gradient settings: truncation degree, step and symmetrisation.

    Returns
    -------
    ndarray, shape (3, 3)
        Gradient in s^-2.
    """
    n = degree if cfg.degree is None else cfg.degree
    h = cfg.step_m
    g = np.empty((3, 3), dtype=float)
    for j in range(3):
        step = np.zeros(3)
        step[j] = h
        a_plus = field.acceleration(r + step, t, n)
        a_minus = field.acceleration(r - step, t, n)
        g[:, j] = (a_plus - a_minus) / (2.0 * h)
    if cfg.symmetrise:
        g = 0.5 * (g + g.T)
    return g
