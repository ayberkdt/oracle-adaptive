"""The omitted tail: its amplitude from the spectrum, its direction from bands.

Everything about the part of the field a truncation leaves out.  The two
halves belong together because they are the two halves of one estimate: the
online probe measures the *direction* of the omitted acceleration from the
first few omitted degree bands, and restores its *magnitude* from the model's
degree variances.  That split is the paper's information statement -- the
spectrum is one-dimensional and costs nothing to carry, the texture is a field
at resolution :math:`\\pi r/N` and has to be measured in flight.

Kept out of :mod:`tda.field`, which is about evaluating the field that is
retained; this module is about the part that is not.

References
----------
.. [Kaula1966] W. M. Kaula, *Theory of Satellite Geodesy*, Blaisdell, 1966,
   §3.2 -- degree variances and the power-law rule that this module
   deliberately does *not* use for the lunar field.
.. [HeiskanenMoritz1967] W. A. Heiskanen and H. Moritz, *Physical Geodesy*,
   Freeman, 1967, §1.16 -- surface spherical harmonics and the mean square of
   their surface gradient, :math:`n(n+1)`.
.. [Rapp1986] R. H. Rapp, "Global geopotential solutions", in *Mathematical
   and Numerical Techniques in Physical Geodesy*, Springer, 1986 -- degree
   variance conventions.
"""

from __future__ import annotations

import math
from typing import Protocol

import numpy as np
from numpy.typing import NDArray

from tda.field import GravityField

__all__ = [
    "BandStack",
    "DifferencingBandStack",
    "acceleration_degree_rms",
    "coefficient_degree_variance",
    "tail_completion_factor",
]

Vec3 = NDArray[np.float64]
Arr = NDArray[np.float64]


# ---------------------------------------------------------------------------
# Spectrum: amplitude of the omitted tail
# ---------------------------------------------------------------------------


def coefficient_degree_variance(c_coeffs: Arr, s_coeffs: Arr) -> Arr:
    """Degree variance of the coefficients.

    :math:`\\sigma_n^2=\\sum_m(\\bar C_{nm}^2+\\bar S_{nm}^2)`, the
    dimensionless power the field carries at degree :math:`n`
    [Rapp1986]_.

    Parameters
    ----------
    c_coeffs, s_coeffs:
        Fully normalised coefficients, shape ``(N+1, N+1)``, lower-triangular
        in order.

    Returns
    -------
    ndarray, shape (N+1,)
        Indexed by degree; entry ``n`` is :math:`\\sigma_n^2`.

    Raises
    ------
    ValueError
        If the two arrays are not square and of equal shape.

    Notes
    -----
    Only the orders :math:`m\\le n` are summed, so a coefficient array padded
    with zeros above the diagonal gives the same answer -- but the slice is
    explicit rather than relying on the padding being zero, because a model
    read from a file is not obliged to zero it.
    """
    c = np.asarray(c_coeffs, dtype=float)
    s = np.asarray(s_coeffs, dtype=float)
    if c.shape != s.shape or c.ndim != 2 or c.shape[0] != c.shape[1]:
        raise ValueError(
            f"expected two equal square arrays, got {c.shape} and {s.shape}")
    return np.array([
        float(np.sum(c[n, : n + 1] ** 2 + s[n, : n + 1] ** 2))
        for n in range(c.shape[0])
    ])


def acceleration_degree_rms(degree_variance: Arr, radius: float,
                            reference_radius: float, mu: float) -> Arr:
    """Per-degree RMS acceleration at a radius, from the degree variances.

    The quantity :math:`\\sigma_a(n;r)` the probe's amplitude completion
    needs, and the reason it is derived here rather than assumed: the ratio in
    :func:`tail_completion_factor` is *not* scale-free.  The factors
    :math:`\\mu/r^2` cancel, but the attenuation :math:`(R/r)^n` and the
    degree-dependent gain do not, and using the bare coefficient variances in
    their place would misweight the tail by orders of magnitude at high
    degree.

    Averaging the squared gradient of a solid harmonic of degree :math:`n`
    over a sphere of radius :math:`r` gives a radial part
    :math:`(n+1)^2` and a tangential part :math:`n(n+1)`
    [HeiskanenMoritz1967]_, so

    .. math::
        \\sigma_a(n;r) \\;=\\; \\frac{\\mu}{r^{2}}
            \\Bigl(\\frac{R}{r}\\Bigr)^{n}
            \\sqrt{(n+1)(2n+1)}\\;\\sigma_n .

    Parameters
    ----------
    degree_variance:
        :math:`\\sigma_n^2`, shape ``(N+1,)``, indexed by degree.
    radius:
        Evaluation radius, m.
    reference_radius:
        Reference radius of the expansion, m.
    mu:
        Gravitational parameter, m^3 s^-2.

    Returns
    -------
    ndarray, shape (N+1,)
        RMS acceleration per degree, m s^-2.  Entry 0 is the central term and
        is returned as computed rather than special-cased; the probe never
        reads it, because it only ever sums degrees above a truncation.

    Raises
    ------
    ValueError
        If any scale is not positive, or a degree variance is negative.

    Notes
    -----
    The attenuation and the coefficient magnitude are combined in logarithms
    and exponentiated once.  Not because the direct form overflows -- at a
    100 km lunar orbit :math:`(R/r)^{1800}\\approx10^{-45}`, comfortably
    representable -- but because it lets a zero degree variance pass through
    as an exact zero without a divide warning, and keeps the two factors from
    being rounded separately before they are multiplied.  Neither is dramatic;
    both are free.
    """
    if radius <= 0.0 or reference_radius <= 0.0 or mu <= 0.0:
        raise ValueError(
            f"radius, reference_radius and mu must be positive; got "
            f"{radius}, {reference_radius}, {mu}")
    sigma2 = np.asarray(degree_variance, dtype=float)
    if np.any(sigma2 < 0.0):
        raise ValueError("degree variances must be non-negative")

    n = np.arange(sigma2.size, dtype=float)
    log_attenuation = n * math.log(reference_radius / radius)
    gain = np.sqrt((n + 1.0) * (2.0 * n + 1.0))
    with np.errstate(divide="ignore"):
        log_sigma = 0.5 * np.log(sigma2)          # -inf where sigma2 == 0
    return (mu / radius**2) * gain * np.exp(log_sigma + log_attenuation)


def tail_completion_factor(sigma_a: Arr, degree: int, depth: int) -> float:
    """Amplitude completion factor :math:`\\gamma(N,r)`.

    The probe measures direction from the first ``depth`` omitted bands;
    :math:`\\gamma` restores the magnitude of the whole omitted tail,

    .. math::
        \\gamma = \\Bigl[\\frac{\\sum_{n>N}\\sigma_a^2(n;r)}
                              {\\sum_{n=N+1}^{N+k}\\sigma_a^2(n;r)}\\Bigr]^{1/2}.

    ``sigma_a`` is the *measured* per-degree acceleration RMS of
    :func:`acceleration_degree_rms`, not a fitted power law.  The lunar
    spectrum is not a single power law -- the previous campaign measured a
    spectral slope of 2.13 against an effective tail exponent of 1.76 for the
    same field -- so a Kaula-type extrapolation [Kaula1966]_ carries a
    systematic error that a table of degree variances does not.  The power-law
    variant is kept as the ablation ``abl-kaula``, not as the default.

    Parameters
    ----------
    sigma_a:
        Per-degree acceleration RMS at the evaluation radius, indexed so that
        ``sigma_a[n]`` is degree ``n``.
    degree:
        Truncation degree :math:`N`, non-negative.
    depth:
        Probe depth :math:`k`, positive.

    Returns
    -------
    float
        The completion factor, :math:`\\ge 1` by construction since the probed
        bands are a subset of the tail.

    Raises
    ------
    ValueError
        If ``depth`` is not positive, ``degree`` is negative, or the probed
        window is empty or carries no power -- each of which makes the ratio
        undefined rather than large.

    Examples
    --------
    A probe that covers the whole tail needs no completion.

    >>> import numpy as np
    >>> tail_completion_factor(np.array([0.0, 1.0, 0.5, 0.25]), 1, 2)
    1.0
    """
    if depth <= 0:
        raise ValueError(f"probe depth must be positive, got {depth}")
    if degree < 0:
        raise ValueError(f"degree must be non-negative, got {degree}")
    sigma_a = np.asarray(sigma_a, dtype=float)

    tail = sigma_a[degree + 1:]
    probed = sigma_a[degree + 1: degree + 1 + depth]
    if probed.size == 0:
        raise ValueError(
            f"probe window is empty: degree={degree}, depth={depth}, "
            f"spectrum has {sigma_a.size} entries")
    denominator = float(probed @ probed)
    if denominator <= 0.0:
        raise ValueError(
            f"probed bands {degree + 1}..{degree + depth} carry no power; "
            "gamma is undefined")
    return math.sqrt(float(tail @ tail) / denominator)


# ---------------------------------------------------------------------------
# Bands: direction of the omitted tail
# ---------------------------------------------------------------------------


class BandStack(Protocol):
    """Supplies the leading omitted degree bands at a point.

    The probe needs :math:`\\sum_{n=N+1}^{N+k}\\mathbf a_n(t)` for every
    candidate degree :math:`N` in a window.  Expressed through cumulative
    truncations that is a difference of two partial sums, so what a caller
    actually needs is :math:`\\mathbf a_{\\le n}` across the window.
    """

    def cumulative(self, r: Vec3, t: float,
                   degrees: NDArray[np.int_]) -> Arr:
        """Return :math:`\\mathbf a_{\\le n}` for each ``n`` in ``degrees``.

        Parameters
        ----------
        r:
            Inertial position, shape ``(3,)``, m.
        t:
            Epoch, seconds from arc start.
        degrees:
            Truncation degrees to evaluate, shape ``(P,)``.

        Returns
        -------
        ndarray, shape (P, 3)
        """

    def syntheses_for(self, n_degrees: int) -> int:
        """Cost model: syntheses one :meth:`cumulative` call would take.

        A pure function of the request, so a design can be priced before it is
        run.  The manuscript's overhead accounting turns on whether this is
        one or one-per-degree.
        """

    @property
    def total_syntheses(self) -> int:
        """Syntheses actually performed; the measured counterpart."""


class DifferencingBandStack:
    """Bands by differencing independent fixed-degree syntheses.

    The *correct but expensive* implementation, and the honest baseline.  The
    manuscript argues that a shared band stack costs one Legendre recursion
    because the accumulator is already summed degree by degree; that is true
    of the mathematics but **not** of the kernel as it stands, which exposes
    fixed-degree and dual-degree entry points and no cumulative-by-degree one.
    Adding that entry point is a small change inside the existing accumulation
    loop, but it is new kernel work and must be validated -- the summed bands
    must reproduce the fixed-degree result to machine precision -- before the
    manuscript's cost model may be quoted (``DECISIONS.md`` D120).

    Until then this supplies the same numbers at one synthesis per requested
    degree and reports that, so the overhead entering :math:`B_+` is the
    overhead of the code that ran.
    """

    __slots__ = ("_field", "_total")

    def __init__(self, field: GravityField) -> None:
        """Wrap ``field`` and start the synthesis counter at zero."""
        self._field = field
        self._total = 0

    def cumulative(self, r: Vec3, t: float,
                   degrees: NDArray[np.int_]) -> Arr:
        """Evaluate the field once per requested degree.

        Raises
        ------
        ValueError
            If ``degrees`` is empty.  An empty stack returns an empty array
            that a caller would sum to zero, which reads as "no omitted
            acceleration" rather than as a bug.
        """
        degrees = np.atleast_1d(np.asarray(degrees, dtype=int))
        if degrees.size == 0:
            raise ValueError("degrees is empty; nothing to evaluate")
        out = np.empty((degrees.size, 3), dtype=float)
        for i, n in enumerate(degrees):
            out[i] = self._field.acceleration(r, t, int(n))
        self._total += degrees.size
        return out

    def syntheses_for(self, n_degrees: int) -> int:
        """One synthesis per degree -- the price of not having a shared stack."""
        return int(n_degrees)

    @property
    def total_syntheses(self) -> int:
        """Full syntheses performed since construction; the measured cost."""
        return self._total
