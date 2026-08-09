"""The band probe: estimating the *direction* of what the truncation omits.

A tail criterion predicts :math:`\\lVert\\Delta\\mathbf a\\rVert`.  The
coupling term needs :math:`\\Delta\\mathbf a` itself, and a controller given
only the magnitude can reach the sensitivity-weighted rung and no further, by
construction.  What is affordable is to evaluate the first few omitted bands
and use their sum as the direction, completing the amplitude from the model's
own degree-variance spectrum:

.. math::
    \\hat{\\mathbf v}(t,N)\\;=\\;-\\,\\gamma(N,r)\\!\\!
    \\sum_{n=N+1}^{N+k}\\!\\!\\mathbf a_n(t) .

That split is the paper's information statement in its sharpest form.  The
spectrum is one-dimensional and costs nothing to upload; the texture is a
field at resolution :math:`\\pi r/N` and has to be measured in flight.

Forward, not backward
---------------------
The probe is taken *ahead* of the boundary, over the interval about to be
committed to.  Accumulating probes inside an interval and applying them at the
next boundary would select a degree from a direction that has already
decorrelated.  Two facts make forward probing a small requirement rather than
a circularity: the probe locations barely depend on the degree, because two
candidates separate by at most
:math:`\\tfrac12\\lVert\\Delta\\mathbf a\\rVert\\Delta t_{\\mathrm{dec}}^2`
over one interval --- millimetres against a kilometre-scale texture --- and
what the predictor must achieve is a position inside that texture, which an
analytic two-body step delivers with no field evaluation at all
(:mod:`tda.kepler`).

What this module does not decide
--------------------------------
The cost.  Whether the shared band stack costs one Legendre recursion or one
per degree is a property of the kernel, not of this code, and the kernel has
no cumulative-by-degree entry point yet (``DECISIONS.md`` D120).  Everything
here asks a :class:`~tda.spectrum.BandStack` and reports what that stack
charged, so the overhead entering :math:`B_+` is the overhead of the code that
ran.  :func:`probe_overhead` makes the degree dependence explicit, because the
overhead is proportional to :math:`N` and was for a while quoted as a single
range (D141).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from tda.kepler import propagate_two_body
from tda.spectrum import BandStack, tail_completion_factor

__all__ = [
    "ProbePlan",
    "band_direction",
    "candidate_window",
    "direction_accuracy",
    "plan_probe_points",
    "probe_overhead",
    "required_degrees",
    "retained_fraction",
]

Arr = NDArray[np.float64]
IntArr = NDArray[np.int_]


# ---------------------------------------------------------------------------
# The candidate window
# ---------------------------------------------------------------------------


def candidate_window(candidates: tuple[int, ...], planned: int,
                     half_width: int) -> tuple[int, ...]:
    """Degrees the online decision searches around the plan's nominal.

    The half-width is in **candidate-grid indices**, not in degrees.  Fixing
    it that way keeps the controller's complexity --- :math:`2\\delta+1`
    candidates per decision --- invariant when the grid is refined or
    coarsened; the induced span in degrees follows from the grid and is
    reported alongside the cost, since it is the span and not the index count
    that sets the band stack.

    Parameters
    ----------
    candidates:
        Sorted, unique candidate degrees.
    planned:
        The offline plan's degree for this interval; must be a candidate.
    half_width:
        :math:`\\delta\\ge0`.

    Returns
    -------
    tuple of int
        Clipped at the ends of the grid, so a window near the edge is shorter
        rather than wrapping or extrapolating.

    Raises
    ------
    ValueError
        If the half-width is negative or the planned degree is not tabulated.

    Examples
    --------
    >>> candidate_window((10, 20, 30, 40, 50), 30, 1)
    (20, 30, 40)
    >>> candidate_window((10, 20, 30, 40, 50), 10, 2)
    (10, 20, 30)
    """
    if half_width < 0:
        raise ValueError(f"half_width must be non-negative, got {half_width}")
    try:
        centre = candidates.index(planned)
    except ValueError:
        raise ValueError(
            f"planned degree {planned} is not a tabulated candidate"
        ) from None
    lo = max(0, centre - half_width)
    hi = min(len(candidates), centre + half_width + 1)
    return candidates[lo:hi]


def required_degrees(window: tuple[int, ...], depth: int) -> IntArr:
    """Cumulative truncations the window's directions can be built from.

    Each candidate :math:`N` needs :math:`\\mathbf a_{\\le N}` and
    :math:`\\mathbf a_{\\le N+k}`, since the probed band sum is their
    difference.  Returning the union, sorted and de-duplicated, is what lets
    one stack serve every candidate: the candidates are not independent
    probes but partial sums of one evaluation.

    Parameters
    ----------
    window:
        Candidate degrees.
    depth:
        Probe depth :math:`k`, positive.

    Returns
    -------
    ndarray of int, ascending.

    Raises
    ------
    ValueError
        If the window is empty or the depth is not positive.

    Examples
    --------
    >>> required_degrees((20, 30), 3).tolist()
    [20, 23, 30, 33]
    """
    if not window:
        raise ValueError("the candidate window is empty")
    if depth <= 0:
        raise ValueError(f"probe depth must be positive, got {depth}")
    degrees = np.asarray(window, dtype=np.int64)
    return np.unique(np.concatenate([degrees, degrees + depth]))


# ---------------------------------------------------------------------------
# Probe geometry
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class ProbePlan:
    """Where and when the probe will evaluate the field.

    Attributes
    ----------
    times:
        Epochs, shape ``(n,)``, inside the interval about to be committed to.
    positions:
        Predicted inertial positions there, shape ``(n, 3)``.
    predictor_error_bound:
        :math:`\\tfrac12\\lVert\\mathbf a_{\\mathrm{pert}}\\rVert\\Delta t^2`
        for the largest step taken, metres --- the analytic bound on how far
        the two-body prediction can be from the flown path, to be read against
        :math:`\\pi r/N`.  Carried rather than recomputed downstream so that
        the number in the manifest is the one the plan was made with.
    """

    times: Arr
    positions: Arr
    predictor_error_bound: float

    def __len__(self) -> int:
        """Number of probe points."""
        return int(self.times.size)


def plan_probe_points(state: Arr, mu: float, epoch: float, span: float,
                      n_points: int,
                      perturbing_acceleration: float = 0.0) -> ProbePlan:
    """Predict where to probe over the coming interval, analytically.

    Parameters
    ----------
    state:
        Current inertial state ``[r, v]``, shape ``(6,)``.
    mu:
        Gravitational parameter.
    epoch:
        Time at the boundary, seconds from arc start.
    span:
        Interval length :math:`\\Delta t_{\\mathrm{dec}}`.
    n_points:
        How many points to place.  They are placed at cell midpoints of a
        uniform split of the interval, not at its ends: an endpoint sample
        would give the boundary twice across consecutive intervals and leave
        the middle unprobed.
    perturbing_acceleration:
        Magnitude used for :attr:`ProbePlan.predictor_error_bound`.  Zero
        leaves the bound at zero, which is honest for a caller that has not
        supplied one rather than a claim that the predictor is exact.

    Returns
    -------
    ProbePlan

    Raises
    ------
    ValueError
        If the span or the point count is not positive.

    Notes
    -----
    Each point is propagated from the *boundary* state rather than from the
    previous point, so the errors do not compound along the interval and the
    largest step --- hence the largest error --- is the last one.
    """
    if span <= 0.0:
        raise ValueError(f"span must be positive, got {span}")
    if n_points <= 0:
        raise ValueError(f"n_points must be positive, got {n_points}")

    offsets = (np.arange(n_points) + 0.5) * (span / n_points)
    positions = np.stack([propagate_two_body(state, mu, float(dt))[0:3]
                          for dt in offsets])
    bound = 0.5 * perturbing_acceleration * float(offsets[-1]) ** 2
    return ProbePlan(times=epoch + offsets, positions=positions,
                     predictor_error_bound=bound)


# ---------------------------------------------------------------------------
# The direction estimate
# ---------------------------------------------------------------------------


def band_direction(stack: BandStack, position: Arr, epoch: float,
                   window: tuple[int, ...], depth: int,
                   sigma_a: Arr) -> dict[int, Arr]:
    """Estimated omitted acceleration for every candidate, from one stack.

    Parameters
    ----------
    stack:
        Supplies :math:`\\mathbf a_{\\le n}`; the cost of that is the stack's
        business and is recorded on it.
    position, epoch:
        Where and when.
    window:
        Candidate degrees.
    depth:
        Probe depth :math:`k`.
    sigma_a:
        Per-degree acceleration RMS at this radius, from
        :func:`tda.spectrum.acceleration_degree_rms`, indexed by degree.

    Returns
    -------
    dict
        Candidate degree to :math:`\\hat{\\mathbf v}`, shape ``(3,)``.

    Raises
    ------
    ValueError
        If a candidate's completion factor is undefined --- which happens when
        the probed bands lie outside the tabulated spectrum, and is a
        configuration error rather than a small number.

    Notes
    -----
    The sign is negative: :math:`\\hat{\\mathbf v}` estimates
    :math:`\\Delta\\mathbf a=\\mathbf a_N-\\mathbf a_{\\mathrm{ref}}`, which is
    minus the omitted tail, so the probed bands enter with a minus.  Getting
    that backwards would leave the magnitude untouched and reverse every
    cancellation decision, which is exactly the failure the ``abl-null``
    control is designed to expose.
    """
    degrees = required_degrees(window, depth)
    cumulative = stack.cumulative(position, epoch, degrees)
    lookup = {int(n): cumulative[i] for i, n in enumerate(degrees)}

    out: dict[int, Arr] = {}
    for candidate in window:
        band_sum = lookup[candidate + depth] - lookup[candidate]
        gamma = tail_completion_factor(sigma_a, candidate, depth)
        out[candidate] = -gamma * band_sum
    return out


# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------


def direction_accuracy(estimate: Arr, truth: Arr) -> float:
    """:math:`\\kappa=\\cos\\angle(\\hat{\\mathbf v},\\Delta\\mathbf a)`.

    A **diagnostic**, not a gate.  It is tempting to read it as the fraction
    of the coupling term the probe retains, and it is not: the retained
    fraction depends on the angles the true and estimated defects make with
    :math:`\\mathbf z=(\\Phi\\mathbf B)^{\\top}\\mathbf c`, not on the angle
    between them, because :math:`\\Phi\\mathbf B` does not preserve angles
    (``DECISIONS.md`` D74).  A hard gate on it would risk an early false
    negative, so the deployability decision rests on the capture fraction
    instead (D106).

    Returns
    -------
    float
        In :math:`[-1,1]`; ``nan`` if either vector vanishes, since the angle
        is then undefined and returning zero would read as "orthogonal".
    """
    estimate = np.asarray(estimate, dtype=float)
    truth = np.asarray(truth, dtype=float)
    norms = np.linalg.norm(estimate) * np.linalg.norm(truth)
    if norms == 0.0:
        return float("nan")
    return float(estimate @ truth / norms)


def retained_fraction(estimate: Arr, truth: Arr, z: Arr) -> float:
    """Fraction of the coupling term the probe actually retains.

    :math:`\\langle\\hat{\\mathbf v},\\mathbf z\\rangle/
    \\langle\\Delta\\mathbf a,\\mathbf z\\rangle` --- the quantity
    :math:`\\kappa` is mistaken for.

    Parameters
    ----------
    estimate, truth:
        :math:`\\hat{\\mathbf v}` and :math:`\\Delta\\mathbf a`, shape
        ``(3,)``.
    z:
        :math:`(\\Phi\\mathbf B)^{\\top}\\mathbf c`, shape ``(3,)``, with
        :math:`\\mathbf c` from the sensitivity-weighted schedule so that the
        two diagnostics are computed in the same run.

    Returns
    -------
    float
        ``nan`` where the true contribution vanishes --- there the term is
        zero and any probe error produces a spurious one, so a ratio would be
        meaningless rather than large.
    """
    denominator = float(np.asarray(truth, dtype=float) @ z)
    if denominator == 0.0:
        return float("nan")
    return float(np.asarray(estimate, dtype=float) @ z) / denominator


def probe_overhead(degree: int, speed: float, radius: float,
                   duration: float, rhs_calls: int,
                   points_per_correlation_time: float = 1.0) -> float:
    """Share of the run's right-hand-side calls the probe adds.

    .. math::
        \\frac{B_{\\mathrm{probe}}}{B_{\\mathrm{tot}}}=\\frac{1}
        {N_{\\mathrm{RHS}}}\\int_0^T\\frac{\\dd t}{\\tau_{\\mathrm{corr}}(t)}
        \\;=\\;\\frac{1}{N_{\\mathrm{RHS}}}\\int_0^T
        \\frac{N(t)v(t)}{\\pi r(t)}\\,\\dd t .

    Proportional to the degree, and that is the point of exposing it as a
    function rather than a constant.  The same integral gives the accumulation
    grid's cell count, so the two are the same quantity up to
    :math:`n_s`; a figure quoted without the degree it belongs to is not a
    property of the arc (D141).

    Parameters
    ----------
    degree, speed, radius:
        Arrays or scalars sampled along the arc; broadcast together.
    duration:
        Arc length, used only to normalise a scalar sample into an integral.
    rhs_calls:
        :math:`N_{\\mathrm{RHS}}`, the realized call count of the run.
    points_per_correlation_time:
        Probe points per correlation time.  One is the default and the
        manuscript's model; ``C-lite`` is the other end of the curve and is
        priced separately, since its single probe is co-located with a call
        the integrator makes anyway.

    Returns
    -------
    float
        A fraction, not a percentage.

    Raises
    ------
    ValueError
        If the call count or the duration is not positive.
    """
    if rhs_calls <= 0:
        raise ValueError(f"rhs_calls must be positive, got {rhs_calls}")
    if duration <= 0.0:
        raise ValueError(f"duration must be positive, got {duration}")
    rate = np.asarray(degree, dtype=float) * np.asarray(speed, dtype=float) / (
        np.pi * np.asarray(radius, dtype=float))
    return float(np.mean(rate) * duration * points_per_correlation_time
                 / rhs_calls)
