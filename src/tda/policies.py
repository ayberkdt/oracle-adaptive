"""Comparator schedules: the constant family and the radial family.

These are the policies the trajectory-aware allocation is measured against.
None of them optimises the objective of :mod:`tda.kernel` --- that is what
makes them comparators --- so all of them are cheap, and all of them are
generated here rather than solved.

* :func:`nearest_constant` is ``F-op``, the degree a user selects in advance.
  It is also the **anchor**: the campaign's work ceiling is its realized work,
  so it is the one policy not constrained by a budget it defines.
* :func:`radial_family` is the one-parameter interpolation between the
  constant and the radial rule.  ``k=0`` returns the constant, ``k=1`` the
  radial endpoint ``R-rad``, and ``k=0.5`` the interior member ``R-int``,
  which the previous campaign found to be the strongest deployable schedule
  available.
* :func:`fixed_family_envelope` is ``F-env``, the post-hoc lower envelope of
  the constant family.  It is not a policy: it reads every orbit's outcome
  before choosing, and exists to bound how much of the constant family's
  potential ``F-op`` leaves unused.

The family, and why it is geometric
-----------------------------------
For a constant degree :math:`N_0` and a radial table :math:`N_A(h)`,

.. math::
    N_k(h) \\;=\\; \\operatorname{round}\\bigl(s_k\\,
    N_0^{1-k}\\,N_A(h)^{k}\\bigr),

with :math:`s_k` calibrated to the budget.  Geometric rather than linear
because it makes the degree *span* multiplicative: :math:`\\operatorname{span}
(k)=\\operatorname{span}(1)^k`, so :math:`k` moves the aggressiveness on a
scale whose endpoints are a flat profile and the rule's own profile.  The
previous campaign introduced this family and reported :math:`k=0.5`; the
construction is reproduced unchanged so that the comparison inherits its
calibration rather than a new one.

Two places where this campaign deliberately differs
---------------------------------------------------
**The budget is time-weighted.**  The previous campaign calibrated
:math:`\\langle N^2\\rangle` as a plain mean over a uniform output grid.  Here
the accumulation grid is refined on the correlation time, so a plain mean
would charge perilune for most of the arc merely because it is sampled most
densely; the work is :math:`\\sum_q W_qN_q^2` with :math:`W_q` the interval's
time weight (``DECISIONS.md`` D60).

**A comparator is calibrated to the ceiling, not below it.**  The budget is a
ceiling and an *optimiser* may rationally leave it unspent, since the
objective is not monotone in degree (D142).  A comparator has no such reason:
it does not optimise the objective, so spending less would only make it
artificially weak.  These policies therefore take the largest scale whose work
still fits, which is the tight end of the same contract rather than a
different one.

References
----------
.. [Atallah2022] A. M. Atallah et al., "Radial adaptive gravity model
   truncation", 2022 -- the radial rule the endpoint reproduces.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

__all__ = [
    "RadialTable",
    "ScheduleCalibration",
    "fixed_family_envelope",
    "interval_altitudes",
    "nearest_constant",
    "radial_family",
    "snap_to_candidates",
]

Arr = NDArray[np.float64]
IntArr = NDArray[np.int_]

BISECTION_STEPS = 60
"""Enough to exhaust a double's mantissa on the scale bracket."""


# ---------------------------------------------------------------------------
# Candidate handling
# ---------------------------------------------------------------------------


def snap_to_candidates(degrees: Arr, candidates: tuple[int, ...]) -> IntArr:
    """Round each degree to the nearest tabulated candidate.

    A schedule that names a degree the defect table was not built for cannot
    be evaluated at all, so the snap happens here, once, rather than failing
    later in a lookup.  It is a rounding of the *policy*, not of the physics:
    the candidate set is declared in advance and contains each comparator's
    own degrees by construction (``NOTATION.md`` §1), so the snap is normally
    exact and the campaign reports the cases where it is not.

    Parameters
    ----------
    degrees:
        Real-valued degrees, any shape.
    candidates:
        Sorted, unique, the table's candidate axis.

    Returns
    -------
    ndarray of int
        Same shape as ``degrees``.

    Raises
    ------
    ValueError
        If the candidate set is empty or unsorted -- an unsorted set would
        make the nearest-neighbour search silently wrong.

    Examples
    --------
    >>> import numpy as np
    >>> snap_to_candidates(np.array([37.0, 44.0]), (30, 40, 50))
    array([40, 40])
    """
    grid = np.asarray(candidates, dtype=float)
    if grid.size == 0:
        raise ValueError("candidate set is empty")
    if np.any(np.diff(grid) <= 0.0):
        raise ValueError("candidate set must be sorted and unique")
    idx = np.clip(np.searchsorted(grid, degrees), 1, grid.size - 1)
    take_lower = np.abs(degrees - grid[idx - 1]) <= np.abs(grid[idx] - degrees)
    return np.where(take_lower, grid[idx - 1], grid[idx]).astype(np.int64)


def interval_altitudes(radius: Arr, widths: Arr, interval_of: IntArr,
                       n_intervals: int, body_radius: float) -> Arr:
    """Time-weighted mean altitude of each decision interval, in metres.

    The radial rule is a function of altitude and the altitude varies inside
    an interval, so reducing the rule to one degree per interval requires a
    representative altitude and the choice has to be declared.  This is the
    time-weighted mean, which is the reduction that leaves the *interval's*
    share of the arc unchanged; a midpoint sample would over-represent
    whichever end the refinement happened to place a node near.

    Parameters
    ----------
    radius:
        Selenocentric radius at each accumulation node, shape ``(M,)``.
    widths:
        Cell widths, shape ``(M,)``.
    interval_of:
        :math:`g(i)`, shape ``(M,)``.
    n_intervals:
        :math:`K_{\\mathrm{dec}}`.
    body_radius:
        Reference radius subtracted to give altitude.

    Returns
    -------
    ndarray, shape (n_intervals,)
    """
    weight = np.bincount(interval_of, weights=widths, minlength=n_intervals)
    if np.any(weight <= 0.0):
        raise ValueError("some decision interval has no accumulation cell")
    weighted = np.bincount(interval_of, weights=widths * radius,
                           minlength=n_intervals)
    return weighted / weight - body_radius


@dataclass(frozen=True, slots=True)
class RadialTable:
    """Altitude-binned degree lookup, the radial rule of [Atallah2022]_.

    Attributes
    ----------
    bin_edges_m:
        Lower edge of each altitude bin, ascending, metres.
    degrees:
        Degree for each bin, same length.

    Notes
    -----
    Altitudes below the first bin take the first degree and those above the
    last take the last, matching the archive's clamped lookup.  The clamp is
    not a convenience: extrapolating a binned table past its range would
    invent a rule the previous campaign never flew.
    """

    bin_edges_m: Arr
    degrees: IntArr

    def __post_init__(self) -> None:
        """Reject a table that cannot be looked up unambiguously."""
        if self.bin_edges_m.shape != self.degrees.shape:
            raise ValueError(
                f"bin_edges_m {self.bin_edges_m.shape} and degrees "
                f"{self.degrees.shape} must match")
        if self.bin_edges_m.size == 0:
            raise ValueError("radial table is empty")
        if np.any(np.diff(self.bin_edges_m) <= 0.0):
            raise ValueError("bin edges must be strictly ascending")

    def at(self, altitude_m: Arr) -> Arr:
        """Degree at each altitude, clamped to the table's range."""
        idx = np.clip(np.searchsorted(self.bin_edges_m, altitude_m,
                                      side="right") - 1,
                      0, self.bin_edges_m.size - 1)
        return self.degrees[idx].astype(float)


# ---------------------------------------------------------------------------
# Work and calibration
# ---------------------------------------------------------------------------


def _work(degrees: IntArr, time_weight: Arr) -> float:
    """:math:`\\sum_q W_qN_q^2`, the quantity the budget bounds."""
    return float(time_weight @ (degrees.astype(float) ** 2))


@dataclass(frozen=True, slots=True)
class ScheduleCalibration:
    """A schedule and the evidence that it was calibrated, not assumed.

    Attributes
    ----------
    degrees:
        One degree per decision interval, shape ``(K_dec,)``.
    scale:
        The multiplier :math:`s_k` the bisection settled on.
    work:
        Realized nominal work :math:`\\sum_q W_qN_q^2`.
    budget:
        The ceiling it was calibrated against.
    scale_limited:
        Whether raising the scale any further would exceed the ceiling.  This
        is the property of the *calibration* and is normally true.
    step_limited:
        Whether raising any single interval to its next candidate would also
        exceed it.  Normally **false**, and that is not a defect: a
        one-parameter family moves every interval together, so on a discrete
        candidate grid it generally cannot land exactly on the ceiling.  The
        gap is reported rather than closed, because closing it would mean
        leaving the family and inventing a schedule the previous campaign
        never flew.  It is also why the matching contract records realized
        work beside every error (D142) instead of assuming equality.
    """

    degrees: IntArr
    scale: float
    work: float
    budget: float
    scale_limited: bool
    step_limited: bool

    @property
    def utilisation(self) -> float:
        """Fraction of the ceiling actually spent."""
        return self.work / self.budget


def nearest_constant(budget: float, duration: float,
                     candidates: tuple[int, ...]) -> int:
    """``F-op``: the candidate degree nearest the budget.

    The degree a user can select in advance, and the campaign's anchor.  It is
    *nearest*, not *largest feasible*: the anchor defines the work ceiling
    rather than obeying one, and rounding it down would move the ceiling every
    time the candidate grid changed.  Every other policy is then calibrated to
    the anchor's realized work.

    Parameters
    ----------
    budget:
        Nominal work :math:`B=B_1T`.
    duration:
        Arc length :math:`T`; a constant degree costs :math:`TN^2`.
    candidates:
        The tabulated degree set.

    Returns
    -------
    int

    Raises
    ------
    ValueError
        If the budget or duration is not positive.

    Examples
    --------
    >>> nearest_constant(4.0e4 * 10.0, 10.0, (100, 200, 300))
    200
    """
    if budget <= 0.0 or duration <= 0.0:
        raise ValueError(
            f"budget and duration must be positive, got {budget}, {duration}")
    target = math.sqrt(budget / duration)
    return int(snap_to_candidates(np.array([target]), candidates)[0])


def radial_family(k: float, constant_degree: int, radial_degrees: Arr,
                  time_weight: Arr, budget: float,
                  candidates: tuple[int, ...],
                  reference_degree: int) -> ScheduleCalibration:
    """Interpolate geometrically between constant and radial, then calibrate.

    .. math::
        N_k \\;=\\; \\operatorname{round}\\bigl(s\\,N_0^{1-k}N_A^{k}\\bigr),

    with :math:`s` bisected to the largest value whose work fits the ceiling.

    Parameters
    ----------
    k:
        Family parameter in :math:`[0,1]`.  ``0`` is the constant endpoint,
        ``1`` the radial endpoint ``R-rad``, ``0.5`` the interior member
        ``R-int``.
    constant_degree:
        :math:`N_0`, normally :func:`nearest_constant`.
    radial_degrees:
        :math:`N_A` per decision interval, shape ``(K_dec,)``, from
        :meth:`RadialTable.at` on :func:`interval_altitudes`.
    time_weight:
        :math:`W_q`, shape ``(K_dec,)``.
    budget:
        The work ceiling.
    candidates:
        Tabulated degree set; the result is snapped to it.
    reference_degree:
        Upper clamp.  A policy that reached the reference degree would have a
        zero defect there by construction, which the campaign censors rather
        than scores, so the clamp keeps the schedule inside the range where
        the comparison means something.

    Returns
    -------
    ScheduleCalibration

    Raises
    ------
    ValueError
        If ``k`` is outside :math:`[0,1]`, the shapes disagree, or no scale in
        the bracket produces a schedule that fits the budget -- which happens
        when even the smallest candidate exceeds it, and is a statement about
        the budget rather than about the policy.

    Notes
    -----
    Bisecting on a rounded, snapped, clamped schedule means the work is a step
    function of the scale, not a continuous one.  The bisection therefore
    tracks the largest *feasible* scale seen rather than converging on a root,
    and the returned schedule is the one that scale produced.  Chasing a root
    that may not exist is how a calibration loop silently returns an
    infeasible schedule.
    """
    if not 0.0 <= k <= 1.0:
        raise ValueError(f"k must lie in [0, 1], got {k}")
    radial_degrees = np.asarray(radial_degrees, dtype=float)
    time_weight = np.asarray(time_weight, dtype=float)
    if radial_degrees.shape != time_weight.shape:
        raise ValueError(
            f"radial_degrees {radial_degrees.shape} and time_weight "
            f"{time_weight.shape} must match")
    if budget <= 0.0:
        raise ValueError(f"budget must be positive, got {budget}")
    if np.any(radial_degrees <= 0.0) or constant_degree <= 0:
        raise ValueError("degrees must be positive for a geometric blend")

    shape = (float(constant_degree) ** (1.0 - k)) * (radial_degrees ** k)
    floor = float(min(candidates))

    def build(scale: float) -> IntArr:
        raw = np.clip(scale * shape, floor, float(reference_degree))
        return snap_to_candidates(raw, candidates)

    lo, hi = 0.0, 1.0
    while _work(build(hi), time_weight) <= budget and hi < 1.0e6:
        hi *= 2.0
    best = build(lo) if _work(build(lo), time_weight) <= budget else None
    best_scale = lo
    for _ in range(BISECTION_STEPS):
        mid = 0.5 * (lo + hi)
        trial = build(mid)
        if _work(trial, time_weight) <= budget:
            best, best_scale, lo = trial, mid, mid
        else:
            hi = mid

    if best is None:
        raise ValueError(
            f"no scale fits the budget {budget:.6g}: even the smallest "
            f"candidate degree {int(floor)} costs "
            f"{_work(build(0.0), time_weight):.6g}")

    work = _work(best, time_weight)
    nudged = build(best_scale * (1.0 + 1.0e-9) + 1.0e-12)
    return ScheduleCalibration(
        degrees=best,
        scale=best_scale,
        work=work,
        budget=budget,
        scale_limited=_work(nudged, time_weight) > budget or hi <= lo,
        step_limited=_is_step_limited(best, time_weight, work, budget,
                                      candidates),
    )


def _is_step_limited(degrees: IntArr, time_weight: Arr, work: float,
                     budget: float, candidates: tuple[int, ...]) -> bool:
    """Report whether raising any single interval one candidate overshoots.

    Not the same as adding one to the degree: with candidates ten apart,
    ``N+1`` snaps back to ``N`` and every schedule would look saturated.  The
    cheapest single step is the one that decides it.
    """
    grid = np.asarray(candidates, dtype=float)
    nxt = np.clip(np.searchsorted(grid, degrees.astype(float), side="right"),
                  0, grid.size - 1)
    raised = grid[nxt]
    increment = time_weight * (raised**2 - degrees.astype(float) ** 2)
    movable = increment > 0.0
    if not movable.any():
        return True                      # every interval is at the top
    return bool(work + increment[movable].min() > budget)


def fixed_family_envelope(degrees: tuple[int, ...],
                          errors: Arr) -> tuple[int, float]:
    """``F-env``: the best constant degree, chosen after seeing the errors.

    Not a policy and never reported as one.  It reads each orbit's outcome
    before choosing, so no deployable rule is expected to match it; it is
    carried to bound how much of the constant family's potential ``F-op``
    leaves unused, and no hypothesis of the campaign requires beating it.

    Parameters
    ----------
    degrees:
        The constant degrees evaluated, shape ``(P,)``.
    errors:
        Their errors on one orbit, same length.  Non-finite entries are
        treated as unavailable rather than as large: a censored cell must not
        be able to win by being unmeasured.

    Returns
    -------
    tuple of (int, float)
        The winning degree and its error.

    Raises
    ------
    ValueError
        If the lengths disagree or every entry is non-finite.
    """
    errors = np.asarray(errors, dtype=float)
    if errors.shape != (len(degrees),):
        raise ValueError(
            f"errors {errors.shape} must match {len(degrees)} degrees")
    usable = np.isfinite(errors)
    if not usable.any():
        raise ValueError("no finite error; the envelope is undefined")
    order = np.where(usable, errors, np.inf)
    winner = int(np.argmin(order))
    return degrees[winner], float(errors[winner])
