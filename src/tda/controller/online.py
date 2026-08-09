"""The online half: one decision per boundary, from a probe taken ahead of it.

Algorithm 3 of the manuscript.  At the boundary :math:`t_q` the controller
predicts where it is about to be, evaluates one shared band stack there, forms
:math:`\\hat{\\mathbf v}` for every candidate in a window around the plan's
nominal degree by partial summation, scores them against the plan's
cancellation target, and holds the winner across the interval.

Nothing from before the boundary is reused for the direction.  A rule that
accumulated probes *inside* an interval and applied them at the *next* boundary
would select from a direction that had already decorrelated over
:math:`\\pi r/N`, which is the failure the resolution-scale argument predicts.

The score is the coordinate update, exactly
-------------------------------------------
For the block :math:`I_q`, with :math:`\\mathbf u_i(N)=\\mathbf S_i
\\hat{\\mathbf v}_i(N)`,

.. math::
    F(N)=\\underbrace{\\frac1T\\sum_{i,k\\in I_q}\\mathbf u_i^{\\top}
      \\mathbf A_{\\max(i,k)}\\mathbf u_k}_{\\text{within block}}
    \\;+\\;2\\sum_{i\\in I_q}\\hat{\\mathbf v}_i(N)^{\\top}\\mathbf z_i
    \\;+\\;\\lambda W_qN^2 ,

which differs from :math:`J` of the whole schedule by a constant that does not
depend on :math:`N`.  The second term is where the plan enters and it is the
only place the offline computation is used online; the first is formed from
data the plan carries but is evaluated against the probe's own estimate.

Two variants, one code path
---------------------------
``C-plan`` places :math:`n_{\\mathrm{probe}}` points at the midpoints of a
uniform split of the coming interval, each predicted analytically from the
boundary state.  ``C-lite`` places a single point at the boundary itself, where
the integrator is about to evaluate the field anyway, so the probe costs the
increment :math:`2(\\Delta_{\\mathrm{span}}+k)/N` rather than a synthesis --- and
the direction is known only at the interval's start.  They differ in
:attr:`OnlineSettings.placement` and :attr:`OnlineSettings.n_probe` and in
nothing else, so the measured difference between them is a difference in
probing and not in implementation (``DECISIONS.md`` D172).

Matching the probe to the plan
------------------------------
The phase index says which *planned* interval the vehicle is in; within it, a
probe point taken at fraction :math:`f` of the flown interval is matched to the
plan cells nearest :math:`f` of the planned one.  Comparing fractions rather
than absolute phase keeps the match insensitive to the two arcs' intervals
having slightly different phase widths, and the nearest-point rule is the
piecewise-constant reconstruction the correlation-time argument already
assumes: the direction is coherent over :math:`\\tau_{\\mathrm{corr}}` and the
probe samples once per :math:`\\tau_{\\mathrm{corr}}`, so there is nothing
finer to reconstruct.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Literal

import numpy as np
from numpy.typing import NDArray

from tda.controller.feedback import BudgetFeedback
from tda.controller.plan import OfflinePlan
from tda.probe import (
    ProbePlan,
    band_direction,
    candidate_window,
    plan_probe_points,
    required_degrees,
)
from tda.spectrum import BandStack

__all__ = [
    "PLACEMENTS",
    "Decision",
    "OnlineController",
    "OnlineSettings",
    "assign_probe_points",
    "colocated_cost_fraction",
    "probe_fractions",
    "score_candidates",
    "select_candidate",
]

Arr = NDArray[np.float64]
IntArr = NDArray[np.int_]

Placement = Literal["midpoint", "leading"]

PLACEMENTS = ("midpoint", "leading")
"""Where the probe points sit inside the interval about to be committed to."""


# ---------------------------------------------------------------------------
# Probe geometry inside an interval
# ---------------------------------------------------------------------------


def probe_fractions(n_points: int, placement: Placement) -> Arr:
    """Positions of the probe points within the interval, in :math:`[0,1)`.

    Parameters
    ----------
    n_points:
        How many points.
    placement:
        ``"midpoint"`` for ``C-plan``, ``"leading"`` for ``C-lite``.

    Returns
    -------
    ndarray, shape (n_points,)

    Raises
    ------
    ValueError
        If the count is not positive, the placement is not declared, or
        ``"leading"`` is asked for more than one point --- a co-located probe
        is co-located with *one* call, and a second leading point would be a
        midpoint probe under the cheaper variant's name.

    Examples
    --------
    >>> probe_fractions(4, "midpoint")
    array([0.125, 0.375, 0.625, 0.875])
    >>> probe_fractions(1, "leading")
    array([0.])
    """
    if placement not in PLACEMENTS:
        raise ValueError(
            f"placement must be one of {PLACEMENTS}, got {placement!r}")
    if n_points <= 0:
        raise ValueError(f"n_points must be positive, got {n_points}")
    if placement == "leading":
        if n_points != 1:
            raise ValueError(
                "the leading placement takes exactly one point; "
                f"got {n_points}. More than one co-located probe is a "
                "midpoint probe under C-lite's name and would be priced with "
                "C-lite's cost model")
        return np.zeros(1, dtype=float)
    return (np.arange(n_points, dtype=float) + 0.5) / n_points


def assign_probe_points(cell_fractions: Arr, point_fractions: Arr) -> IntArr:
    """Nearest probe point to each plan cell, by position in the interval.

    Parameters
    ----------
    cell_fractions:
        Where each plan cell sits, shape ``(m,)``.
    point_fractions:
        Where each probe point sits, ascending, shape ``(n,)``.

    Returns
    -------
    ndarray of int, shape (m,)

    Raises
    ------
    ValueError
        If either array is empty or the probe fractions are not ascending.

    Examples
    --------
    >>> import numpy as np
    >>> assign_probe_points(np.array([0.1, 0.4, 0.9]), np.array([0.25, 0.75]))
    array([0, 0, 1])
    """
    cell_fractions = np.asarray(cell_fractions, dtype=float)
    point_fractions = np.asarray(point_fractions, dtype=float)
    if cell_fractions.size == 0 or point_fractions.size == 0:
        raise ValueError("both cells and probe points are required")
    if np.any(np.diff(point_fractions) < 0.0):
        raise ValueError("point_fractions must be ascending")
    if point_fractions.size == 1:
        # The single-point case has to be taken out before the bracket search:
        # clipping an index into [1, 0] is an empty range, and the fallback
        # would hand back -1 for any cell at or before the point. It would
        # still address the right entry here, being the only one, and would be
        # wrong the moment anything read the index itself.
        return np.zeros(cell_fractions.shape, dtype=np.int64)
    idx = np.clip(np.searchsorted(point_fractions, cell_fractions), 1,
                  point_fractions.size - 1)
    take_lower = (cell_fractions - point_fractions[idx - 1]
                  <= point_fractions[idx] - cell_fractions)
    return np.where(take_lower, idx - 1, idx).astype(np.int64)


def colocated_cost_fraction(degree: int, span: int, depth: int) -> float:
    """Cost of a co-located probe, as a share of the synthesis it rides on.

    :math:`2(\\Delta_{\\mathrm{span}}+k)/N`, from the manuscript's incremental
    estimate :math:`2N(\\Delta_{\\mathrm{span}}+k)` against :math:`N^2`.  Valid
    **only** where the field is already being evaluated at that point: at a
    look-ahead point the associated Legendre recursion must run from degree
    zero and the probe costs a full synthesis, not an increment.

    Parameters
    ----------
    degree:
        :math:`N`, positive.
    span:
        :math:`\\Delta_{\\mathrm{span}}=N_{j+\\delta}-N_{j-\\delta}`, the
        window's width in degrees --- not its width in candidate indices.
    depth:
        Probe depth :math:`k`, positive.

    Returns
    -------
    float

    Raises
    ------
    ValueError
        If the degree or depth is not positive, or the span is negative.

    Examples
    --------
    >>> round(colocated_cost_fraction(120, 0, 3), 4)
    0.05
    """
    if degree <= 0:
        raise ValueError(f"degree must be positive, got {degree}")
    if depth <= 0:
        raise ValueError(f"depth must be positive, got {depth}")
    if span < 0:
        raise ValueError(f"span must be non-negative, got {span}")
    return 2.0 * (span + depth) / degree


# ---------------------------------------------------------------------------
# The decision
# ---------------------------------------------------------------------------


def score_candidates(plan: OfflinePlan, interval: int, estimates: Arr,
                     multiplier: float, window: tuple[int, ...]) -> Arr:
    """Penalised score of every candidate for one interval.

    Parameters
    ----------
    plan:
        Supplies :math:`\\mathbf S_i`, :math:`\\mathbf A_i`, :math:`\\mathbf
        z_i`, :math:`T` and :math:`W_q`.
    interval:
        :math:`q`.
    estimates:
        :math:`\\hat{\\mathbf v}_i(N)` for every cell of the block and every
        candidate, shape ``(m, P, 3)``.
    multiplier:
        :math:`\\lambda`, from the plan or from the budget feedback.
    window:
        The candidate degrees, length ``P``, matching ``estimates``.

    Returns
    -------
    ndarray, shape (P,)
        Differs from :math:`J` of the corresponding full schedule by a
        constant independent of the candidate, plus the work penalty.

    Raises
    ------
    ValueError
        If the estimate tensor does not match the block or the window.

    Notes
    -----
    The within-block quadratic is accumulated with one running prefix, which is
    the same identity :func:`tda.allocate.descent.sweep_once` uses:
    :math:`\\sum_{i,k}\\mathbf u_i^{\\top}\\mathbf A_{\\max(i,k)}\\mathbf u_k
    =\\sum_i(\\mathbf u_i+2\\sum_{k<i}\\mathbf u_k)^{\\top}\\mathbf A_i
    \\mathbf u_i`.
    """
    lo, hi = (int(v) for v in plan.cell_slices[interval])
    estimates = np.asarray(estimates, dtype=float)
    if estimates.shape != (hi - lo, len(window), 3):
        raise ValueError(
            f"estimates must have shape ({hi - lo}, {len(window)}, 3), got "
            f"{estimates.shape}")

    transport = plan.transport[lo:hi]
    suffix = plan.kernel.suffix[lo:hi]
    target = plan.target[lo:hi]

    u = np.einsum("mab,mpb->mpa", transport, estimates)
    au = np.einsum("mab,mpb->mpa", suffix, u)
    diagonal = np.einsum("mpa,mpa->p", u, au)
    running = np.cumsum(u, axis=0) - u
    cross_in = 2.0 * np.einsum("mpa,mpa->p", running, au)
    within = (diagonal + cross_in) / plan.kernel.duration

    between = 2.0 * np.einsum("mpa,ma->p", estimates, target)
    degrees = np.asarray(window, dtype=float)
    penalty = multiplier * plan.time_weight[interval] * degrees**2
    return within + between + penalty


def select_candidate(scores: Arr, window: tuple[int, ...]) -> int:
    """Pick the winning degree, breaking ties toward the cheaper one.

    The same lexicographic rule as the offline sweep.  Sharing it is what makes
    a difference between the controller and the benchmark a difference in
    information rather than in tie-breaking.

    Examples
    --------
    >>> import numpy as np
    >>> select_candidate(np.array([1.0, 1.0, 2.0]), (30, 40, 50))
    30
    """
    scores = np.asarray(scores, dtype=float)
    cost = np.asarray(window, dtype=float) ** 2
    return int(window[int(np.lexsort((cost, scores))[0])])


@dataclass(frozen=True, slots=True)
class Decision:
    """One boundary's outcome, with the evidence for what it cost.

    Attributes
    ----------
    interval:
        The plan interval the vehicle was matched to.
    degree:
        :math:`N_q`, held across the coming interval.
    window:
        The candidates searched.
    scores:
        Their penalised scores, aligned with ``window``.
    multiplier:
        The :math:`\\lambda` actually used, after any budget feedback.
    syntheses:
        Field evaluations the band stack performed for this decision,
        *measured* from the stack rather than predicted.  With no
        cumulative-by-degree kernel entry point this is one per requested
        truncation (D120), and the accounting reflects the code that ran.
    probe_points:
        :math:`n_{\\mathrm{probe}}`.
    co_located:
        Whether the probe rode on a call the integrator makes anyway, in which
        case :func:`colocated_cost_fraction` prices it and ``syntheses``
        overstates it.
    span_degrees:
        :math:`\\Delta_{\\mathrm{span}}`, the window's width in degrees --- the
        quantity that sets the band stack, as against the index half-width that
        sets the candidate count.
    predictor_error_bound:
        The two-body predictor's analytic bound, metres; exactly zero for the
        co-located variant, which uses no predictor at all.
    """

    interval: int
    degree: int
    nominal: int
    window: tuple[int, ...]
    scores: tuple[float, ...]
    multiplier: float
    syntheses: int
    probe_points: int
    co_located: bool
    span_degrees: int
    predictor_error_bound: float

    @property
    def moved(self) -> bool:
        """Whether the online search left the plan's nominal degree.

        Compared against the nominal rather than against the middle of the
        window: near the ends of the candidate grid the window is clipped and
        its middle element is not the plan's degree, so the obvious test would
        report a move where none happened.
        """
        return self.degree != self.nominal


@dataclass(frozen=True, slots=True)
class OnlineSettings:
    """Declared parameters of the online decision.

    Attributes
    ----------
    half_width:
        :math:`\\delta`, in candidate-grid *indices*.
    depth:
        Probe depth :math:`k`.
    n_probe:
        Probe points per interval.
    placement:
        ``"midpoint"`` (``C-plan``) or ``"leading"`` (``C-lite``).
    """

    half_width: int = 2
    depth: int = 3
    n_probe: int = 4
    placement: Placement = "midpoint"

    def __post_init__(self) -> None:
        """Reject a combination the cost model would misprice."""
        # Validated here rather than at the first decision: an inconsistent
        # setting must fail when the run is configured, not a thousand
        # boundaries in.
        probe_fractions(self.n_probe, self.placement)
        if self.half_width < 0:
            raise ValueError(
                f"half_width must be non-negative, got {self.half_width}")
        if self.depth <= 0:
            raise ValueError(f"depth must be positive, got {self.depth}")


@dataclass(slots=True)
class OnlineController:
    """``C-plan`` / ``C-lite``, optionally closed on the realized budget.

    Attributes
    ----------
    plan:
        The phase-indexed plan.
    settings:
        Declared online parameters.
    feedback:
        ``C-fb`` if present; without it the multiplier stays at the plan's
        :math:`\\lambda_0` and the controller is open loop, which is the
        control the campaign measures the loop against.
    decisions:
        How many boundaries have been decided.
    syntheses:
        Total probe syntheses, measured.
    """

    plan: OfflinePlan
    settings: OnlineSettings = field(default_factory=OnlineSettings)
    feedback: BudgetFeedback | None = None
    decisions: int = 0
    syntheses: int = 0

    def charge(self, degree: int, n_calls: int = 1) -> None:
        """Record realized gravity work, if the budget loop is closed.

        Raises
        ------
        RuntimeError
            If there is no feedback loop.  Silently discarding the charge
            would leave an open-loop run looking like a closed one whose error
            happened to stay at zero.
        """
        if self.feedback is None:
            raise RuntimeError(
                "this controller is open loop; attach a BudgetFeedback to "
                "account realized work")
        self.feedback.charge(degree, n_calls)

    def multiplier_at(self, coordinate: float) -> float:
        """:math:`\\lambda` for the interval starting at this phase."""
        if self.feedback is None:
            return self.plan.multiplier
        return self.feedback.multiplier_at(coordinate)

    def decide(self, coordinate: float, state: Arr, epoch: float, span: float,
               stack: BandStack, sigma_a_of: Callable[[float], Arr],
               mu: float, perturbing_acceleration: float = 0.0) -> Decision:
        """Choose the degree for the interval starting at ``epoch``.

        Parameters
        ----------
        coordinate:
            The vehicle's own ``revolution + phase``, from its
            :class:`~tda.controller.phase.RevolutionIndex`.
        state:
            Current inertial state ``[r, v]``, shape ``(6,)``.
        epoch:
            :math:`t_q`, seconds from arc start.
        span:
            :math:`\\Delta t_{\\mathrm{dec}}` of the coming interval, on the
            flown arc.
        stack:
            The band stack; its synthesis counter is read before and after.
        sigma_a_of:
            Radius to per-degree acceleration RMS, from
            :func:`tda.spectrum.acceleration_degree_rms`.  Passed as a callable
            because the completion factor depends on where the probe point is,
            and evaluating it at the boundary radius for a point half an
            interval away would misweight the tail by the attenuation over
            that arc.
        mu:
            Gravitational parameter, for the analytic predictor.
        perturbing_acceleration:
            Magnitude used for the predictor's error bound.

        Returns
        -------
        Decision

        Raises
        ------
        IndexError
            If the phase coordinate is outside the plan's coverage.
        """
        settings = self.settings
        interval = self.plan.locate(coordinate)
        window = candidate_window(self.plan.candidate_degrees,
                                  int(self.plan.nominal_degrees[interval]),
                                  settings.half_width)
        fractions = probe_fractions(settings.n_probe, settings.placement)

        if settings.placement == "leading":
            probe = ProbePlan(times=np.array([epoch]),
                              positions=np.asarray(state, dtype=float)[None,
                                                                       0:3],
                              predictor_error_bound=0.0)
        else:
            probe = plan_probe_points(state, mu, epoch, span,
                                      settings.n_probe,
                                      perturbing_acceleration)

        before = stack.total_syntheses
        per_point = [
            band_direction(stack, position, float(time), window,
                           settings.depth,
                           sigma_a_of(float(np.linalg.norm(position))))
            for position, time in zip(probe.positions, probe.times,
                                      strict=True)
        ]
        spent = stack.total_syntheses - before

        assignment = assign_probe_points(self.plan.cell_fractions(interval),
                                         fractions)
        estimates = np.stack([
            np.stack([per_point[point][degree] for degree in window])
            for point in assignment
        ])

        multiplier = self.multiplier_at(coordinate)
        scores = score_candidates(self.plan, interval, estimates, multiplier,
                                  window)
        degree = select_candidate(scores, window)

        self.decisions += 1
        self.syntheses += spent
        return Decision(
            interval=interval,
            degree=degree,
            nominal=int(self.plan.nominal_degrees[interval]),
            window=window,
            scores=tuple(float(v) for v in scores),
            multiplier=multiplier,
            syntheses=int(spent),
            probe_points=len(probe),
            co_located=settings.placement == "leading",
            span_degrees=int(window[-1] - window[0]),
            predictor_error_bound=probe.predictor_error_bound,
        )

    def stack_depth(self, interval: int) -> int:
        """Distinct truncations the shared stack must supply for an interval.

        The union over the window, not the sum over candidates: that the two
        differ is the whole content of "the candidates are partial sums of one
        stack" (``paper`` §7.2), and a controller that reported the sum would
        be describing a design nobody proposed.
        """
        window = candidate_window(self.plan.candidate_degrees,
                                  int(self.plan.nominal_degrees[interval]),
                                  self.settings.half_width)
        return int(required_degrees(window, self.settings.depth).size)
