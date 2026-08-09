"""The budget contract, shared by every solver in the package.

One place, because getting it wrong in one solver and right in another would
make their comparison meaningless.

The constraint is a **ceiling**, :math:`\\sum_qW_qN_q^2\\le B`, and a solution
that leaves work unspent is admissible.  That is not a corner case here: the
objective is not monotone in degree, because the whole content of the signed
formulation is that contributions cancel, so raising a degree can shrink one
term of a cancelling pair and *increase* :math:`J` (``DECISIONS.md`` D142).

The multiplier is a search device, not an optimality condition
--------------------------------------------------------------
The degrees are integers, so this is an integer program and the
Karush--Kuhn--Tucker conditions are **not** necessary for its optimum.
Complementary slackness in particular fails outright.  Take two candidates
with :math:`(W,J)` of :math:`(5,1)` and :math:`(10,0)` under :math:`B=6`: the
constrained optimum is the first, at :math:`J^\\star=1`, leaving a unit of
budget unspent, and no multiplier makes :math:`\\lambda(B-W^\\star)` vanish.
Stopping a sweep on :math:`\\lambda(B-W)<\\epsilon` would reject the right
answer (D182).

The deeper limit is the same one [Shoham1988]_ describes for separable
allocation: minimising :math:`J+\\lambda W` is a weighted-sum scalarisation,
and over a discrete set that recovers only the *supported* efficient points.
A constrained optimum sitting in a non-convex notch of the trade-off curve is
unreachable from any :math:`\\lambda`.

So :func:`solve_to_budget` is a **Lagrangian continuation**: it sweeps the
multiplier to generate good feasible candidates cheaply, and returns the one
with the smallest objective among *every* feasible schedule the sweep
produced.  It is not claimed to be optimal.  What establishes how close it
gets is :func:`tda.allocate.exhaustive.verify_schedule`, which enumerates
where enumeration is affordable, and the hard-budget polish that each solver
applies to the continuation's answer.

Keeping the best rather than the last is not a refinement.  A sweep that
returns whichever schedule the final bracket happened to land on can return a
*worse* schedule when given more starting points --- measured on a pilot arc,
where a four-start descent came back above a single-start one because the two
bisections stopped at different multipliers.

Selection among several starts *at a fixed multiplier* is on the penalized
objective :math:`F_\\lambda=J+\\lambda W`, never on :math:`J` alone: that is
the Lagrangian subproblem and solving it on :math:`J` would prefer whichever
basin bought a lower objective by spending more (D131).  Selection *across*
multipliers is on :math:`J` among the feasible, which is a different question
-- which feasible point is best -- and has a different answer.

References
----------
.. [Everett1963] H. Everett III, "Generalized Lagrange multiplier method for
   solving problems of optimum allocation of resources", *Operations
   Research* 11(3), 1963 -- the multiplier sweep over a discrete set.
.. [Shoham1988] Y. Shoham and A. Gersho, "Efficient bit allocation for an
   arbitrary set of quantizers", *IEEE Trans. ASSP* 36(9), 1988 -- and the
   observation that over a discrete set the sweep recovers only supported
   efficient solutions, so a duality gap may remain.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

__all__ = [
    "BUDGET_BISECTION_STEPS",
    "ScheduleSolution",
    "expand_to_cells",
    "schedule_work",
    "solve_to_budget",
]

Arr = NDArray[np.float64]
IntArr = NDArray[np.int_]

BUDGET_BISECTION_STEPS = 64
"""Bisection depth on the multiplier bracket."""


def expand_to_cells(schedule: IntArr, interval_of: IntArr) -> IntArr:
    """Spread a per-interval schedule onto the accumulation cells.

    Parameters
    ----------
    schedule:
        Degree per decision interval, shape ``(K,)``.
    interval_of:
        :math:`g(i)`, shape ``(M,)``.

    Returns
    -------
    ndarray of int, shape (M,)
    """
    return np.asarray(schedule)[interval_of]


def schedule_work(schedule: IntArr, time_weight: Arr) -> float:
    """:math:`\\sum_qW_qN_q^2`, the quantity the ceiling bounds."""
    return float(time_weight @ np.asarray(schedule, dtype=float) ** 2)


@dataclass(frozen=True, slots=True)
class ScheduleSolution:
    """A schedule together with the evidence about how it met its budget.

    Attributes
    ----------
    degrees:
        Degree per decision interval, shape ``(K,)``.
    objective:
        :math:`J` at that schedule.
    work:
        Realized nominal work :math:`\\sum_qW_qN_q^2`.
    budget:
        The ceiling.
    multiplier:
        The multiplier whose subproblem produced the returned schedule.  Not
        :math:`\\lambda^\\star` and not an optimality certificate --- the
        problem is integer and has no such multiplier in general (D182) ---
        but a record of where in the continuation the winner was found.
        Exactly zero means the unconstrained answer already fitted, which is a
        legitimate outcome the campaign reports rather than calibrates away.
    feasible:
        Whether the work fits.  Always true on return; carried so that a
        serialized record answers the question without recomputing.
    spread:
        Objective across the multi-start set at the winning multiplier,
        relative to the best.  Reported rather than suppressed: the descent
        is not convex in the integers.
    diagnostics:
        Free-form measurements a solver wants on the record --- sweep counts,
        monotonicity verdicts, and so on.
    """

    degrees: IntArr
    objective: float
    work: float
    budget: float
    multiplier: float
    feasible: bool
    spread: tuple[float, ...] = ()
    diagnostics: dict[str, float] = field(default_factory=dict)

    @property
    def utilisation(self) -> float:
        """Fraction of the ceiling spent."""
        return self.work / self.budget

    @property
    def ceiling_binds(self) -> bool:
        """Whether the multiplier is positive, i.e. the ceiling is active."""
        return self.multiplier > 0.0


def solve_to_budget(
    solve_at: Callable[[float], tuple[IntArr, float]],
    time_weight: Arr,
    budget: float,
    max_multiplier: float = 1.0e12,
) -> ScheduleSolution:
    """Find the schedule that meets the ceiling, testing ``lambda = 0`` first.

    Parameters
    ----------
    solve_at:
        Given a multiplier, return ``(degrees, objective)`` for the schedule
        that minimizes :math:`J+\\lambda W`.  The callable owns whatever
        multi-start or sweep logic it needs; this function only handles the
        multiplier.
    time_weight:
        :math:`W_q`, shape ``(K,)``.
    budget:
        The ceiling :math:`B`.
    max_multiplier:
        Where to give up bracketing.  Reaching it means even the smallest
        candidate degree overruns the budget, which is a statement about the
        budget.

    Returns
    -------
    ScheduleSolution

    Raises
    ------
    ValueError
        If the budget is not positive, or no multiplier in the bracket
        produces a feasible schedule.

    Notes
    -----
    The bisection is a search accelerator and nothing more.  It narrows the
    bracket so that few subproblems are solved; the answer is the best
    feasible schedule *seen anywhere*, and it does not have to have come from
    the last bracket.  Over a discrete degree set the attained work is a step
    function of the multiplier and no root need exist; converging on one is
    how a calibration loop silently returns an infeasible schedule
    [Shoham1988]_.

    There is no complementary-slackness stopping test, and there must not be:
    the problem is integer, and its optimum can leave budget unspent at every
    positive multiplier (D182).

    Monotonicity of the attained work in :math:`\\lambda` is a property of a
    *global* minimizer [Everett1963]_.  A separable solve has one; a
    coordinate descent does not, and its caller checks for the violation
    rather than assuming it away (D33).
    """
    if budget <= 0.0:
        raise ValueError(f"budget must be positive, got {budget}")
    time_weight = np.asarray(time_weight, dtype=float)

    # Every feasible schedule the continuation produces is a candidate; the
    # smallest objective among them wins, whichever multiplier found it.
    best: tuple[IntArr, float, float, float] | None = None
    seen = 0.0

    def offer(schedule: IntArr, value: float, multiplier: float) -> float:
        """Record a schedule if it fits, and keep it if it is the best so far."""
        nonlocal best, seen
        work = schedule_work(schedule, time_weight)
        if work > budget:
            return work
        seen += 1.0
        if best is None or value < best[1]:
            best = (schedule, value, work, multiplier)
        return work

    degrees, objective = solve_at(0.0)
    if offer(degrees, objective, 0.0) <= budget:
        # The unconstrained answer already fits, and nothing a positive
        # multiplier returns can beat it: the penalty only ever moves the
        # subproblem away from the objective's own minimiser.
        schedule, value, work, _ = best
        return ScheduleSolution(
            degrees=schedule, objective=value, work=work, budget=budget,
            multiplier=0.0, feasible=True,
            diagnostics={"bracket_doublings": 0.0, "ceiling_binds": 0.0,
                         "feasible_seen": seen},
        )

    lo, hi, doublings = 0.0, 1.0, 0.0
    while hi <= max_multiplier:
        trial, value = solve_at(hi)
        if offer(trial, value, hi) <= budget:
            break
        lo, hi, doublings = hi, hi * 2.0, doublings + 1.0
    if best is None:
        raise ValueError(
            f"no multiplier up to {max_multiplier:.3g} fits the budget "
            f"{budget:.6g}; the smallest schedule available costs "
            f"{schedule_work(solve_at(max_multiplier)[0], time_weight):.6g}")

    for _ in range(BUDGET_BISECTION_STEPS):
        mid = 0.5 * (lo + hi)
        trial, value = solve_at(mid)
        if offer(trial, value, mid) <= budget:
            hi = mid
        else:
            lo = mid

    degrees, objective, work, multiplier = best
    return ScheduleSolution(
        degrees=degrees, objective=objective, work=work, budget=budget,
        multiplier=multiplier, feasible=True,
        diagnostics={"bracket_doublings": doublings, "ceiling_binds": 1.0,
                     "feasible_seen": seen},
    )
