"""The budget contract, shared by every solver in the package.

One place, because getting it wrong in one solver and right in another would
make their comparison meaningless.

The constraint is a **ceiling**, :math:`\\sum_qW_qN_q^2\\le B`, and the
Karush--Kuhn--Tucker conditions for it are :math:`\\lambda\\ge0`,
:math:`B-W\\ge0` and :math:`\\lambda(B-W)=0`.  A solution that leaves work
unspent is admissible whenever :math:`\\lambda^\\star=0`, and that is not a
corner case here: the objective is not monotone in degree, because the whole
content of the signed formulation is that contributions cancel, so raising a
degree can shrink one term of a cancelling pair and *increase* :math:`J`.
:func:`solve_to_budget` therefore tests :math:`\\lambda=0` before it brackets
anything (``DECISIONS.md`` D142).

Selection among several starts is on the *penalized* objective
:math:`F_\\lambda=J+\\lambda W`, never on :math:`J` alone.  At a fixed
multiplier different starts land in basins that do not spend the same work, so
selecting on :math:`J` would prefer whichever basin bought a lower objective
by spending more --- which is not the solution of the Lagrangian subproblem,
and leaves the attained work free to rise with :math:`\\lambda` for reasons
that have nothing to do with the problem (D131).

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
        :math:`\\lambda^\\star`.  Exactly zero means the ceiling did not bind
        and the schedule is the unconstrained optimum of its own criterion --
        a legitimate outcome, and one the campaign reports rather than
        calibrates away.
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
    The bisection tracks the best *feasible* iterate rather than converging on
    a root.  Over a discrete degree set the attained work is a step function
    of the multiplier and the root need not exist; chasing one is how a
    calibration loop silently returns an infeasible schedule
    [Shoham1988]_.

    Monotonicity of the attained work in :math:`\\lambda` is a property of a
    *global* minimizer [Everett1963]_.  A separable solve has one; a
    coordinate descent does not, and its caller checks for the violation
    rather than assuming it away (D33).
    """
    if budget <= 0.0:
        raise ValueError(f"budget must be positive, got {budget}")
    time_weight = np.asarray(time_weight, dtype=float)

    degrees, objective = solve_at(0.0)
    work = schedule_work(degrees, time_weight)
    if work <= budget:
        return ScheduleSolution(
            degrees=degrees, objective=objective, work=work, budget=budget,
            multiplier=0.0, feasible=True,
            diagnostics={"bracket_doublings": 0.0, "ceiling_binds": 0.0},
        )

    lo, hi, doublings = 0.0, 1.0, 0.0
    best: tuple[IntArr, float, float] | None = None
    while hi <= max_multiplier:
        trial, value = solve_at(hi)
        trial_work = schedule_work(trial, time_weight)
        if trial_work <= budget:
            best = (trial, value, trial_work)
            break
        lo, hi, doublings = hi, hi * 2.0, doublings + 1.0
    if best is None:
        raise ValueError(
            f"no multiplier up to {max_multiplier:.3g} fits the budget "
            f"{budget:.6g}; the smallest schedule available costs "
            f"{schedule_work(solve_at(max_multiplier)[0], time_weight):.6g}")

    best_multiplier = hi
    for _ in range(BUDGET_BISECTION_STEPS):
        mid = 0.5 * (lo + hi)
        trial, value = solve_at(mid)
        trial_work = schedule_work(trial, time_weight)
        if trial_work <= budget:
            best, best_multiplier, hi = (trial, value, trial_work), mid, mid
        else:
            lo = mid

    degrees, objective, work = best
    return ScheduleSolution(
        degrees=degrees, objective=objective, work=work, budget=budget,
        multiplier=best_multiplier, feasible=True,
        diagnostics={"bracket_doublings": doublings, "ceiling_binds": 1.0},
    )
