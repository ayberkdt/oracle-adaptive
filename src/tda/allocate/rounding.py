"""``S-round``: a second solver, out of the relaxation the certificate builds.

Frank--Wolfe already produces the relaxed optimum :math:`\\theta^{\\star}`,
and the campaign currently uses it only for a lower bound.  Three cheap steps
turn it into a *schedule* as well, at no extra field evaluation:

1. per-interval :math:`\\argmax_N\\theta_{qN}`;
2. **sum-up rounding**, which chooses at each interval the candidate whose
   accumulated relaxed mass has run furthest ahead of what it has been given
   [Sager2012]_, and so tracks the relaxed solution's integral rather than its
   pointwise mode;
3. one polishing sweep of the coordinate descent.

If ``S-round`` beats the multi-start descent it becomes the reported schedule
and the descent becomes the control (``DECISIONS.md`` D132/T2).  Either way
the comparison is worth having: the certified gap conflates how far the
descent sits above the optimum with how far the relaxation sits below it, and
a second solver on the same relaxation separates the two.

Why sum-up rounding and not just the mode
-----------------------------------------
Taking the largest weight interval by interval discards the fact that the
relaxation may split mass deliberately --- spending part of an interval at a
high degree is its way of expressing a budget it cannot meet with either
neighbour.  Sum-up rounding preserves the *running total* of that split, which
is the quantity the objective integrates, and its integer trajectory is known
to stay within a bounded distance of the relaxed one for a control that enters
the state linearly [Sager2012]_.  That bound does not transfer verbatim here,
because the degree enters the defect nonlinearly; it is the motivation, not a
guarantee, and the guarantee is the measured objective.

References
----------
.. [Sager2012] S. Sager, H. G. Bock, M. Diehl, "The integer approximation
   error in mixed-integer optimal control", *Mathematical Programming* 133,
   2012.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from tda.allocate.budget import ScheduleSolution, schedule_work
from tda.allocate.descent import DescentProblem, sweep_once

__all__ = ["argmax_rounding", "round_and_polish", "sum_up_rounding"]

Arr = NDArray[np.float64]
IntArr = NDArray[np.int_]


def argmax_rounding(theta: Arr, candidate_degrees: tuple[int, ...]) -> IntArr:
    """Per-interval mode of the relaxed solution.

    Parameters
    ----------
    theta:
        Relaxed weights, shape ``(K, P)``.
    candidate_degrees:
        The candidate axis.

    Returns
    -------
    ndarray of int, shape (K,)
    """
    return np.asarray(candidate_degrees)[np.argmax(theta, axis=1)]


def sum_up_rounding(theta: Arr, candidate_degrees: tuple[int, ...]) -> IntArr:
    """Sum-up rounding of the relaxed weights [Sager2012]_.

    At each interval the candidate chosen is the one whose accumulated relaxed
    mass most exceeds the mass already awarded to it, so the integer schedule
    tracks the relaxed solution's running total rather than its mode.

    Parameters
    ----------
    theta:
        Relaxed weights, shape ``(K, P)``.  Rows are expected to sum to one;
        they are not renormalised, because a row that does not is a sign the
        relaxation did not solve and should surface rather than be repaired.
    candidate_degrees:
        The candidate axis.

    Returns
    -------
    ndarray of int, shape (K,)

    Raises
    ------
    ValueError
        If the shapes disagree or a row is not a probability vector to within
        a loose tolerance.
    """
    theta = np.asarray(theta, dtype=float)
    degrees = np.asarray(candidate_degrees)
    if theta.ndim != 2 or theta.shape[1] != degrees.size:
        raise ValueError(
            f"theta {theta.shape} must be (K, {degrees.size})")
    if not np.allclose(theta.sum(axis=1), 1.0, atol=1e-6):
        raise ValueError(
            "each row of theta must sum to one; a row that does not means the "
            "relaxation did not solve, which is not something to repair here")

    awarded = np.zeros(degrees.size)
    accumulated = np.zeros(degrees.size)
    picks = np.empty(theta.shape[0], dtype=np.int64)
    for q in range(theta.shape[0]):
        accumulated += theta[q]
        choice = int(np.argmax(accumulated - awarded))
        awarded[choice] += 1.0
        picks[q] = choice
    return degrees[picks]


def round_and_polish(problem: DescentProblem, theta: Arr, budget: float,
                     multiplier: float, max_sweeps: int = 4
                     ) -> ScheduleSolution:
    """Round the relaxed solution both ways, polish, and keep the best.

    Parameters
    ----------
    problem:
        The same data the descent uses.
    theta:
        Relaxed weights from :mod:`tda.allocate.frankwolfe`, shape ``(K, P)``.
    budget:
        The work ceiling.
    multiplier:
        The multiplier to polish at --- normally the one the descent settled
        on, so that the polish enforces the same budget rather than a new one.
    max_sweeps:
        Polishing sweeps per candidate rounding.

    Returns
    -------
    ScheduleSolution
        The best **feasible** rounding.  ``diagnostics`` records which of the
        two roundings won and whether either was infeasible before polishing,
        because "the mode was infeasible and sum-up was not" is a fact about
        the relaxation worth keeping.

    Raises
    ------
    ValueError
        If neither rounding is feasible after polishing.  Reported rather than
        repaired: a schedule forced under the ceiling by an ad-hoc descent
        would no longer be the relaxation's rounding, and quoting it as one
        would misattribute the result.
    """
    index = {d: p for p, d in enumerate(problem.candidate_degrees)}
    degrees = np.asarray(problem.candidate_degrees)

    outcomes: dict[str, tuple[IntArr, float, float]] = {}
    infeasible_before = 0.0
    for name, schedule in (("argmax", argmax_rounding(theta,
                                                      problem.candidate_degrees)),
                           ("sum_up", sum_up_rounding(theta,
                                                      problem.candidate_degrees))):
        if schedule_work(schedule, problem.time_weight) > budget:
            infeasible_before += 1.0
        columns = np.array([index[int(n)] for n in schedule], dtype=np.int64)
        for _ in range(max_sweeps):
            columns, changed = sweep_once(problem, columns, multiplier)
            if not changed:
                break
        polished = degrees[columns]
        work = schedule_work(polished, problem.time_weight)
        if work <= budget:
            value = problem.kernel.objective(problem.gather(polished))
            outcomes[name] = (polished, value, work)

    if not outcomes:
        raise ValueError(
            "neither rounding of the relaxed solution is feasible after "
            f"polishing at multiplier {multiplier:.6g}; no S-round schedule "
            "is reported for this orbit")

    winner = min(outcomes, key=lambda k: outcomes[k][1])
    schedule, value, work = outcomes[winner]
    return ScheduleSolution(
        degrees=schedule,
        objective=value,
        work=work,
        budget=budget,
        multiplier=multiplier,
        feasible=True,
        diagnostics={
            "sum_up_won": 1.0 if winner == "sum_up" else 0.0,
            "infeasible_before_polish": infeasible_before,
            "roundings_feasible": float(len(outcomes)),
        },
    )
