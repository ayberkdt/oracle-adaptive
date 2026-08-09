"""Hard-budget polish: improve a feasible schedule under the constraint itself.

The multiplier continuation of :mod:`tda.allocate.budget` searches a penalized
objective, and over a discrete candidate set that reaches only the *supported*
points of the trade-off curve.  A constrained optimum sitting in a non-convex
notch is unreachable from any multiplier, however finely it is swept
(``DECISIONS.md`` D182).

What reaches it is search under the constraint as written.  This module takes
the continuation's answer and improves it against :math:`J` subject to
:math:`\\sum_qW_qN_q^2\\le B`, with no penalty term anywhere: a move is
accepted when it lowers the objective and still fits, and for no other reason.

Two move classes, and the second is why the first is not enough
---------------------------------------------------------------
A single-interval move can be blocked by the budget even when the schedule is
far from optimal: raising the interval that most wants raising may need work
that only becomes available if some *other* interval comes down.  A descent
restricted to one coordinate at a time therefore stalls at a point that is not
a local optimum of the constrained problem, only of the unconstrained
neighbourhood intersected with feasibility.

The exchange step supplies the missing direction: pair a raise with a lower
so the net work fits.  Enumerating all pairs is :math:`\\mathcal O(K^2P^2)`,
which at campaign scale is out of the question, so candidates are proposed by
marginal rate --- the best objective reduction per unit of work released or
consumed --- and that proposal is a *heuristic*.  Its acceptance is not: the
pair is applied and the objective re-evaluated exactly before the move is
kept, because the two intervals' cross terms depend on each other and the sum
of two independently computed changes is not the change of applying both.

What this is not
----------------
It is not a claim of optimality.  It removes a class of failure the
continuation has by construction; how close the result then sits to the true
integer optimum is measured, not argued, by
:func:`tda.allocate.exhaustive.verify_schedule` wherever enumeration is
affordable.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from tda.allocate.budget import schedule_work
from tda.allocate.descent import DescentProblem

__all__ = ["PolishReport", "block_deltas", "polish_to_budget"]

Arr = NDArray[np.float64]
IntArr = NDArray[np.int_]

_EXCHANGE_WIDTH = 5
"""How many leading raises and lowers the exchange step pairs against each other.

Twenty-five exact re-evaluations per pass at worst.  One each is too few --- the
best-rated lower frequently fails to release enough work for the best-rated
raise, and the pass then ends on a single miss.
"""


def _leaders(score: Arr, mask: NDArray[np.bool_], count: int,
             width: int) -> list[tuple[int, int]]:
    """The ``count`` smallest-scoring ``(interval, candidate)`` pairs in ``mask``."""
    flat = np.flatnonzero(mask.ravel())
    if flat.size == 0:
        return []
    order = flat[np.argsort(score.ravel()[flat])][:count]
    return [divmod(int(i), width) for i in order]


@dataclass(frozen=True, slots=True)
class PolishReport:
    """What the polish did, so that "nothing" and "not run" stay distinct.

    Attributes
    ----------
    single_moves, exchange_moves:
        Accepted moves of each class.  An exchange count above zero is direct
        evidence that the single-coordinate neighbourhood was blocked by the
        budget rather than exhausted.
    exchanges_rejected:
        Proposals whose exact re-evaluation did not improve the objective.
        Reported because it measures how good the marginal-rate heuristic is;
        if it were near zero the proposal rule could be trusted without the
        exact check, and it is not.
    objective_before, objective_after, work_after:
        The numbers a manifest carries.
    passes:
        Sweeps over the interval set before nothing moved.
    """

    single_moves: int
    exchange_moves: int
    exchanges_rejected: int
    objective_before: float
    objective_after: float
    work_after: float
    passes: int

    @property
    def improved(self) -> bool:
        """Whether the polish found anything the continuation had missed."""
        return self.objective_after < self.objective_before


def block_deltas(problem: DescentProblem, columns: IntArr) -> tuple[Arr, Arr]:
    """Objective and work change of every single-interval move, shape ``(K, P)``.

    Parameters
    ----------
    problem:
        The instance.
    columns:
        Current candidate column per interval, shape ``(K,)``.

    Returns
    -------
    tuple of (ndarray, ndarray)
        ``(delta_objective, delta_work)``, both ``(K, P)`` and both zero on
        the currently held column.

    Notes
    -----
    Exact for one move in isolation.  The block score is the same quantity
    :func:`tda.allocate.descent.sweep_once` minimises with the multiplier set
    to zero, so the difference against the held column is the exact change in
    :math:`J` from moving that one interval and nothing else.  Applying two of
    them together is *not* the sum of the two, which is why the exchange step
    re-evaluates.
    """
    kernel = problem.kernel
    contributions = problem.contributions
    suffix = kernel.suffix
    duration = kernel.duration
    degrees = np.asarray(problem.candidate_degrees, dtype=float)
    rows = np.arange(problem.n_intervals)

    u = contributions[np.arange(contributions.shape[0]),
                      columns[problem.interval_of]]
    au = np.einsum("iab,ib->ia", suffix, u)
    right = np.zeros_like(au)
    right[:-1] = np.cumsum(au[::-1], axis=0)[::-1][1:]
    prefix_all = np.zeros((u.shape[0] + 1, 6))
    prefix_all[1:] = np.cumsum(u, axis=0)

    score = np.empty((problem.n_intervals, degrees.size))
    for q in range(problem.n_intervals):
        lo, hi = (int(v) for v in problem.cell_slices[q])
        block = contributions[lo:hi]
        blocks_a = suffix[lo:hi]
        prefix = prefix_all[lo]
        r_right = right[hi - 1]

        a_u = np.einsum("mab,mpb->mpa", blocks_a, block)
        diagonal = np.einsum("mpa,mpa->p", block, a_u)
        running = np.cumsum(block, axis=0) - block
        cross_in = 2.0 * np.einsum("mpa,mpa->p", running, a_u)
        between = 2.0 * (np.einsum("mpa,mab,b->p", block, blocks_a, prefix)
                         + np.einsum("mpa,a->p", block, r_right))
        score[q] = (diagonal + cross_in + between) / duration

    delta_objective = score - score[rows, columns][:, None]
    cost = problem.time_weight[:, None] * degrees[None, :] ** 2
    delta_work = cost - cost[rows, columns][:, None]
    return delta_objective, delta_work


def polish_to_budget(problem: DescentProblem, schedule: IntArr, budget: float,
                     max_passes: int = 20, exchange: bool = True
                     ) -> tuple[IntArr, PolishReport]:
    """Improve a feasible schedule against ``J`` under the ceiling.

    Parameters
    ----------
    problem:
        The instance.
    schedule:
        A feasible starting schedule, normally the continuation's answer.
    budget:
        The work ceiling.
    max_passes:
        Cap on sweeps.  Each accepted move strictly lowers :math:`J`, so the
        loop cannot cycle; the cap bounds the cost on a pathological instance.
    exchange:
        Whether to run the paired raise/lower step.  ``False`` is the control
        that shows what the pairing was worth.

    Returns
    -------
    tuple of (ndarray, PolishReport)

    Raises
    ------
    ValueError
        If the starting schedule overspends, or names an untabulated degree.
        The polish improves a feasible point; it cannot rescue an infeasible
        one, and pretending otherwise would let an overspending schedule be
        reported as polished.
    """
    index = {int(d): p for p, d in enumerate(problem.candidate_degrees)}
    try:
        columns = np.array([index[int(n)] for n in np.asarray(schedule)],
                           dtype=np.int64)
    except KeyError as exc:
        raise ValueError(
            f"the schedule names degree {exc.args[0]}, which is not a "
            "tabulated candidate") from None
    degrees = np.asarray(problem.candidate_degrees)
    work = schedule_work(degrees[columns], problem.time_weight)
    if work > budget * (1.0 + 1.0e-9):
        raise ValueError(
            f"the starting schedule spends {work:.6g} against a ceiling of "
            f"{budget:.6g}; the polish improves a feasible point")

    start_objective = problem.kernel.objective(problem.gather(degrees[columns]))
    objective = start_objective
    singles = exchanges = rejected = passes = 0

    for _ in range(max_passes):
        passes += 1
        moved = False
        delta_j, delta_w = block_deltas(problem, columns)
        room = budget - work

        # --- single moves: strictly better and affordable -------------------
        affordable = (delta_j < 0.0) & (delta_w <= room * (1.0 + 1.0e-12))
        if affordable.any():
            masked = np.where(affordable, delta_j, 0.0)
            flat = int(np.argmin(masked))
            q, column = divmod(flat, delta_j.shape[1])
            work += float(delta_w[q, column])
            objective += float(delta_j[q, column])
            columns[q] = column
            singles += 1
            moved = True

        # --- exchange: pay for a raise by lowering somewhere else -----------
        elif exchange:
            # Rank raises by objective bought per unit of work needed and
            # lowers by objective given up per unit released, then try the
            # leading few against each other. One pair is not enough: the
            # best-rated lower often does not release enough work for the
            # best-rated raise, and a single attempt per pass ends the polish
            # on the first miss.
            raises = (delta_j < 0.0) & (delta_w > 0.0)
            lowers = (delta_j >= 0.0) & (delta_w < 0.0)
            if raises.any() and lowers.any():
                gain = np.where(raises, delta_j / np.where(delta_w > 0.0,
                                                           delta_w, 1.0), 0.0)
                loss = np.where(lowers, delta_j / np.where(delta_w < 0.0,
                                                           -delta_w, 1.0),
                                np.inf)
                best_raise = _leaders(gain, raises, _EXCHANGE_WIDTH,
                                      delta_j.shape[1])
                best_lower = _leaders(loss, lowers, _EXCHANGE_WIDTH,
                                      delta_j.shape[1])
                for rq, rc in best_raise:
                    for lq, lc in best_lower:
                        if rq == lq:
                            continue
                        trial = columns.copy()
                        trial[rq], trial[lq] = rc, lc
                        trial_work = schedule_work(degrees[trial],
                                                   problem.time_weight)
                        if trial_work > budget * (1.0 + 1.0e-9):
                            continue
                        # Exact, because the two intervals' cross terms see
                        # each other and the independent deltas do not.
                        trial_objective = problem.kernel.objective(
                            problem.gather(degrees[trial]))
                        if trial_objective < objective:
                            columns, work, objective = (trial, trial_work,
                                                        trial_objective)
                            exchanges += 1
                            moved = True
                            break
                        rejected += 1
                    if moved:
                        break
        if not moved:
            break

    final = degrees[columns]
    return final, PolishReport(
        single_moves=singles,
        exchange_moves=exchanges,
        exchanges_rejected=rejected,
        objective_before=start_objective,
        objective_after=problem.kernel.objective(problem.gather(final)),
        work_after=schedule_work(final, problem.time_weight),
        passes=passes,
    )
