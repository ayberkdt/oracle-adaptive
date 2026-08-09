"""``A-sign``: budgeted coordinate descent on the coupled objective.

The attainable schedule of the ladder's top rung, and the only solver here
that sees the off-diagonal coupling.  Holding every other interval fixed, the
terms involving :math:`N_q` are

.. math::
    \\underbrace{\\sum_{i,k\\in I_q}\\mathbf u_i(N)^{\\top}\\mathbf Q_{ik}
      \\mathbf u_k(N)}_{\\text{within block}}
    \\;+\\;2\\sum_{i\\in I_q}\\mathbf u_i(N)^{\\top}\\mathbf c_i^{(-q)}
    \\;+\\;\\lambda W_qN^2 ,

and the second group is the whole difference between this and a
sensitivity-weighted allocation: it rewards a defect pointing *against* the
accumulated displacement.

The bookkeeping that makes a sweep linear
-----------------------------------------
Excluding the block being updated collapses to two vectors.  For
:math:`i\\in I_q`, every :math:`k\\le i` outside the block lies to the left of
it and every :math:`k>i` outside it lies to the right, so

.. math::
    \\mathbf c_i^{(-q)}=\\tfrac1T\\bigl[\\mathbf A_i\\,\\mathbf P_{\\text{left}}
    +\\mathbf R_{\\text{right}}\\bigr],

with :math:`\\mathbf P_{\\text{left}}=\\sum_{k<\\min I_q}\\mathbf u_k` and
:math:`\\mathbf R_{\\text{right}}=\\sum_{k>\\max I_q}\\mathbf A_k\\mathbf u_k`.
The second does not depend on :math:`i` at all.

Sweeping in temporal order, the prefix is one vector carried forward and the
suffix is only ever read to the right of the block being updated --- a region
the sweep has not touched --- so the array built once at the start of the
sweep stays valid throughout it.  Rewriting either would cost
:math:`\\mathcal O(MK_{\\mathrm{dec}})` per sweep instead of
:math:`\\mathcal O(M\\lvert\\mathcal N\\rvert)`.  The stale suffix is not an
approximation but the definition: it is what makes the update Gauss--Seidel in
the intervals and Jacobi across the sweep (``DECISIONS.md`` D155).

What is not claimed
-------------------
The objective is not convex in the integer variables, so this returns a good
feasible point and not an optimum.  How far it can be from one is bounded
separately by :mod:`tda.allocate.frankwolfe`, and the spread across starts is
reported rather than suppressed.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from tda.allocate.budget import (
    ScheduleSolution,
    schedule_work,
    solve_to_budget,
)
from tda.kernel import CouplingKernel

__all__ = [
    "DescentProblem",
    "monotonicity_report",
    "solve_descent",
    "sweep_once",
]

Arr = NDArray[np.float64]
IntArr = NDArray[np.int_]


@dataclass(frozen=True, slots=True)
class DescentProblem:
    """Everything the sweep needs, assembled once per orbit.

    Attributes
    ----------
    contributions:
        :math:`\\mathbf u_i(N)` for every cell and candidate, shape
        ``(M, P, 6)``.  Built once from the defect table; the sweep performs
        no field evaluation.
    kernel:
        The coupling structure.
    interval_of:
        :math:`g(i)`, shape ``(M,)``.
    time_weight:
        :math:`W_q`, shape ``(K,)``.
    candidate_degrees:
        The candidate axis, shape ``(P,)``.
    cell_slices:
        Start and stop index of each interval's contiguous cell block, shape
        ``(K, 2)``.  Contiguity is guaranteed by the grid construction, and
        relying on it is what keeps the inner loop free of fancy indexing.
    """

    contributions: Arr
    kernel: CouplingKernel
    interval_of: IntArr
    time_weight: Arr
    candidate_degrees: tuple[int, ...]
    cell_slices: IntArr

    @classmethod
    def from_table(cls, table, kernel: CouplingKernel,
                   interval_of: IntArr, time_weight: Arr) -> DescentProblem:
        """Build the contribution tensor from a defect table.

        Parameters
        ----------
        table:
            A :class:`tda.tables.DefectTable`.
        kernel:
            Built from the same arc.
        interval_of, time_weight:
            The decision grid.

        Raises
        ------
        ValueError
            If the intervals are not contiguous in cell order, which would
            mean the decision grid and the accumulation grid were not built
            together.
        """
        contributions = np.einsum(
            "iab,ipb,i->ipa", table.node_transport, np.asarray(table.defect),
            table.widths)
        interval_of = np.asarray(interval_of)
        if np.any(np.diff(interval_of) < 0):
            raise ValueError(
                "interval_of must be non-decreasing; the decision grid and "
                "the accumulation grid were not built together")
        n_intervals = int(interval_of[-1]) + 1
        starts = np.searchsorted(interval_of, np.arange(n_intervals), "left")
        stops = np.searchsorted(interval_of, np.arange(n_intervals), "right")
        if np.any(stops <= starts):
            raise ValueError("some decision interval contains no cell")
        return cls(
            contributions=contributions,
            kernel=kernel,
            interval_of=interval_of,
            time_weight=np.asarray(time_weight, dtype=float),
            candidate_degrees=tuple(table.schema.candidate_degrees),
            cell_slices=np.stack([starts, stops], axis=1),
        )

    @property
    def n_intervals(self) -> int:
        """:math:`K_{\\mathrm{dec}}`."""
        return int(self.cell_slices.shape[0])

    def gather(self, schedule: IntArr) -> Arr:
        """:math:`\\mathbf u` for a per-interval schedule, shape ``(M, 6)``."""
        index = {d: p for p, d in enumerate(self.candidate_degrees)}
        columns = np.array([index[int(n)] for n in schedule])[self.interval_of]
        return self.contributions[np.arange(self.contributions.shape[0]),
                                  columns]


def sweep_once(problem: DescentProblem, columns: IntArr,
           multiplier: float) -> tuple[IntArr, bool]:
    """Run one Gauss--Seidel pass over the intervals, in temporal order.

    Shared with :mod:`tda.allocate.rounding`, whose polish is exactly this
    sweep applied to a rounded relaxation -- using the same code is what makes
    "polished" mean the same thing in both places.

    Returns
    -------
    tuple of (ndarray, bool)
        The updated candidate columns and whether any interval moved.
    """
    kernel = problem.kernel
    contributions = problem.contributions
    suffix = kernel.suffix
    duration = kernel.duration
    degrees = np.asarray(problem.candidate_degrees, dtype=float)

    u = contributions[np.arange(contributions.shape[0]),
                      columns[problem.interval_of]]
    # Frozen for this sweep: sum_{k > i} A_k u_k, read only to the right of
    # the block being updated.
    au = np.einsum("iab,ib->ia", suffix, u)
    right = np.zeros_like(au)
    right[:-1] = np.cumsum(au[::-1], axis=0)[::-1][1:]

    prefix = np.zeros(6)
    changed = False
    for q in range(problem.n_intervals):
        lo, hi = problem.cell_slices[q]
        block = contributions[lo:hi]                       # (m, P, 6)
        blocks_a = suffix[lo:hi]                           # (m, 6, 6)
        r_right = right[hi - 1]

        # Within-block term, Eq. (within-block): one running prefix over the
        # block, evaluated for every candidate at once.
        a_u = np.einsum("mab,mpb->mpa", blocks_a, block)   # (m, P, 6)
        diagonal = np.einsum("mpa,mpa->p", block, a_u)
        running = np.cumsum(block, axis=0) - block         # sum_{p'<p}
        cross_in = 2.0 * np.einsum("mpa,mpa->p", running, a_u)
        within = (diagonal + cross_in) / duration

        # Between-block term: 2 u^T c^(-q) with c^(-q) = (A_i P + R)/T.
        between = 2.0 * (np.einsum("mpa,mab,b->p", block, blocks_a, prefix)
                         + np.einsum("mpa,a->p", block, r_right)) / duration

        penalty = multiplier * problem.time_weight[q] * degrees**2
        total = within + between + penalty
        cost = problem.time_weight[q] * degrees**2
        pick = int(np.lexsort((cost, total))[0])
        if pick != columns[q]:
            changed = True
            columns[q] = pick
        prefix = prefix + block[:, pick, :].sum(axis=0)
    return columns, changed


def solve_descent(problem: DescentProblem, budget: float,
                  starts: Sequence[IntArr], max_sweeps: int = 12,
                  polish: bool = True) -> ScheduleSolution:
    """Multi-start budgeted coordinate descent.

    Parameters
    ----------
    problem:
        Assembled by :meth:`DescentProblem.from_table`.
    budget:
        The work ceiling.
    starts:
        Starting schedules, each shape ``(K,)`` in degrees.  The campaign
        declares them in advance --- the comparators, the separable solution
        and a fixed set of random seeds --- and reports the spread.
    max_sweeps:
        Cap on inner sweeps per start.
    polish:
        Whether to run the hard-budget polish of
        :func:`tda.allocate.polish.polish_to_budget` on the continuation's
        answer.  On by default: the multiplier sweep reaches only the
        supported points of a discrete trade-off curve, and the polish is
        what searches under the constraint as written (D182).  ``False`` is
        the control that measures what it was worth.

    Returns
    -------
    ScheduleSolution
        ``objective`` is :math:`J`; ``spread`` holds every start's :math:`J`
        at the winning multiplier, relative to the best.

    Raises
    ------
    ValueError
        If no start is given, or a start names an untabulated degree.

    Notes
    -----
    The multi-start loop sits **inside** the multiplier search rather than
    outside it: each multiplier is solved from every start.  Warm-starting
    from the previous multiplier would make the returned schedule depend on
    the search path, which is precisely the failure the monotonicity check is
    meant to detect.

    Selection among starts is on :math:`J+\\lambda W`, never on :math:`J`
    alone (D131) --- that is the Lagrangian subproblem.  Selection *across*
    multipliers, and the polish that follows, work on :math:`J` under the
    ceiling, which is the actual problem.
    """
    if not starts:
        raise ValueError("at least one starting schedule is required")
    index = {d: p for p, d in enumerate(problem.candidate_degrees)}
    try:
        start_columns = [np.array([index[int(n)] for n in s], dtype=np.int64)
                         for s in starts]
    except KeyError as exc:
        raise ValueError(
            f"a starting schedule names degree {exc.args[0]}, which is not a "
            "tabulated candidate") from None

    degrees = np.asarray(problem.candidate_degrees)
    last_spread: list[float] = []

    def solve_at(multiplier: float) -> tuple[IntArr, float]:
        best_penalized = np.inf
        best_columns: IntArr | None = None
        objectives: list[float] = []
        for initial in start_columns:
            columns = initial.copy()
            for _ in range(max_sweeps):
                columns, changed = sweep_once(problem, columns, multiplier)
                if not changed:
                    break
            schedule = degrees[columns]
            value = problem.kernel.objective(problem.gather(schedule))
            work = schedule_work(schedule, problem.time_weight)
            objectives.append(value)
            penalized = value + multiplier * work
            if penalized < best_penalized:
                best_penalized, best_columns = penalized, columns
        last_spread.clear()
        floor = min(objectives)
        last_spread.extend(v / floor if floor > 0 else 1.0 for v in objectives)
        return degrees[best_columns], problem.kernel.objective(
            problem.gather(degrees[best_columns]))

    solution = solve_to_budget(solve_at, problem.time_weight, budget)
    if not polish:
        return ScheduleSolution(
            degrees=solution.degrees, objective=solution.objective,
            work=solution.work, budget=solution.budget,
            multiplier=solution.multiplier, feasible=solution.feasible,
            spread=tuple(last_spread), diagnostics=solution.diagnostics,
        )

    from tda.allocate.polish import polish_to_budget

    polished, report = polish_to_budget(problem, solution.degrees, budget)
    diagnostics = dict(solution.diagnostics)
    diagnostics.update({
        "polish_single_moves": float(report.single_moves),
        "polish_exchange_moves": float(report.exchange_moves),
        "polish_exchanges_rejected": float(report.exchanges_rejected),
        "polish_passes": float(report.passes),
    })
    return ScheduleSolution(
        degrees=polished,
        objective=report.objective_after,
        work=report.work_after,
        budget=solution.budget,
        multiplier=solution.multiplier,
        feasible=True,
        spread=tuple(last_spread),
        diagnostics=diagnostics,
    )


def monotonicity_report(problem: DescentProblem, starts: Sequence[IntArr],
                        multipliers: Arr, max_sweeps: int = 12
                        ) -> dict[str, float]:
    """Report whether the attained work falls with the multiplier.

    Monotonicity is a property of a *global* minimiser.  What the descent
    returns is a best-of-starts local solution, and nothing guarantees that
    the basin it lands in varies monotonically --- so the campaign sweeps a
    dense logarithmic grid before any run and verifies it, rather than
    assuming it (D33).  Where it fails, the bisection is abandoned for that
    orbit and the schedule is taken as the feasible grid point with the
    smallest objective; the orbits on which this happens are counted and
    reported.

    Parameters
    ----------
    problem, starts, max_sweeps:
        As for :func:`solve_descent`.
    multipliers:
        Ascending grid to sweep.

    Returns
    -------
    dict
        ``violations`` --- how many adjacent pairs rise; ``worst_rise`` ---
        the largest relative increase; ``monotone`` --- 1.0 if none.
    """
    index = {d: p for p, d in enumerate(problem.candidate_degrees)}
    degrees = np.asarray(problem.candidate_degrees)
    works = []
    for multiplier in np.asarray(multipliers, dtype=float):
        best_penalized, best = np.inf, None
        for start in starts:
            columns = np.array([index[int(n)] for n in start], dtype=np.int64)
            for _ in range(max_sweeps):
                columns, changed = sweep_once(problem, columns, float(multiplier))
                if not changed:
                    break
            schedule = degrees[columns]
            work = schedule_work(schedule, problem.time_weight)
            value = problem.kernel.objective(problem.gather(schedule))
            if value + multiplier * work < best_penalized:
                best_penalized, best = value + multiplier * work, work
        works.append(best)

    work_array = np.array(works, dtype=float)
    rises = np.diff(work_array)
    violations = int(np.sum(rises > 0.0))
    worst = float(np.max(rises / work_array[:-1])) if rises.size else 0.0
    return {
        "violations": float(violations),
        "worst_rise": max(worst, 0.0),
        "monotone": 1.0 if violations == 0 else 0.0,
    }
