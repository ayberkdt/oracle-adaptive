"""Exact optimisation by enumeration, where enumeration is affordable.

The campaign's benchmark is only a benchmark if something says how far it is
from the best schedule available.  That was to be the Frank--Wolfe certificate
of :mod:`tda.allocate.frankwolfe`, and on the pilot arc the certificate turned
out to be bounded not by the solver but by the relaxation: with the relaxed
problem solved to convergence the remaining integrality gap was 1.00, 1.00,
1.00, 2.13 and 1.32 at four, six, eight, ten and twelve decision intervals.
It does not shrink with the interval count, it is not smooth, and there is no
basis for extrapolating it to the eight hundred intervals of a real arc
(``DECISIONS.md`` D178).

What the same instances *did* show is that the declared coordinate descent
returned the exact integer optimum in every one of them.  That is a stronger
statement than any bound would have made --- a certificate at
:math:`g_E=0.13` says "within thirteen per cent in error"; exhaustion says
"optimal" --- and it is obtained by comparing against truth rather than
against a relaxation.

So the verification of the solver is a **panel**: enumerate every feasible
schedule on instances small enough to allow it, over several arcs and several
ceilings, and report how often the solver is exactly optimal and how far off
it is when it is not.  The certificate stays, as a secondary diagnostic on
the full-size problem where enumeration is impossible.

What this does and does not establish
-------------------------------------
It establishes optimality *for the instances enumerated*.  It says nothing
about :math:`K_{\\mathrm{dec}}=810` beyond the ordinary inductive weight of a
solver that has never been caught suboptimal, and the manuscript must say so
in those words.  A panel that reported "optimal" and let the reader supply
the extrapolation would be worse than no panel.

Cost
----
:math:`\\lvert\\mathcal N\\rvert^{K}` schedules, so the panel lives at
:math:`K\\le12` with three candidates or :math:`K\\le9` with four.  The
enumeration refuses rather than runs when the count exceeds its cap: a
verification that silently became a four-hour job is a verification nobody
runs.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from tda.allocate.budget import ScheduleSolution, schedule_work
from tda.allocate.descent import DescentProblem

__all__ = [
    "MAX_SCHEDULES",
    "VerificationRecord",
    "panel_summary",
    "schedule_count",
    "solve_exhaustive",
    "verify_schedule",
]

Arr = NDArray[np.float64]
IntArr = NDArray[np.int_]

MAX_SCHEDULES = 2_000_000
"""Enumeration cap.  Above this the panel refuses rather than runs."""

_REFINE = 64
"""How many leaders are re-scored with the compensated objective."""


def schedule_count(problem: DescentProblem) -> int:
    """:math:`\\lvert\\mathcal N\\rvert^{K}`, the size of the enumeration."""
    return len(problem.candidate_degrees) ** problem.n_intervals


def _all_columns(problem: DescentProblem) -> IntArr:
    """Every candidate assignment, shape ``(P**K, K)``, in odometer order."""
    n_candidates = len(problem.candidate_degrees)
    n_intervals = problem.n_intervals
    codes = np.arange(n_candidates**n_intervals, dtype=np.int64)
    places = n_candidates ** np.arange(n_intervals - 1, -1, -1, dtype=np.int64)
    return (codes[:, None] // places[None, :]) % n_candidates


def solve_exhaustive(problem: DescentProblem, budget: float,
                     max_schedules: int = MAX_SCHEDULES,
                     chunk: int = 256) -> ScheduleSolution:
    """The true integer optimum, by enumerating every feasible schedule.

    Parameters
    ----------
    problem:
        The instance.
    budget:
        The work ceiling.
    max_schedules:
        Refuse if the enumeration is larger than this.
    chunk:
        Schedules evaluated per batch.  Sets the peak memory, which is
        ``chunk * M * 6`` floats.

    Returns
    -------
    ScheduleSolution
        ``diagnostics`` carries ``n_total`` and ``n_feasible``.

    Raises
    ------
    ValueError
        If the enumeration exceeds ``max_schedules``, or no schedule fits the
        ceiling.

    Notes
    -----
    Ranking runs on a plain cumulative sum so that a whole batch goes through
    one einsum; the leaders are then re-scored with
    :meth:`~tda.kernel.CouplingKernel.objective`, whose compensated prefix sum
    is what the rest of the campaign uses.  Keeping the top
    :data:`_REFINE` rather than only the batch leader is what makes the two
    stages agree: with a cancellation ratio of order ten the plain sum can
    reorder schedules that are within a few ulp of each other, and taking one
    winner per batch would let that reordering decide the answer.
    """
    total = schedule_count(problem)
    if total > max_schedules:
        raise ValueError(
            f"enumeration would visit {total:,} schedules "
            f"({len(problem.candidate_degrees)} candidates over "
            f"{problem.n_intervals} intervals), above the cap of "
            f"{max_schedules:,}; verify on a coarser decision grid or a "
            "smaller candidate set rather than raising the cap")

    kernel = problem.kernel
    degrees = np.asarray(problem.candidate_degrees, dtype=float)
    cost = problem.time_weight[:, None] * degrees[None, :] ** 2
    columns = _all_columns(problem)
    rows = np.arange(problem.n_intervals)
    work = cost[rows[None, :], columns].sum(axis=1)
    columns = columns[work <= budget]
    if columns.shape[0] == 0:
        raise ValueError(
            f"no schedule fits a ceiling of {budget:.6g}; the cheapest costs "
            f"{float(work.min()):.6g}")

    n_cells = problem.contributions.shape[0]
    cells = np.arange(n_cells)
    leaders: list[tuple[float, IntArr]] = []
    for start in range(0, columns.shape[0], chunk):
        block = columns[start:start + chunk]
        u = problem.contributions[cells[None, :], block[:, problem.interval_of]]
        prefix = np.zeros((u.shape[0], n_cells + 1, 6))
        np.cumsum(u, axis=1, out=prefix[:, 1:])
        displacement = np.einsum("jai,bji->bja", kernel.edge_rows, prefix)
        value = (np.einsum("bja,bja->bj", displacement, displacement)
                 @ kernel.outer_weights) / kernel.duration
        keep = np.argsort(value)[:_REFINE]
        leaders.extend((float(value[i]), block[i]) for i in keep)
        leaders.sort(key=lambda pair: pair[0])
        del leaders[_REFINE:]

    best_value, best_columns = np.inf, None
    for _, candidate in leaders:
        schedule = degrees[candidate].astype(np.int64)
        exact = kernel.objective(problem.gather(schedule))
        if exact < best_value:
            best_value, best_columns = exact, candidate

    schedule = degrees[best_columns].astype(np.int64)
    return ScheduleSolution(
        degrees=schedule,
        objective=best_value,
        work=schedule_work(schedule, problem.time_weight),
        budget=budget,
        multiplier=float("nan"),
        feasible=True,
        diagnostics={"n_total": float(total),
                     "n_feasible": float(columns.shape[0])},
    )


@dataclass(frozen=True, slots=True)
class VerificationRecord:
    """What one enumerable instance says about a solver.

    Attributes
    ----------
    optimum:
        :math:`J` of the true integer optimum.
    attained:
        :math:`J` of the schedule under test.
    ratio:
        ``attained / optimum``.  One means exactly optimal; the campaign
        reports the distribution of this over the panel, not its mean.
    optimal:
        Whether the ratio is one to within tolerance.
    n_intervals, n_candidates, n_feasible:
        The size of the instance and how much of it the ceiling admitted.
    """

    optimum: float
    attained: float
    ratio: float
    optimal: bool
    n_intervals: int
    n_candidates: int
    n_feasible: int


def verify_schedule(problem: DescentProblem, budget: float,
                    schedule: IntArr, tolerance: float = 1.0e-9,
                    max_schedules: int = MAX_SCHEDULES) -> VerificationRecord:
    """Compare a solver's schedule against the enumerated optimum.

    Parameters
    ----------
    problem, budget:
        The instance.
    schedule:
        What the solver returned.
    tolerance:
        Relative slack allowed before ``optimal`` is false.
    max_schedules:
        Passed through to :func:`solve_exhaustive`.

    Returns
    -------
    VerificationRecord

    Raises
    ------
    ValueError
        If the schedule overspends the ceiling.  A schedule that beats the
        optimum by breaking the constraint is not evidence about the solver,
        it is evidence about the constraint, and conflating the two is how a
        verification panel comes to endorse an infeasible answer.
    """
    schedule = np.asarray(schedule, dtype=np.int64)
    work = schedule_work(schedule, problem.time_weight)
    if work > budget * (1.0 + 1.0e-9):
        raise ValueError(
            f"the schedule under test spends {work:.6g} against a ceiling of "
            f"{budget:.6g}; it is not a feasible point and cannot be compared "
            "with the constrained optimum")

    reference = solve_exhaustive(problem, budget, max_schedules)
    attained = problem.kernel.objective(problem.gather(schedule))
    ratio = attained / reference.objective if reference.objective > 0.0 else 1.0
    return VerificationRecord(
        optimum=reference.objective,
        attained=attained,
        ratio=ratio,
        optimal=ratio <= 1.0 + tolerance,
        n_intervals=problem.n_intervals,
        n_candidates=len(problem.candidate_degrees),
        n_feasible=int(reference.diagnostics["n_feasible"]),
    )


def panel_summary(records: list[VerificationRecord]) -> dict[str, float]:
    """Reduce a panel of records to the numbers a results table carries.

    Returns
    -------
    dict
        ``n``, ``n_optimal``, ``fraction_optimal``, ``worst_ratio`` and
        ``median_ratio``.  The worst case is reported beside the fraction
        because a panel that is optimal nine times in ten and 40 per cent off
        the tenth is a different solver from one that is never more than a
        per cent out.

    Raises
    ------
    ValueError
        If the panel is empty.
    """
    if not records:
        raise ValueError("an empty panel summarises to nothing, not to zero")
    ratios = np.array([r.ratio for r in records], dtype=float)
    optimal = sum(1 for r in records if r.optimal)
    return {
        "n": float(len(records)),
        "n_optimal": float(optimal),
        "fraction_optimal": optimal / len(records),
        "worst_ratio": float(ratios.max()),
        "median_ratio": float(np.median(ratios)),
    }
