"""Verification of the verifier.

The panel is now what the campaign's claim about solver quality rests on, so
the enumeration itself has to be beyond doubt. Three properties carry it:

* it finds the true optimum, checked against an independent brute force
  written the obvious slow way -- if the batched evaluation and the plain
  loop disagree, the fast one is wrong;
* it never returns an infeasible schedule, and never accepts one for
  comparison;
* it refuses, loudly, when the enumeration is too large, rather than becoming
  a job nobody runs.
"""

from __future__ import annotations

import itertools

import numpy as np
import pytest

from tda.allocate import (
    DescentProblem,
    panel_summary,
    schedule_count,
    schedule_work,
    solve_descent,
    solve_exhaustive,
    verify_schedule,
)
from tda.allocate.exhaustive import VerificationRecord
from tda.kernel import CouplingKernel

DEGREES = (10, 20, 30)
CELLS_PER_INTERVAL = 2
INTERVALS = 5
M_CELLS = INTERVALS * CELLS_PER_INTERVAL


class _Table:
    def __init__(self, defect, node_transport, widths, degrees):
        self.defect = defect
        self.node_transport = node_transport
        self.widths = widths
        self.schema = type("S", (), {"candidate_degrees": degrees})()


@pytest.fixture
def problem():
    """A tiny instance, matching the one the solver tests use."""
    rng = np.random.default_rng(4242)
    edge_rows = rng.normal(size=(M_CELLS + 1, 3, 6))
    weights = rng.uniform(0.5, 1.5, M_CELLS + 1)
    kernel = CouplingKernel.from_arc(edge_rows, weights, float(weights.sum()))
    defect = rng.normal(size=(M_CELLS, len(DEGREES), 3))
    defect *= np.array([1.0, 0.4, 0.15])[None, :, None]
    node_transport = rng.normal(size=(M_CELLS, 6, 3))
    widths = rng.uniform(0.5, 1.5, M_CELLS)
    interval_of = np.repeat(np.arange(INTERVALS), CELLS_PER_INTERVAL)
    time_weight = np.add.reduceat(widths,
                                  np.arange(0, M_CELLS, CELLS_PER_INTERVAL))
    return DescentProblem.from_table(
        _Table(defect, node_transport, widths, DEGREES), kernel, interval_of,
        time_weight)


def _slow_optimum(problem, budget):
    """The obvious loop, written to disagree with the fast path if it can."""
    best, argbest = np.inf, None
    for combo in itertools.product(DEGREES, repeat=problem.n_intervals):
        schedule = np.array(combo, dtype=np.int64)
        if schedule_work(schedule, problem.time_weight) > budget:
            continue
        value = problem.kernel.objective(problem.gather(schedule))
        if value < best:
            best, argbest = value, schedule
    return best, argbest


def _mid_budget(problem):
    lo = schedule_work(np.full(INTERVALS, DEGREES[0]), problem.time_weight)
    hi = schedule_work(np.full(INTERVALS, DEGREES[-1]), problem.time_weight)
    return 0.5 * (lo + hi)


def test_the_batched_enumeration_agrees_with_the_obvious_loop(problem) -> None:
    """One is fast and one is readable; they must return the same schedule."""
    for fraction in (0.35, 0.5, 0.8):
        budget = fraction * schedule_work(np.full(INTERVALS, DEGREES[-1]),
                                          problem.time_weight)
        slow, slow_schedule = _slow_optimum(problem, budget)
        fast = solve_exhaustive(problem, budget)
        assert fast.objective == pytest.approx(slow, rel=1e-12)
        assert fast.degrees.tolist() == slow_schedule.tolist()


def test_the_optimum_obeys_the_ceiling(problem) -> None:
    budget = _mid_budget(problem)
    result = solve_exhaustive(problem, budget)
    assert result.work <= budget * (1.0 + 1e-12)
    assert result.feasible
    assert result.diagnostics["n_total"] == schedule_count(problem)
    assert 0 < result.diagnostics["n_feasible"] <= schedule_count(problem)


def test_a_chunk_boundary_does_not_change_the_answer(problem) -> None:
    """The leader board is kept across batches, not per batch."""
    budget = _mid_budget(problem)
    answers = [solve_exhaustive(problem, budget, chunk=c).degrees.tolist()
               for c in (1, 7, 64, 10_000)]
    assert all(a == answers[0] for a in answers)


def test_an_unaffordable_enumeration_is_refused_by_name(problem) -> None:
    with pytest.raises(ValueError, match="above the cap"):
        solve_exhaustive(problem, _mid_budget(problem), max_schedules=10)


def test_a_ceiling_below_the_cheapest_schedule_is_an_error(problem) -> None:
    with pytest.raises(ValueError, match="no schedule fits"):
        solve_exhaustive(problem, 1.0e-9)


def test_the_declared_descent_is_checked_against_truth(problem) -> None:
    """The panel's actual job."""
    budget = _mid_budget(problem)
    descent = solve_descent(problem, budget,
                            [np.full(INTERVALS, DEGREES[0], dtype=np.int64)])
    record = verify_schedule(problem, budget, descent.degrees)
    assert record.ratio >= 1.0 - 1e-12
    assert record.optimal == (record.ratio <= 1.0 + 1e-9)
    assert record.n_intervals == INTERVALS
    assert record.n_candidates == len(DEGREES)


def test_an_infeasible_schedule_is_not_compared(problem) -> None:
    """Beating the optimum by overspending is evidence about the constraint."""
    budget = _mid_budget(problem)
    with pytest.raises(ValueError, match="not a feasible point"):
        verify_schedule(problem, budget, np.full(INTERVALS, DEGREES[-1]))


def test_a_suboptimal_schedule_is_reported_as_suboptimal(problem) -> None:
    budget = _mid_budget(problem)
    optimum = solve_exhaustive(problem, budget)
    worse = np.full(INTERVALS, DEGREES[0], dtype=np.int64)
    if worse.tolist() == optimum.degrees.tolist():
        pytest.skip("the cheapest schedule is the optimum on this instance")
    record = verify_schedule(problem, budget, worse)
    assert not record.optimal
    assert record.ratio > 1.0


def test_the_panel_reports_the_worst_case_beside_the_fraction() -> None:
    """Nine optimal and one badly off is a different solver from ten close."""
    def record(ratio: float) -> VerificationRecord:
        return VerificationRecord(optimum=1.0, attained=ratio, ratio=ratio,
                                  optimal=ratio <= 1.0 + 1e-9, n_intervals=4,
                                  n_candidates=3, n_feasible=10)

    summary = panel_summary([record(1.0)] * 9 + [record(1.4)])
    assert summary["n"] == 10.0
    assert summary["n_optimal"] == 9.0
    assert summary["fraction_optimal"] == pytest.approx(0.9)
    assert summary["worst_ratio"] == pytest.approx(1.4)
    assert summary["median_ratio"] == pytest.approx(1.0)

    with pytest.raises(ValueError, match="empty panel"):
        panel_summary([])
