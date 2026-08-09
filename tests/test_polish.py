"""The hard-budget polish, and the integer-program failure it exists for.

A multiplier sweep minimises ``J + lambda W``. Over a discrete candidate set
that is a weighted-sum scalarisation, and it reaches only the *supported*
points of the trade-off curve; a constrained optimum in a non-convex notch is
unreachable from every multiplier. Complementary slackness fails for the same
reason -- the optimum can leave budget unspent at every positive lambda.

These tests pin the consequences: the continuation is not claimed to be
optimal, the polish searches under the constraint as written, and where the
instance is small enough the two together are checked against the enumerated
truth rather than against a bound.
"""

from __future__ import annotations

import numpy as np
import pytest

from tda.allocate import (
    DescentProblem,
    block_deltas,
    polish_to_budget,
    schedule_work,
    solve_descent,
    solve_exhaustive,
    solve_to_budget,
    verify_schedule,
)
from tda.kernel import CouplingKernel

DEGREES = (10, 20, 30)
CELLS_PER_INTERVAL = 2
INTERVALS = 6
M_CELLS = INTERVALS * CELLS_PER_INTERVAL


class _Table:
    def __init__(self, defect, node_transport, widths, degrees):
        self.defect = defect
        self.node_transport = node_transport
        self.widths = widths
        self.schema = type("S", (), {"candidate_degrees": degrees})()


def _problem(seed: int) -> DescentProblem:
    rng = np.random.default_rng(seed)
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


@pytest.fixture
def problem():
    return _problem(90210)


def _mid_budget(problem):
    lo = schedule_work(np.full(INTERVALS, DEGREES[0]), problem.time_weight)
    hi = schedule_work(np.full(INTERVALS, DEGREES[-1]), problem.time_weight)
    return 0.5 * (lo + hi)


# ---------------------------------------------------------------------------
# The premise: complementary slackness is not an optimality condition here
# ---------------------------------------------------------------------------


def test_the_optimum_can_leave_budget_unspent_at_every_multiplier() -> None:
    """The advisor's two-candidate counterexample, as a property of the solver.

    One interval, candidates costing 5 and 10 with objectives 1 and 0, ceiling
    6. The constrained optimum is the cheap one: it spends 5 of 6 and no
    multiplier makes ``lambda (B - W)`` vanish. A solver that stopped on
    complementary slackness would reject the right answer.
    """
    # One interval, unit time weight, so work is the squared "degree": the two
    # candidates cost 5 and 10 exactly.
    candidates = np.array([np.sqrt(5.0), np.sqrt(10.0)])
    values = np.array([1.0, 0.0])

    def solve_at(multiplier: float):
        picked = int(np.argmin(values + multiplier * candidates**2))
        return candidates[picked:picked + 1], float(values[picked])

    solution = solve_to_budget(solve_at, np.array([1.0]), 6.0)

    assert solution.objective == pytest.approx(1.0)
    assert solution.work == pytest.approx(5.0)
    slack = 6.0 - solution.work
    assert slack > 0.0
    assert solution.multiplier > 0.0
    assert solution.multiplier * slack > 0.0     # slackness does NOT vanish


def test_the_continuation_keeps_the_best_feasible_not_the_last(problem) -> None:
    """More starting points must never make the answer worse.

    Returning whichever schedule the final bracket landed on lets a larger
    start set come back above a smaller one, which was measured on a pilot arc
    before it was fixed.
    """
    budget = _mid_budget(problem)
    one = solve_descent(problem, budget,
                        [np.full(INTERVALS, DEGREES[0], dtype=np.int64)])
    many = solve_descent(problem, budget,
                         [np.full(INTERVALS, d, dtype=np.int64)
                          for d in DEGREES])
    assert many.objective <= one.objective * (1.0 + 1e-12)
    assert many.work <= budget * (1.0 + 1e-9)


# ---------------------------------------------------------------------------
# The polish itself
# ---------------------------------------------------------------------------


def test_the_single_move_deltas_are_exact(problem) -> None:
    """Each entry must be the true change in J from moving one interval."""
    budget = _mid_budget(problem)
    columns = np.zeros(INTERVALS, dtype=np.int64)
    degrees = np.asarray(problem.candidate_degrees)
    delta_j, delta_w = block_deltas(problem, columns)
    base = problem.kernel.objective(problem.gather(degrees[columns]))

    for q in range(INTERVALS):
        for p in range(len(DEGREES)):
            trial = columns.copy()
            trial[q] = p
            moved = problem.kernel.objective(problem.gather(degrees[trial]))
            assert delta_j[q, p] == pytest.approx(moved - base,
                                                  rel=1e-9, abs=1e-12)
            assert delta_w[q, p] == pytest.approx(
                schedule_work(degrees[trial], problem.time_weight)
                - schedule_work(degrees[columns], problem.time_weight),
                rel=1e-9, abs=1e-12)
    assert budget > 0.0


def test_the_polish_never_overspends_and_never_worsens(problem) -> None:
    budget = _mid_budget(problem)
    for seed_degree in DEGREES:
        start = np.full(INTERVALS, seed_degree, dtype=np.int64)
        if schedule_work(start, problem.time_weight) > budget:
            continue
        polished, report = polish_to_budget(problem, start, budget)
        assert report.work_after <= budget * (1.0 + 1e-9)
        assert report.objective_after <= report.objective_before + 1e-12
        assert schedule_work(polished, problem.time_weight) == pytest.approx(
            report.work_after)


def test_the_exchange_step_reaches_what_single_moves_cannot(problem) -> None:
    """A raise blocked by the budget becomes affordable once something drops.

    If the exchange never mattered the two runs would agree; they are asserted
    to differ on at least one of a small set of instances, so the step is
    justified by measurement rather than by the argument for it.
    """
    improved_somewhere = False
    for seed in (90210, 7, 31337, 555):
        instance = _problem(seed)
        budget = _mid_budget(instance)
        start = np.full(INTERVALS, DEGREES[0], dtype=np.int64)
        _, with_exchange = polish_to_budget(instance, start, budget)
        _, without = polish_to_budget(instance, start, budget, exchange=False)
        assert with_exchange.objective_after <= without.objective_after + 1e-12
        if with_exchange.objective_after < without.objective_after - 1e-12:
            improved_somewhere = True
            assert with_exchange.exchange_moves > 0
    assert improved_somewhere


def test_an_infeasible_start_is_refused(problem) -> None:
    with pytest.raises(ValueError, match="the polish improves a feasible"):
        polish_to_budget(problem, np.full(INTERVALS, DEGREES[-1]), 1.0)


def test_an_untabulated_degree_is_refused(problem) -> None:
    with pytest.raises(ValueError, match="not a tabulated candidate"):
        polish_to_budget(problem, np.full(INTERVALS, 17), _mid_budget(problem))


# ---------------------------------------------------------------------------
# Against truth
# ---------------------------------------------------------------------------


def _pipeline_ratios() -> list[float]:
    """``J_pipeline / J*`` on four enumerable random instances."""
    ratios = []
    for seed in (90210, 7, 31337, 555):
        instance = _problem(seed)
        budget = _mid_budget(instance)
        solution = solve_descent(
            instance, budget,
            [np.full(INTERVALS, d, dtype=np.int64) for d in DEGREES])
        record = verify_schedule(instance, budget, solution.degrees)
        ratios.append(record.ratio)
    return ratios


def test_the_pipeline_never_beats_the_enumerated_optimum() -> None:
    """Feasible, and never below truth. The correctness half."""
    for ratio in _pipeline_ratios():
        assert ratio >= 1.0 - 1e-9


def test_the_pipeline_is_not_optimal_in_general_and_the_tests_say_so() -> None:
    """Pinned so that nobody later writes that the solver is optimal.

    On the *structured* instances the pilot enumerated -- real field, real arc
    -- the descent returned the exact optimum five times out of five. On
    unstructured random contributions it does not: at least one of these four
    lands strictly above, and on the worst of them by tens of per cent. Both
    facts are real and the campaign reports the distribution rather than a
    bound, which is what the panel of :mod:`tda.allocate.exhaustive` is for.
    """
    ratios = _pipeline_ratios()
    assert max(ratios) > 1.0 + 1e-6
    assert min(ratios) == pytest.approx(1.0, abs=1e-9)


def test_the_polish_moves_toward_the_optimum_when_the_sweep_stops_short(
        problem) -> None:
    budget = _mid_budget(problem)
    starts = [np.full(INTERVALS, d, dtype=np.int64) for d in DEGREES]
    raw = solve_descent(problem, budget, starts, polish=False)
    done = solve_descent(problem, budget, starts)
    optimum = solve_exhaustive(problem, budget)

    assert done.objective <= raw.objective + 1e-12
    assert done.objective >= optimum.objective - 1e-9
    assert "polish_single_moves" in done.diagnostics
    assert "polish_single_moves" not in raw.diagnostics
