"""Verification of the allocation solvers.

The properties that carry the campaign, in order of how badly a silent failure
would hurt:

* the certificate really is a **lower** bound -- checked against an exhaustive
  optimum on instances small enough to enumerate;
* the descent's objective really is the kernel's, evaluated on the schedule it
  returns, not a bookkeeping quantity that drifted from it;
* the budget is a ceiling and is never exceeded;
* selection among starts is on the penalized objective, which is the one
  place a plausible implementation is wrong.

Everything is built on random contributions rather than a propagated arc: an
exhaustive optimum is only affordable on a tiny instance, and random data is
stronger than physical data for an algebraic claim.
"""

from __future__ import annotations

import itertools

import numpy as np
import pytest

from tda.allocate import (
    DescentProblem,
    argmax_rounding,
    certify,
    force_values,
    linear_minimisation,
    monotonicity_report,
    round_and_polish,
    schedule_work,
    sensitivity_values,
    solve_descent,
    solve_separable,
    solve_to_budget,
    sum_up_rounding,
)
from tda.kernel import CouplingKernel

DEGREES = (10, 20, 30)
CELLS_PER_INTERVAL = 2
INTERVALS = 5
M_CELLS = INTERVALS * CELLS_PER_INTERVAL


class _Table:
    """The three arrays :class:`DescentProblem.from_table` reads."""

    def __init__(self, defect, node_transport, widths, degrees):
        self.defect = defect
        self.node_transport = node_transport
        self.widths = widths
        self.schema = type("S", (), {"candidate_degrees": degrees})()


@pytest.fixture
def problem():
    """A tiny instance whose optimum can be found by enumeration."""
    rng = np.random.default_rng(4242)
    edge_rows = rng.normal(size=(M_CELLS + 1, 3, 6))
    weights = rng.uniform(0.5, 1.5, M_CELLS + 1)
    kernel = CouplingKernel.from_arc(edge_rows, weights, float(weights.sum()))

    defect = rng.normal(size=(M_CELLS, len(DEGREES), 3))
    # Higher degree, smaller defect -- the physical direction, so that a
    # budget actually bites.
    defect *= np.array([1.0, 0.4, 0.15])[None, :, None]
    node_transport = rng.normal(size=(M_CELLS, 6, 3))
    widths = rng.uniform(0.5, 1.5, M_CELLS)
    table = _Table(defect, node_transport, widths, DEGREES)

    interval_of = np.repeat(np.arange(INTERVALS), CELLS_PER_INTERVAL)
    time_weight = np.add.reduceat(widths,
                                  np.arange(0, M_CELLS, CELLS_PER_INTERVAL))
    return DescentProblem.from_table(table, kernel, interval_of, time_weight)


def _exhaustive(problem, budget):
    """Brute-force optimum over every integer schedule that fits."""
    best, best_schedule = np.inf, None
    for combo in itertools.product(DEGREES, repeat=problem.n_intervals):
        schedule = np.array(combo, dtype=np.int64)
        if schedule_work(schedule, problem.time_weight) > budget:
            continue
        value = problem.kernel.objective(problem.gather(schedule))
        if value < best:
            best, best_schedule = value, schedule
    return best, best_schedule


def _mid_budget(problem):
    """A ceiling that binds: between the cheapest and the most expensive."""
    lo = schedule_work(np.full(INTERVALS, DEGREES[0]), problem.time_weight)
    hi = schedule_work(np.full(INTERVALS, DEGREES[-1]), problem.time_weight)
    return 0.5 * (lo + hi)


# ---------------------------------------------------------------------------
# The certificate is a lower bound
# ---------------------------------------------------------------------------


def test_the_bound_never_exceeds_the_true_optimum(problem) -> None:
    """The one property the certificate exists for.

    A bound above the optimum would make every reported gap an understatement
    and nothing else in the campaign would notice.
    """
    budget = _mid_budget(problem)
    optimum, _ = _exhaustive(problem, budget)
    result = certify(problem, optimum, budget, iterations=40)
    assert result.lower_bound <= optimum * (1.0 + 1e-9)


def test_the_bound_is_not_trivially_zero(problem) -> None:
    """A vacuous bound is honest but useless.

    On this instance it should actually say something.
    """
    budget = _mid_budget(problem)
    optimum, _ = _exhaustive(problem, budget)
    result = certify(problem, optimum, budget, iterations=60)
    assert not result.vacuous
    assert result.lower_bound > 0.0


def test_the_error_gap_is_the_square_root_of_the_objective_gap(problem) -> None:
    """``g_E = 1 - sqrt(1 - g_J)``: ten per cent in error is nineteen in J."""
    budget = _mid_budget(problem)
    optimum, _ = _exhaustive(problem, budget)
    result = certify(problem, optimum, budget, iterations=30)
    assert result.gap_error == pytest.approx(
        1.0 - np.sqrt(1.0 - result.gap_objective), rel=1e-12)


def test_the_naming_rule_is_bound_to_the_gap_not_the_result(problem) -> None:
    budget = _mid_budget(problem)
    optimum, _ = _exhaustive(problem, budget)
    result = certify(problem, optimum, budget, iterations=60)
    assert result.earns_the_name(1.0) is not result.earns_the_name(0.0)


def test_the_subproblem_solution_is_a_feasible_vertex(problem) -> None:
    """Simplex rows and the knapsack, both respected by the LP's answer."""
    budget = _mid_budget(problem)
    schedule = np.full(INTERVALS, DEGREES[1])
    gradient = problem.kernel.gradient(problem.gather(schedule))
    theta = linear_minimisation(problem, gradient, budget)

    assert theta.shape == (INTERVALS, len(DEGREES))
    assert np.all(theta > -1e-12)
    assert np.allclose(theta.sum(axis=1), 1.0, atol=1e-9)
    work = float((problem.time_weight[:, None]
                  * np.array(DEGREES, dtype=float) ** 2 * theta).sum())
    assert work <= budget * (1.0 + 1e-9)


def test_the_subproblem_is_solved_over_the_set_not_its_boundary(
        problem) -> None:
    """A ceiling far above anything reachable must not pull the answer onto it.

    Forcing the budget to bind is the failure that would silently invalidate
    the bound (D142), so the case where it must not bind is pinned.
    """
    generous = 1.0e6 * _mid_budget(problem)
    schedule = np.full(INTERVALS, DEGREES[0])
    gradient = problem.kernel.gradient(problem.gather(schedule))
    theta = linear_minimisation(problem, gradient, generous)
    work = float((problem.time_weight[:, None]
                  * np.array(DEGREES, dtype=float) ** 2 * theta).sum())
    assert work < generous * 0.5


# ---------------------------------------------------------------------------
# The descent
# ---------------------------------------------------------------------------


def test_descent_reports_the_kernel_objective_of_its_own_schedule(
        problem) -> None:
    """The bookkeeping and the objective must not have drifted apart."""
    budget = _mid_budget(problem)
    starts = [np.full(INTERVALS, d) for d in DEGREES]
    result = solve_descent(problem, budget, starts)
    recomputed = problem.kernel.objective(problem.gather(result.degrees))
    assert result.objective == pytest.approx(recomputed, rel=1e-12)


def test_descent_respects_the_ceiling(problem) -> None:
    budget = _mid_budget(problem)
    starts = [np.full(INTERVALS, d) for d in DEGREES]
    result = solve_descent(problem, budget, starts)
    assert result.work <= budget
    assert result.feasible


def test_descent_is_close_to_the_exhaustive_optimum(problem) -> None:
    """A local method may fall short, but not by much on an instance this small.

    Exactness is not required -- this is a local method on a non-convex integer
    problem -- but a large gap would mean the sweep is broken rather than
    merely local.
    """
    budget = _mid_budget(problem)
    optimum, _ = _exhaustive(problem, budget)
    starts = [np.full(INTERVALS, d) for d in DEGREES]
    result = solve_descent(problem, budget, starts)
    assert result.objective >= optimum * (1.0 - 1e-12)
    assert result.objective <= optimum * 1.20


def test_an_unbinding_ceiling_leaves_the_multiplier_at_zero(problem) -> None:
    """KKT complementary slackness, and the case that motivated it.

    With a ceiling nothing can reach, the answer is the unconstrained optimum
    and the multiplier must be exactly zero -- not a small positive number
    that a bisection happened to stop at.
    """
    generous = 1.0e6 * _mid_budget(problem)
    starts = [np.full(INTERVALS, DEGREES[0])]
    result = solve_descent(problem, generous, starts)
    assert result.multiplier == 0.0
    assert not result.ceiling_binds
    assert result.utilisation < 1.0


def test_selection_among_starts_uses_the_penalized_objective(problem) -> None:
    """Selection must be on the penalized objective, which is the D131 failure.

    At a fixed multiplier, choosing on ``J`` alone prefers whichever basin
    bought a lower objective by spending more.  Checked at the level that
    matters: the returned schedule must minimise ``J + lambda W`` among the
    starts, not ``J``.
    """
    budget = _mid_budget(problem)
    starts = [np.full(INTERVALS, d) for d in DEGREES]
    result = solve_descent(problem, budget, starts)
    chosen = result.objective + result.multiplier * result.work
    for start in starts:
        value = problem.kernel.objective(problem.gather(start))
        work = schedule_work(start, problem.time_weight)
        if work <= budget:
            assert chosen <= value + result.multiplier * work + 1e-9


def test_the_spread_across_starts_is_reported(problem) -> None:
    budget = _mid_budget(problem)
    starts = [np.full(INTERVALS, d) for d in DEGREES]
    result = solve_descent(problem, budget, starts)
    assert len(result.spread) == len(starts)
    assert min(result.spread) == pytest.approx(1.0)


def test_an_untabulated_start_is_refused(problem) -> None:
    with pytest.raises(ValueError, match="not a tabulated candidate"):
        solve_descent(problem, _mid_budget(problem),
                      [np.full(INTERVALS, 17)])


def test_no_start_at_all_is_refused(problem) -> None:
    with pytest.raises(ValueError, match="at least one starting schedule"):
        solve_descent(problem, _mid_budget(problem), [])


def test_monotonicity_is_measured_rather_than_assumed(problem) -> None:
    """Monotonicity holds for a global minimiser, and this is not one.

    So the campaign measures it rather than assuming it (D33).
    """
    starts = [np.full(INTERVALS, d) for d in DEGREES]
    report = monotonicity_report(problem, starts, np.logspace(-14, -6, 12))
    assert set(report) == {"violations", "worst_rise", "monotone"}
    assert report["worst_rise"] >= 0.0


# ---------------------------------------------------------------------------
# The separable rungs
# ---------------------------------------------------------------------------


def test_separable_solve_is_exact_at_a_fixed_multiplier(problem) -> None:
    """Each interval decouples, so the solve is a global minimum.

    It is compared against brute force on the penalized criterion, which is only
    affordable because the problem separates -- which is the point.
    """
    rng = np.random.default_rng(8)
    values = rng.uniform(0.0, 1.0, (INTERVALS, len(DEGREES)))
    budget = _mid_budget(problem)
    result = solve_separable(values, problem.time_weight, DEGREES, budget)

    penalized = (values + result.multiplier * problem.time_weight[:, None]
                 * np.array(DEGREES, dtype=float) ** 2)
    expected = np.array(DEGREES)[np.argmin(penalized, axis=1)]
    assert np.array_equal(result.degrees, expected)


def test_separable_respects_the_ceiling(problem) -> None:
    rng = np.random.default_rng(9)
    values = rng.uniform(0.0, 1.0, (INTERVALS, len(DEGREES)))
    budget = _mid_budget(problem)
    result = solve_separable(values, problem.time_weight, DEGREES, budget)
    assert result.work <= budget


def test_force_and_sensitivity_values_have_the_right_shape(problem) -> None:
    rng = np.random.default_rng(10)
    defect = rng.normal(size=(M_CELLS, len(DEGREES), 3))
    widths = np.ones(M_CELLS)
    local = np.stack([np.eye(3)] * M_CELLS)
    force = force_values(defect, widths, problem.interval_of, INTERVALS)
    sens = sensitivity_values(defect, widths, local, problem.interval_of,
                              INTERVALS)
    assert force.shape == sens.shape == (INTERVALS, len(DEGREES))
    # With K = I the two criteria coincide, which pins the contraction order.
    assert np.allclose(force, sens, rtol=1e-12, atol=0.0)


def test_sensitivity_rejects_a_misshaped_kernel(problem) -> None:
    defect = np.zeros((M_CELLS, len(DEGREES), 3))
    with pytest.raises(ValueError, match=r"shape \(10, 3, 3\)"):
        sensitivity_values(defect, np.ones(M_CELLS), np.zeros((3, 3, 3)),
                           problem.interval_of, INTERVALS)


# ---------------------------------------------------------------------------
# Rounding
# ---------------------------------------------------------------------------


def test_sum_up_rounding_tracks_the_running_total() -> None:
    """A relaxation that splits evenly must alternate, not pick one mode."""
    theta = np.tile(np.array([0.5, 0.5, 0.0]), (4, 1))
    got = sum_up_rounding(theta, (10, 20, 30))
    assert sorted(got.tolist()) == [10, 10, 20, 20]


def test_argmax_rounding_takes_the_mode() -> None:
    theta = np.array([[0.6, 0.4, 0.0], [0.1, 0.2, 0.7]])
    assert list(argmax_rounding(theta, (10, 20, 30))) == [10, 30]


def test_sum_up_refuses_a_row_that_is_not_a_distribution() -> None:
    """A row that does not sum to one means the relaxation did not solve."""
    with pytest.raises(ValueError, match="sum to one"):
        sum_up_rounding(np.array([[0.5, 0.2, 0.0]]), (10, 20, 30))


def test_round_and_polish_returns_a_feasible_schedule(problem) -> None:
    budget = _mid_budget(problem)
    starts = [np.full(INTERVALS, d) for d in DEGREES]
    descent = solve_descent(problem, budget, starts)
    gradient = problem.kernel.gradient(problem.gather(descent.degrees))
    theta = linear_minimisation(problem, gradient, budget)

    result = round_and_polish(problem, theta, budget, descent.multiplier)
    assert result.work <= budget
    assert result.objective == pytest.approx(
        problem.kernel.objective(problem.gather(result.degrees)), rel=1e-12)
    assert "sum_up_won" in result.diagnostics


# ---------------------------------------------------------------------------
# The shared budget contract
# ---------------------------------------------------------------------------


def test_solve_to_budget_reports_an_impossible_ceiling(problem) -> None:
    """Below the cheapest schedule there is nothing, and it must say so."""
    def solve_at(_):
        return np.full(INTERVALS, DEGREES[0]), 1.0

    with pytest.raises(ValueError, match="no multiplier up to"):
        solve_to_budget(solve_at, problem.time_weight, budget=1e-12)


def test_solve_to_budget_rejects_a_non_positive_ceiling(problem) -> None:
    with pytest.raises(ValueError, match="budget must be positive"):
        solve_to_budget(lambda _: (np.zeros(INTERVALS, int), 0.0),
                        problem.time_weight, budget=0.0)
