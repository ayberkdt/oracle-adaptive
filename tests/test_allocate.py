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
from dataclasses import replace

import numpy as np
import pytest

from tda.allocate import (
    DescentProblem,
    argmax_rounding,
    certify,
    force_values,
    greedy_fill,
    linear_minimisation,
    monotonicity_report,
    objective_coefficients,
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


def test_away_steps_do_not_change_what_the_bound_is(problem) -> None:
    """Validity is the away step's one non-negotiable.

    It changes how fast the iterate reaches the relaxed optimum and nothing
    about what the linearisation at that iterate means, so the bound must
    still sit under the true optimum.
    """
    budget = _mid_budget(problem)
    optimum, _ = _exhaustive(problem, budget)
    result = certify(problem, optimum, budget, iterations=60, away_steps=True)
    assert result.lower_bound <= optimum * (1.0 + 1e-9)
    assert result.away_steps > 0
    assert result.active_atoms >= 1


def test_away_steps_close_the_relaxation_and_classical_steps_do_not(
        problem) -> None:
    """What the away step actually buys: termination.

    On this instance the away-step iteration drives the relaxation's own
    duality gap to zero, so its bound *is* the relaxed optimum and no further
    computation can improve it. Classical Frank--Wolfe is still short of that
    after more than twice the steps.
    """
    budget = _mid_budget(problem)
    optimum, _ = _exhaustive(problem, budget)
    away = certify(problem, optimum, budget, iterations=800, away_steps=True)
    plain = certify(problem, optimum, budget, iterations=2000,
                    away_steps=False)

    # Converged to the requested tolerance, which is relative to the scale
    # being certified rather than absolute.
    assert away.relaxed_gap <= 1e-9 * optimum
    assert away.lower_bound == pytest.approx(away.relaxed_objective, rel=1e-8)
    assert plain.relaxed_gap > 1e-6
    assert away.lower_bound >= plain.lower_bound


def test_starting_from_the_schedule_being_certified_is_worth_more_than_steps(
        problem) -> None:
    """Where the iteration starts dominates how long it runs.

    From the cheapest vertex the iterate has the whole objective to descend
    before the bound means anything, and the forward direction wins every
    comparison on the way -- so the away step never fires and the run is
    classical Frank--Wolfe with extra bookkeeping. Started at the schedule
    under test it begins where the answer is.
    """
    budget = _mid_budget(problem)
    optimum, schedule = _exhaustive(problem, budget)
    cold = certify(problem, optimum, budget, iterations=20)
    warm = certify(problem, optimum, budget, iterations=20,
                   start_schedule=schedule)
    assert warm.relaxed_gap < cold.relaxed_gap
    assert warm.lower_bound > cold.lower_bound
    assert warm.lower_bound <= optimum * (1.0 + 1e-9)


def test_a_start_outside_the_ceiling_is_refused(problem) -> None:
    """A bound is only a bound if the iterate began inside the feasible set."""
    budget = _mid_budget(problem)
    with pytest.raises(ValueError, match="against a ceiling"):
        certify(problem, 1.0, budget, iterations=5,
                start_schedule=np.full(INTERVALS, DEGREES[-1]))
    with pytest.raises(ValueError, match="not a tabulated candidate"):
        certify(problem, 1.0, budget, iterations=5,
                start_schedule=np.full(INTERVALS, 17))


def test_an_emptied_vertex_does_not_end_the_iteration(problem) -> None:
    """The away step must not be offered a vertex carrying no weight.

    A full Frank--Wolfe step leaves every previous atom at zero weight. If the
    away search then picks one, its ceiling is zero, the step is zero, and the
    loop reads that as convergence -- producing a certificate that looks
    settled after three steps and is worth nothing. It stops early on real
    data and not on the random fixture, so the property is asserted directly
    rather than left to whichever instance happens to trigger it.
    """
    budget = _mid_budget(problem)
    optimum, _ = _exhaustive(problem, budget)
    for iterations in (50, 200, 800):
        result = certify(problem, optimum, budget, iterations=iterations)
        # Either it used the steps it was given, or it stopped because the
        # relaxation was solved -- never because it ran out of directions.
        settled = result.relaxed_gap <= 1e-9 * optimum
        assert result.iterations == iterations or settled
        assert result.active_atoms >= 1


def test_a_solved_relaxation_stops_instead_of_spending_its_budget(
        problem) -> None:
    """Once the gap is zero no further step can raise the bound.

    Worth stopping on rather than iterating through: the campaign certifies
    one orbit after another, and iterations spent after convergence are
    iterations the next orbit does not get.
    """
    budget = _mid_budget(problem)
    optimum, _ = _exhaustive(problem, budget)
    generous = certify(problem, optimum, budget, iterations=5000)
    exact = certify(problem, optimum, budget, iterations=800)
    assert generous.iterations < 5000
    assert generous.lower_bound == pytest.approx(exact.lower_bound, rel=1e-9)


def test_the_bound_is_a_maximum_over_a_sequence_that_is_not_monotone(
        problem) -> None:
    """Why the away step is not uniformly better at a fixed step count.

    The reported bound is the best linearisation seen along the path, and the
    two variants take different paths, so at low step counts the classical one
    can happen to pass through a better point. It stops mattering once the
    relaxation is solved, but a run that compared the two at sixty steps and
    concluded the away step had hurt would be reading path noise.
    """
    budget = _mid_budget(problem)
    optimum, _ = _exhaustive(problem, budget)
    early_away = certify(problem, optimum, budget, iterations=60)
    early_plain = certify(problem, optimum, budget, iterations=60,
                          away_steps=False)
    assert early_away.relaxed_gap < early_plain.relaxed_gap
    assert early_away.lower_bound < early_plain.lower_bound
    for result in (early_away, early_plain):
        assert result.lower_bound <= optimum * (1.0 + 1e-9)


def test_a_solved_relaxation_still_leaves_an_integrality_gap(problem) -> None:
    """The certificate's ceiling is the relaxation, not the solver.

    Once the away-step iteration has converged, the bound is exactly the
    relaxed optimum, and on this instance that sits well below the true
    integer optimum. No amount of further computation closes the remainder:
    it is the price of relaxing an integer schedule to a convex mixture, and
    it is what the naming rule of D1/D29 is really testing.
    """
    budget = _mid_budget(problem)
    optimum, _ = _exhaustive(problem, budget)
    solved = certify(problem, optimum, budget, iterations=800)
    assert solved.relaxed_gap <= 1e-9 * optimum
    assert solved.lower_bound < optimum
    # The best gap this relaxation can ever certify on this instance.
    assert solved.gap_error == pytest.approx(
        1.0 - np.sqrt(solved.lower_bound / optimum), rel=1e-12)


def test_a_loose_relaxation_is_named_rather_than_read_as_slow(
        problem) -> None:
    """The predicate that decides whether more compute could ever help.

    Tested on constructed records rather than on a run: the point is the
    judgement it encodes, and a run that happened to exhibit both cases would
    be a coincidence rather than a test.
    """
    budget = _mid_budget(problem)
    optimum, _ = _exhaustive(problem, budget)
    solved = certify(problem, optimum, budget, iterations=60)
    assert not solved.vacuous
    assert not solved.structurally_vacuous

    # Bound at zero, relaxation essentially converged: no run is long enough.
    loose = replace(solved, lower_bound=0.0, vacuous=True,
                    relaxed_gap=1e-6 * optimum)
    assert loose.structurally_vacuous
    # Bound at zero, relaxation nowhere near converged: run it longer.
    slow = replace(solved, lower_bound=0.0, vacuous=True,
                   relaxed_gap=10.0 * optimum)
    assert not slow.structurally_vacuous


def test_the_linearised_costs_match_a_direct_sum(problem) -> None:
    """The away step scores stored vertices with these, so they must be exact."""
    schedule = np.full(INTERVALS, DEGREES[1])
    gradient = problem.kernel.gradient(problem.gather(schedule))
    coefficients = objective_coefficients(problem, gradient)

    expected = np.zeros_like(coefficients)
    for i in range(problem.contributions.shape[0]):
        q = int(problem.interval_of[i])
        for p in range(len(DEGREES)):
            expected[q, p] += float(problem.contributions[i, p] @ gradient[i])
    assert coefficients == pytest.approx(expected, rel=1e-12, abs=1e-15)


def _staircase():
    """A criterion whose multiplier sweep provably stops short of the ceiling.

    Four intervals, three candidates, unit weights. The middle candidate is
    only slightly better than the cheapest, the expensive one much better --
    so every multiplier that admits the expensive candidate admits it
    everywhere and overruns, and the bisection falls back to a schedule with
    room to spare.
    """
    values = np.array([[10.0, 9.5, 1.0]] * 4)
    cost = np.array([[1.0, 4.0, 9.0]] * 4)
    return values, cost


def test_the_fill_spends_room_the_multiplier_sweep_cannot_reach() -> None:
    values, cost = _staircase()
    budget = 21.0                      # 9 + 9 + 1 + 1 fits; 9*3 does not
    swept = np.zeros(4, dtype=np.int64)
    filled, moves = greedy_fill(values, cost, swept.copy(), budget)

    assert moves > 0
    work = float(cost[np.arange(4), filled].sum())
    assert work <= budget
    assert work > float(cost[np.arange(4), swept].sum())
    assert (float(values[np.arange(4), filled].sum())
            < float(values[np.arange(4), swept].sum()))


def test_the_fill_never_makes_the_criterion_worse() -> None:
    """It is marginal analysis, not a rule that empties the budget.

    The criterion is not monotone in degree --- the omitted tail is the norm
    of a partial sum and dropping a band can leave a larger residual --- so a
    fill that raised degrees until the money ran out would sometimes buy a
    worse comparator, silently.
    """
    rng = np.random.default_rng(11)
    for _ in range(40):
        values = rng.uniform(0.1, 10.0, size=(6, 4))
        cost = np.sort(rng.uniform(0.5, 6.0, size=(6, 4)), axis=1)
        picked = rng.integers(0, 4, size=6)
        budget = float(cost[np.arange(6), picked].sum()) * 1.5
        before = float(values[np.arange(6), picked].sum())
        filled, _ = greedy_fill(values, cost, picked.copy(), budget)
        assert float(values[np.arange(6), filled].sum()) <= before + 1e-12
        assert float(cost[np.arange(6), filled].sum()) <= budget * (1 + 1e-12)


def test_a_schedule_at_its_best_returns_unchanged_with_zero_moves() -> None:
    """Zero moves and no fill attempted have to stay distinguishable."""
    values = np.array([[1.0, 5.0, 9.0]] * 3)      # cheapest is also best
    cost = np.array([[1.0, 4.0, 9.0]] * 3)
    picked = np.zeros(3, dtype=np.int64)
    filled, moves = greedy_fill(values, cost, picked.copy(), 100.0)
    assert moves == 0
    assert filled.tolist() == picked.tolist()


def test_the_fill_refuses_an_infeasible_starting_schedule() -> None:
    values, cost = _staircase()
    with pytest.raises(ValueError, match="the fill improves a feasible point"):
        greedy_fill(values, cost, np.full(4, 2, dtype=np.int64), 10.0)


def test_the_separable_solver_reports_what_the_fill_did(problem) -> None:
    budget = _mid_budget(problem)
    values = force_values(np.asarray(problem.contributions[:, :, :3]),
                          np.ones(problem.contributions.shape[0]),
                          problem.interval_of, INTERVALS)
    plain = solve_separable(values, problem.time_weight, DEGREES, budget,
                            fill=False)
    filled = solve_separable(values, problem.time_weight, DEGREES, budget)

    assert "fill_moves" not in plain.diagnostics
    assert "fill_moves" in filled.diagnostics
    assert filled.objective <= plain.objective + 1e-12
    assert filled.work <= budget * (1.0 + 1e-9)


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
