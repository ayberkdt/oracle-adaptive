"""Verification of the deployable controller.

Ordered by how badly a silent failure would hurt.

* The online score is the coordinate update *exactly*. Checked against the
  kernel's own objective on an instance small enough to enumerate: for every
  candidate, the score and the full schedule's :math:`J` differ by a constant.
  A scorer that is merely plausible --- right shape, right sign, right
  ordering most of the time --- would move the controller off the benchmark
  for reasons no comparison could then attribute.
* The cancellation target excludes the interval being decided. Including it
  double-counts the block's own planned contribution, changes no norm and no
  shape, and biases every decision.
* The band stack is shared. The measured synthesis count is the size of the
  union over the window, not the sum over candidates; if it were the sum, the
  probe overhead in the manuscript would be wrong by the window size.
* The phase index anchors at apolune, so a perilune passage lands near phase
  one half rather than on the wrap.
* The budget feedback compares realized work against a realized profile. The
  two profile bases are pinned to disagree, because that disagreement is the
  previous campaign's leak and it must stay visible.
"""

from __future__ import annotations

import numpy as np
import pytest

from conftest import MU_MOON, R_MOON
from tda.allocate import DescentProblem, schedule_work
from tda.controller import (
    BudgetFeedback,
    OfflinePlan,
    OnlineController,
    OnlineSettings,
    PhaseIndexError,
    ReferenceProfile,
    RevolutionIndex,
    WorkTracker,
    apsis_epochs,
    assign_probe_points,
    build_plan,
    cancellation_target,
    colocated_cost_fraction,
    osculating_period,
    overspend_ratio,
    probe_fractions,
    radial_rate,
    score_candidates,
    select_candidate,
)
from tda.kepler import propagate_two_body
from tda.kernel import CouplingKernel
from tda.probe import required_degrees

DEGREES = (10, 20, 30)
CELLS_PER_INTERVAL = 3
INTERVALS = 4
M_CELLS = INTERVALS * CELLS_PER_INTERVAL


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


class _Table:
    """The three arrays :meth:`DescentProblem.from_table` reads."""

    def __init__(self, defect, node_transport, widths, degrees):
        self.defect = defect
        self.node_transport = node_transport
        self.widths = widths
        self.schema = type("S", (), {"candidate_degrees": degrees})()


class _Stack:
    """A band stack with a deterministic field and an honest counter.

    The values are arbitrary; what matters is that ``a_{<=n}`` converges in
    ``n`` so that consecutive bands are non-zero, and that the counter records
    one synthesis per requested degree, as
    :class:`tda.spectrum.DifferencingBandStack` does.
    """

    def __init__(self):
        self.total = 0
        self.requests = []

    def cumulative(self, r, t, degrees):
        degrees = np.asarray(degrees, dtype=int)
        self.total += degrees.size
        self.requests.append((np.asarray(r, dtype=float).copy(), float(t),
                              degrees.copy()))
        direction = np.asarray(r, dtype=float)
        direction = direction / np.linalg.norm(direction)
        swirl = np.array([np.cos(t * 1.0e-3), np.sin(t * 1.0e-3), 0.3])
        base = direction + 0.2 * swirl
        return np.stack([base * (1.0 - 1.0 / (1.0 + float(n)))
                         for n in degrees])

    def syntheses_for(self, n_degrees):
        return int(n_degrees)

    @property
    def total_syntheses(self):
        return self.total


@pytest.fixture
def pieces():
    """A tiny instance, plus the arrays the plan needs alongside it."""
    rng = np.random.default_rng(20260809)
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

    table = _Table(defect, node_transport, widths, DEGREES)
    problem = DescentProblem.from_table(table, kernel, interval_of,
                                        time_weight)
    # S_i = dt_i * Phi(t_0, m_i) B: the width folded in, so that
    # contributions[i, p] == transport[i] @ defect[i, p].
    transport = node_transport * widths[:, None, None]

    # Two revolutions over the arc; the values only have to be increasing.
    edges = np.concatenate([[0.0], np.cumsum(widths)])
    span = edges[-1]
    coordinate = 2.0 * 0.5 * (edges[:-1] + edges[1:]) / span
    edge_coordinate = 2.0 * edges[::CELLS_PER_INTERVAL] / span
    return {
        "problem": problem,
        "defect": defect,
        "transport": transport,
        "node_transport": node_transport,
        "widths": widths,
        "coordinate": coordinate,
        "edge_coordinate": edge_coordinate,
        "budget": 0.75 * schedule_work(np.full(INTERVALS, DEGREES[-1]),
                                       time_weight),
    }


@pytest.fixture
def plan(pieces) -> OfflinePlan:
    """A plan built from the tiny instance, single-start as declared."""
    starts = [np.full(INTERVALS, DEGREES[0], dtype=np.int64)]
    return build_plan(pieces["problem"], pieces["node_transport"],
                      pieces["widths"], pieces["coordinate"],
                      pieces["edge_coordinate"], pieces["budget"], starts,
                      iterations=3)


def _elliptic_arc(n_samples=4001, revolutions=3.0):
    """A sampled Keplerian ellipse: perilune at t=0, apolune at half a period."""
    perilune = R_MOON + 5.0e4
    apolune = R_MOON + 2.0e6
    axis = 0.5 * (perilune + apolune)
    eccentricity = (apolune - perilune) / (apolune + perilune)
    speed = np.sqrt(MU_MOON / axis * (1.0 + eccentricity)
                    / (1.0 - eccentricity))
    state = np.array([perilune, 0.0, 0.0, 0.0, speed, 0.0])
    period = 2.0 * np.pi * np.sqrt(axis**3 / MU_MOON)
    times = np.linspace(0.0, revolutions * period, n_samples)
    states = np.stack([propagate_two_body(state, MU_MOON, float(t))
                       for t in times])
    return times, states, period


# ---------------------------------------------------------------------------
# Phase index
# ---------------------------------------------------------------------------


def test_radial_rate_vanishes_on_a_circular_orbit(circular_state):
    assert abs(radial_rate(circular_state[None, :])[0]) < 1e-9


def test_radial_rate_is_the_speed_in_purely_radial_motion():
    state = np.array([3.0, 4.0, 0.0, 0.6, 0.8, 0.0])   # v parallel to r
    assert radial_rate(state[None, :])[0] == pytest.approx(1.0)


def test_apsis_epoch_is_interpolated_not_snapped():
    times = np.arange(6.0)
    rate = times - 3.7                       # zero at 3.7, between samples
    assert apsis_epochs(times, rate, falling=False)[0] == pytest.approx(3.7)
    assert apsis_epochs(times, -rate, falling=True)[0] == pytest.approx(3.7)


def test_a_sample_exactly_at_zero_counts_once():
    times = np.arange(5.0)
    rate = np.array([1.0, 0.0, -1.0, -2.0, -3.0])
    assert apsis_epochs(times, rate, falling=True) == pytest.approx([1.0])


def test_osculating_period_matches_the_closed_form():
    axis = R_MOON + 3.0e5
    speed = np.sqrt(MU_MOON / axis)
    state = np.array([axis, 0.0, 0.0, 0.0, speed, 0.0])
    expected = 2.0 * np.pi * np.sqrt(axis**3 / MU_MOON)
    assert osculating_period(state, MU_MOON) == pytest.approx(expected)


def test_an_unbound_state_has_no_revolutions():
    radius = R_MOON + 1.0e5
    escape = np.sqrt(2.0 * MU_MOON / radius)
    state = np.array([radius, 0.0, 0.0, 0.0, 1.1 * escape, 0.0])
    with pytest.raises(PhaseIndexError, match="unbound"):
        osculating_period(state, MU_MOON)


def test_perilune_lands_at_mid_phase():
    """The design point: anchoring at apolune keeps perilune off the wrap."""
    times, states, period = _elliptic_arc()
    index = RevolutionIndex.from_arc(times, states, MU_MOON)

    assert len(index) == 2
    assert index.anchors[0] == pytest.approx(0.5 * period, rel=1e-4)
    assert index.spans == pytest.approx([period, period], rel=1e-4)
    # Perilune of the second revolution, at t = 2P.
    assert float(index.coordinate(2.0 * period)) == pytest.approx(1.5,
                                                                  abs=1e-3)


def test_coordinate_and_epoch_invert_each_other():
    times, states, period = _elliptic_arc()
    index = RevolutionIndex.from_arc(times, states, MU_MOON)
    probe = np.array([-0.3 * period, 0.9 * period, 2.4 * period])
    assert index.epoch(index.coordinate(probe)) == pytest.approx(probe)


def test_the_coordinate_stays_monotone_through_the_partial_revolutions():
    times, states, period = _elliptic_arc()
    index = RevolutionIndex.from_arc(times, states, MU_MOON)
    sweep = index.coordinate(np.linspace(-period, 4.0 * period, 500))
    assert np.all(np.diff(sweep) > 0.0)
    assert sweep[0] < 0.0                     # before the first anchor
    assert not index.covers(sweep[0])
    assert index.covers(1.0)


def test_a_near_circular_arc_is_refused_rather_than_indexed_from_noise():
    radius = R_MOON + 1.0e5
    speed = np.sqrt(MU_MOON / radius)
    state = np.array([radius, 0.0, 0.0, 0.0, speed, 0.0])
    period = 2.0 * np.pi * np.sqrt(radius**3 / MU_MOON)
    times = np.linspace(0.0, 3.0 * period, 601)
    states = np.stack([propagate_two_body(state, MU_MOON, float(t))
                       for t in times])
    with pytest.raises(PhaseIndexError, match="harmonics rather than"):
        RevolutionIndex.from_arc(times, states, MU_MOON)


def test_spurious_apsides_are_dropped_and_counted():
    times, states, period = _elliptic_arc()
    assert apsis_epochs(times, radial_rate(states), falling=True).size == 3

    # Reverse the velocity over three samples just after the first apolune:
    # the radial rate goes negative -> positive -> negative and manufactures a
    # fourth falling crossing a few samples later.
    first = int(np.argmin(np.abs(times - 0.5 * period)))
    states = states.copy()
    states[first + 2:first + 5, 3:6] *= -1.0
    assert apsis_epochs(times, radial_rate(states), falling=True).size == 4

    index = RevolutionIndex.from_arc(times, states, MU_MOON)
    assert index.spurious == 1
    assert len(index) == 2
    assert index.anchors[0] == pytest.approx(0.5 * period, rel=1e-4)


def test_an_arc_shorter_than_a_revolution_cannot_be_indexed():
    times, states, _ = _elliptic_arc(n_samples=400, revolutions=0.9)
    with pytest.raises(PhaseIndexError, match="shorter than one revolution"):
        RevolutionIndex.from_arc(times, states, MU_MOON)


# ---------------------------------------------------------------------------
# The plan
# ---------------------------------------------------------------------------


def _dense_q(kernel):
    """:math:`\\mathbf Q_{ik}=\\mathbf A_{\\max(i,k)}/T`, formed explicitly."""
    n = len(kernel)
    out = np.empty((n, n, 6, 6))
    for i in range(n):
        for k in range(n):
            out[i, k] = kernel.suffix[max(i, k)] / kernel.duration
    return out


def test_cancellation_target_matches_its_definition(pieces, plan):
    """Brute force against ``z_i = S_i^T sum_{k outside I_q} Q_ik u_k``."""
    problem = pieces["problem"]
    u = problem.gather(plan.nominal_degrees)
    dense = _dense_q(problem.kernel)
    expected = np.empty((M_CELLS, 3))
    for i in range(M_CELLS):
        q = int(problem.interval_of[i])
        outside = problem.interval_of != q
        c = np.einsum("kab,kb->a", dense[i, outside], u[outside])
        expected[i] = pieces["transport"][i].T @ c
    assert plan.target == pytest.approx(expected, rel=1e-10, abs=1e-14)
    # And the function reproduces what build_plan stored.
    direct = cancellation_target(problem.kernel, pieces["transport"], u,
                                 problem.cell_slices, problem.interval_of)
    assert direct == pytest.approx(plan.target, rel=1e-12, abs=1e-15)


def test_forgetting_the_exclusion_would_change_the_decision(pieces, plan):
    """The block's own contribution must not appear in its own target."""
    problem = pieces["problem"]
    u = problem.gather(plan.nominal_degrees)
    unexcluded = np.einsum("iba,ib->ia", pieces["transport"],
                           problem.kernel.cross_terms(u))
    assert not np.allclose(unexcluded, plan.target, rtol=1e-6, atol=1e-12)

    leaky = OfflinePlan(**{**{f: getattr(plan, f) for f in plan.__slots__},
                           "target": unexcluded})
    changed = 0
    for q in range(len(plan)):
        lo, hi = plan.cell_slices[q]
        estimates = pieces["defect"][lo:hi]
        honest = score_candidates(plan, q, estimates, plan.multiplier, DEGREES)
        wrong = score_candidates(leaky, q, estimates, plan.multiplier, DEGREES)
        # The damage is not a constant offset: it reorders the candidates.
        offset = honest - wrong
        assert offset.max() - offset.min() > 1e-12
        changed += select_candidate(honest, DEGREES) != select_candidate(
            wrong, DEGREES)
    assert changed > 0


def test_the_plan_records_whether_another_sweep_would_move_anything(pieces):
    starts = [np.full(INTERVALS, DEGREES[0], dtype=np.int64)]
    converged = build_plan(pieces["problem"], pieces["node_transport"],
                           pieces["widths"], pieces["coordinate"],
                           pieces["edge_coordinate"], pieces["budget"],
                           starts, iterations=8)
    assert converged.converged
    assert converged.iterations == 8


def test_the_plan_declares_where_its_information_came_from(pieces):
    starts = [np.full(INTERVALS, DEGREES[0], dtype=np.int64)]
    with pytest.raises(ValueError, match="source must be one of"):
        build_plan(pieces["problem"], pieces["node_transport"],
                   pieces["widths"], pieces["coordinate"],
                   pieces["edge_coordinate"], pieces["budget"], starts,
                   source="whatever")


def test_a_non_monotone_phase_coordinate_is_refused(pieces):
    starts = [np.full(INTERVALS, DEGREES[0], dtype=np.int64)]
    broken = pieces["coordinate"].copy()
    broken[3], broken[4] = broken[4], broken[3]
    with pytest.raises(ValueError, match="strictly increasing"):
        build_plan(pieces["problem"], pieces["node_transport"],
                   pieces["widths"], broken, pieces["edge_coordinate"],
                   pieces["budget"], starts)


def test_zero_iterations_is_not_a_plan(pieces):
    starts = [np.full(INTERVALS, DEGREES[0], dtype=np.int64)]
    with pytest.raises(ValueError, match="iterations must be positive"):
        build_plan(pieces["problem"], pieces["node_transport"],
                   pieces["widths"], pieces["coordinate"],
                   pieces["edge_coordinate"], pieces["budget"], starts,
                   iterations=0)


def test_a_transport_from_another_arc_is_caught(pieces):
    """It would leave every shape right and every number finite."""
    starts = [np.full(INTERVALS, DEGREES[0], dtype=np.int64)]
    stranger = np.random.default_rng(7).normal(size=(M_CELLS, 6, 3))
    with pytest.raises(ValueError, match="column space"):
        build_plan(pieces["problem"], stranger, pieces["widths"],
                   pieces["coordinate"], pieces["edge_coordinate"],
                   pieces["budget"], starts)


def test_the_cell_width_is_folded_by_the_plan_not_by_the_caller(pieces, plan):
    """A missing width passes every structural check, so it is not offered."""
    assert plan.transport == pytest.approx(pieces["transport"])
    assert not np.allclose(plan.transport, pieces["node_transport"])


def test_the_plan_obeys_the_ceiling(plan):
    assert plan.work <= plan.budget
    assert plan.multiplier >= 0.0
    assert plan.nbytes > 0


def test_drifting_past_the_plan_raises_rather_than_clamping(plan):
    with pytest.raises(IndexError, match="outside the plan's coverage"):
        plan.locate(float(plan.edge_coordinate[-1]) + 1e-6)
    assert plan.locate(float(plan.edge_coordinate[-1])) == len(plan) - 1
    assert plan.locate(float(plan.edge_coordinate[0])) == 0


def test_cell_fractions_partition_the_interval(plan):
    for q in range(len(plan)):
        fractions = plan.cell_fractions(q)
        assert fractions.size == CELLS_PER_INTERVAL
        assert np.all(np.diff(fractions) > 0.0)
        assert fractions[0] > 0.0
        assert fractions[-1] < 1.0


# ---------------------------------------------------------------------------
# The online score: the axiom
# ---------------------------------------------------------------------------


def test_the_online_score_is_the_objective_up_to_a_constant(pieces, plan):
    """Every candidate, every interval, against the kernel's own ``J``."""
    problem = pieces["problem"]
    multiplier = 3.0e-4
    for q in range(len(plan)):
        lo, hi = plan.cell_slices[q]
        scores = score_candidates(plan, q, pieces["defect"][lo:hi],
                                  multiplier, DEGREES)
        penalized = []
        for degree in DEGREES:
            schedule = plan.nominal_degrees.copy()
            schedule[q] = degree
            penalized.append(
                problem.kernel.objective(problem.gather(schedule))
                + multiplier * schedule_work(schedule, problem.time_weight))
        scale = max(1.0, float(np.abs(scores).max()))
        # The candidates really do separate, so the identity is not holding
        # because everything collapsed to one number.
        assert scores.max() - scores.min() > 1.0
        offset = scores - np.array(penalized)
        assert offset.max() - offset.min() < 1e-12 * scale


def test_the_plan_is_a_fixed_point_of_its_own_online_decision(pieces, plan):
    """Given the exact defect, the online rule reproduces the offline plan.

    The two halves are written independently --- one sweeps the whole arc, the
    other sees a single block and a pre-contracted target --- so agreeing on
    every interval is evidence that the exclusion, the widths, the
    :math:`1/T` and the tie rule survived the contraction. It holds exactly
    because the plan converged; a plan capped short of convergence is not
    required to be its own fixed point.
    """
    assert plan.converged
    for q in range(len(plan)):
        lo, hi = plan.cell_slices[q]
        scores = score_candidates(plan, q, pieces["defect"][lo:hi],
                                  plan.multiplier, DEGREES)
        assert select_candidate(scores, DEGREES) == plan.nominal_degrees[q]


def test_the_score_rejects_an_estimate_tensor_of_the_wrong_shape(plan):
    with pytest.raises(ValueError, match="estimates must have shape"):
        score_candidates(plan, 0, np.zeros((2, 3, 3)), 0.0, DEGREES)


def test_a_tie_goes_to_the_cheaper_degree():
    assert select_candidate(np.array([1.0, 1.0, 1.0]), (10, 20, 30)) == 10
    assert select_candidate(np.array([2.0, 1.0, 1.0]), (10, 20, 30)) == 20


# ---------------------------------------------------------------------------
# Probe geometry and cost
# ---------------------------------------------------------------------------


def test_probe_points_avoid_the_boundary():
    fractions = probe_fractions(3, "midpoint")
    assert fractions == pytest.approx([1 / 6, 0.5, 5 / 6])
    assert fractions[0] > 0.0
    assert fractions[-1] < 1.0


def test_the_co_located_variant_takes_exactly_one_point():
    assert probe_fractions(1, "leading") == pytest.approx([0.0])
    with pytest.raises(ValueError, match="exactly one point"):
        probe_fractions(2, "leading")


def test_an_undeclared_placement_is_refused():
    with pytest.raises(ValueError, match="placement must be one of"):
        probe_fractions(1, "wherever")


def test_cells_take_the_nearest_probe_point():
    cells = np.array([0.05, 0.3, 0.5, 0.95])
    points = np.array([0.25, 0.75])
    assert assign_probe_points(cells, points).tolist() == [0, 0, 0, 1]
    # A single point serves the whole interval: this is what C-lite does. The
    # cell sitting exactly on the point is the case a bracket search gets
    # wrong, so it is in the list.
    assert assign_probe_points(np.array([0.0, *cells]),
                               np.array([0.0])).tolist() == [0] * 5


def test_the_incremental_cost_is_the_manuscripts_formula():
    assert colocated_cost_fraction(120, 0, 3) == pytest.approx(0.05)
    # A wider window costs more, and it is the span in degrees that counts.
    assert colocated_cost_fraction(120, 20, 3) == pytest.approx(2 * 23 / 120)
    with pytest.raises(ValueError, match="span must be non-negative"):
        colocated_cost_fraction(120, -1, 3)


# ---------------------------------------------------------------------------
# The controller end to end
# ---------------------------------------------------------------------------


def _sigma_a(_radius):
    """A decaying spectrum, long enough for the deepest probe used here."""
    n = np.arange(64, dtype=float)
    return 1.0 / (n + 1.0) ** 2


def _boundary_state():
    radius = R_MOON + 1.0e5
    speed = np.sqrt(MU_MOON / radius)
    return np.array([radius, 0.0, 0.0, 0.0, 0.8 * speed, 0.4 * speed])


def test_one_shared_stack_serves_every_candidate(plan):
    stack = _Stack()
    controller = OnlineController(plan, OnlineSettings(half_width=1, depth=3,
                                                       n_probe=4))
    decision = controller.decide(0.1, _boundary_state(), 0.0, 120.0, stack,
                                 _sigma_a, MU_MOON)

    per_point = required_degrees(decision.window, 3).size
    assert decision.syntheses == 4 * per_point
    # One stack request per probe point, not one per candidate: the candidates
    # come out of it by partial summation.
    assert len(stack.requests) == 4
    assert controller.stack_depth(decision.interval) == per_point


def test_the_union_shrinks_only_when_the_depth_meets_the_spacing():
    """What the shared stack is worth on the campaign's own candidate grid.

    At spacing ten and depth three no candidate's ``N+k`` is another
    candidate's ``N``, so the union is the sum and a differencing stack pays
    two syntheses per candidate. The manuscript's stronger claim --- one
    Legendre recursion for the whole window --- needs the cumulative-by-degree
    kernel entry point that does not exist yet (``DECISIONS.md`` D120), and
    this test is what keeps the two from being quoted for one another.
    """
    assert required_degrees((10, 20, 30), 3).size == 2 * 3
    assert required_degrees((10, 20, 30), 10).size == 4


def test_the_probe_is_taken_ahead_of_the_boundary_at_its_own_radius(plan):
    seen = []

    def recording(radius):
        seen.append(radius)
        return _sigma_a(radius)

    stack = _Stack()
    controller = OnlineController(plan, OnlineSettings(n_probe=4))
    state = _boundary_state()
    decision = controller.decide(0.1, state, 0.0, 600.0, stack, recording,
                                 MU_MOON)

    assert decision.probe_points == 4
    assert len(seen) == 4
    # Distinct radii, and none of them the boundary's: the completion factor
    # is evaluated where the probe is, not where the decision is taken.
    assert len(set(seen)) == 4
    assert min(abs(r - np.linalg.norm(state[0:3])) for r in seen) > 1.0
    epochs = [t for _, t, _ in stack.requests]
    assert min(epochs) > 0.0
    assert max(epochs) < 600.0


def test_c_lite_uses_no_predictor_and_says_so(plan):
    stack = _Stack()
    controller = OnlineController(
        plan, OnlineSettings(n_probe=1, placement="leading"))
    state = _boundary_state()
    decision = controller.decide(0.1, state, 17.0, 120.0, stack, _sigma_a,
                                 MU_MOON, perturbing_acceleration=1e-3)

    assert decision.co_located
    assert decision.probe_points == 1
    assert decision.predictor_error_bound == 0.0
    position, epoch, _ = stack.requests[0]
    assert position == pytest.approx(state[0:3])
    assert epoch == 17.0


def test_the_predictor_bound_is_carried_not_recomputed(plan):
    stack = _Stack()
    controller = OnlineController(plan, OnlineSettings(n_probe=2))
    decision = controller.decide(0.1, _boundary_state(), 0.0, 120.0, stack,
                                 _sigma_a, MU_MOON,
                                 perturbing_acceleration=2.0e-3)
    # Half a milli-g over the last probe step, three quarters of the interval.
    assert decision.predictor_error_bound == pytest.approx(
        0.5 * 2.0e-3 * (0.75 * 120.0) ** 2)


def test_a_clipped_window_does_not_report_a_spurious_move(plan):
    """Near the ends of the grid the window's middle is not the nominal."""
    edge = OfflinePlan(**{**{f: getattr(plan, f) for f in plan.__slots__},
                          "nominal_degrees": np.full(len(plan), DEGREES[0])})
    stack = _Stack()
    controller = OnlineController(edge, OnlineSettings(half_width=2))
    decision = controller.decide(0.1, _boundary_state(), 0.0, 120.0, stack,
                                 _sigma_a, MU_MOON)
    assert decision.window == (10, 20, 30)
    assert decision.nominal == 10
    assert decision.moved == (decision.degree != 10)


def test_a_decision_outside_the_plan_is_an_error(plan):
    controller = OnlineController(plan)
    with pytest.raises(IndexError):
        controller.decide(99.0, _boundary_state(), 0.0, 120.0, _Stack(),
                          _sigma_a, MU_MOON)


def test_an_open_loop_controller_refuses_to_be_charged(plan):
    controller = OnlineController(plan)
    with pytest.raises(RuntimeError, match="open loop"):
        controller.charge(30)
    assert controller.multiplier_at(0.1) == plan.multiplier


# ---------------------------------------------------------------------------
# Budget feedback
# ---------------------------------------------------------------------------


def test_the_tracker_counts_calls_and_squared_degree():
    tracker = WorkTracker()
    tracker.charge(10, 3)
    tracker.charge(20)
    assert tracker.spent == pytest.approx(3 * 100 + 400)
    assert tracker.calls == 4
    assert tracker.fraction(1400.0) == pytest.approx(0.5)
    with pytest.raises(ValueError, match="n_calls must be positive"):
        tracker.charge(10, 0)


def test_the_two_profile_bases_disagree_when_the_calls_are_not_uniform():
    """The previous campaign's leak, made visible rather than assumed away."""
    # One decision interval of each degree, equal time weight, but the
    # integrator takes four times as many steps in the high-degree one.
    edge_coordinate = np.array([0.0, 0.5, 1.0])
    nominal = ReferenceProfile.from_nominal_schedule(
        edge_coordinate, np.array([1.0, 1.0]), np.array([10.0, 20.0]))
    calls = np.concatenate([np.linspace(0.0, 0.5, 4, endpoint=False),
                            np.linspace(0.5, 1.0, 16, endpoint=False)])
    degrees = np.concatenate([np.full(4, 10.0), np.full(16, 20.0)])
    realized = ReferenceProfile.from_realized_calls(calls, degrees)

    assert nominal.basis == "nominal"
    assert realized.basis == "realized"
    # By the halfway point the nominal profile says a fifth of the budget is
    # gone; the realized one says an eighth. Driving on the nominal profile
    # would command a slowdown against an overspend that is not happening.
    assert nominal.at(0.5) == pytest.approx(100 / 500)
    assert realized.at(0.5) == pytest.approx(800 / 6800)
    assert nominal.at(0.5) > 1.5 * realized.at(0.5)


def test_calls_sharing_an_epoch_are_collapsed_not_rejected():
    """A multi-stage integrator can evaluate twice at the same node."""
    profile = ReferenceProfile.from_realized_calls(
        np.array([0.0, 0.0, 0.5, 1.0]), np.array([10.0, 10.0, 10.0, 10.0]))
    assert profile.coordinate == pytest.approx([0.0, 0.5, 1.0])
    assert profile.fraction == pytest.approx([0.5, 0.75, 1.0])
    with pytest.raises(ValueError, match="non-decreasing"):
        ReferenceProfile.from_realized_calls(np.array([0.0, 1.0, 0.5]),
                                             np.array([10.0, 10.0, 10.0]))


def test_a_profile_that_falls_is_refused():
    with pytest.raises(ValueError, match="non-decreasing"):
        ReferenceProfile(np.array([0.0, 1.0]), np.array([0.5, 0.2]),
                         "realized")
    with pytest.raises(ValueError, match="basis must be one of"):
        ReferenceProfile(np.array([0.0, 1.0]), np.array([0.0, 1.0]), "made up")


def _feedback(**kwargs):
    profile = ReferenceProfile(np.array([0.0, 1.0]), np.array([0.0, 1.0]),
                               "realized")
    settings = {"ceiling": 1000.0, "nominal_multiplier": 2.0,
                "reference_multiplier": 5.0, "profile": profile}
    settings.update(kwargs)
    return BudgetFeedback(**settings)


def test_on_plan_spending_leaves_the_multiplier_alone():
    feedback = _feedback()
    feedback.charge(10, 5)                      # 500 of 1000, at phase 0.5
    assert feedback.error(0.5) == pytest.approx(0.0)
    assert feedback.multiplier_at(0.5) == pytest.approx(2.0)


def test_an_overspend_raises_the_price_and_an_underspend_lowers_it():
    feedback = _feedback(gain=1.0)
    feedback.charge(10, 8)                      # 800 of 1000, at phase 0.5
    assert feedback.multiplier_at(0.5) == pytest.approx(2.0 + 5.0 * 0.3)
    lenient = _feedback(gain=1.0)
    lenient.charge(10, 2)
    assert lenient.multiplier_at(0.5) == pytest.approx(2.0 - 5.0 * 0.3)


def test_the_loop_can_raise_the_price_from_a_slack_ceiling():
    """A multiplicative law could not; ``lambda_0 = 0`` is legitimate (D142)."""
    feedback = _feedback(nominal_multiplier=0.0, gain=1.0)
    feedback.charge(10, 8)
    assert feedback.multiplier_at(0.5) == pytest.approx(5.0 * 0.3)


def test_the_price_never_goes_negative():
    feedback = _feedback(nominal_multiplier=0.1, gain=10.0)
    assert feedback.multiplier_at(1.0) == 0.0   # nothing spent, all planned


def test_zero_gain_is_the_open_loop_control():
    feedback = _feedback(gain=0.0)
    feedback.charge(30, 100)
    assert feedback.multiplier_at(0.1) == pytest.approx(2.0)


def test_saturation_is_counted_rather_than_hidden():
    feedback = _feedback(gain=1.0e9, max_factor=10.0)
    feedback.charge(30, 100)
    assert feedback.multiplier_at(0.0) == pytest.approx(50.0)
    assert feedback.saturations == 1


def test_a_negative_gain_is_refused():
    with pytest.raises(ValueError, match="gain must be non-negative"):
        _feedback(gain=-1.0)
    with pytest.raises(ValueError, match="reference_multiplier must be"):
        _feedback(reference_multiplier=0.0)


def test_the_controller_prices_through_its_feedback_loop(plan):
    feedback = _feedback(gain=1.0, nominal_multiplier=plan.multiplier)
    controller = OnlineController(plan, OnlineSettings(n_probe=2),
                                  feedback=feedback)
    controller.charge(30, 30)                   # 27000 against a 1000 ceiling
    decision = controller.decide(0.1, _boundary_state(), 0.0, 120.0, _Stack(),
                                 _sigma_a, MU_MOON)
    assert decision.multiplier > plan.multiplier
    assert decision.multiplier == pytest.approx(feedback.multiplier_at(0.1))


def test_the_overspend_ratio_is_one_when_the_calls_are_uniform():
    assert overspend_ratio(400.0, 100, 200.0, 50.0) == pytest.approx(1.0)
    # Twice as many calls where the degree is doubled: the leak.
    spent = 100 * 10.0**2 + 200 * 20.0**2
    nominal = 50.0 * 10.0**2 + 50.0 * 20.0**2
    assert overspend_ratio(spent, 300, nominal, 100.0) > 1.0
