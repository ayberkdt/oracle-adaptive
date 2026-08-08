"""Verification of the accumulation and decision grids.

Two properties carry the campaign and are asserted rather than assumed: the
inner weights sum to the arc length (so the budget constraint means what it
says), and the refinement actually equidistributes the correlation time (so
the convergence parameter :math:`n_s` is the quantity it is swept as).
"""

from __future__ import annotations

import numpy as np
import pytest

from tda.config import GridConfig
from tda.grids import (
    build_accumulation_grid,
    build_decision_grid,
    correlation_time,
    trapezoid_weights,
)

DURATION = 7200.0


def _eccentric_sample(n: int = 4001) -> tuple[np.ndarray, np.ndarray,
                                              np.ndarray]:
    """A sampled arc whose radius swings by a factor of three.

    Not a real orbit: an analytic stand-in with the property that matters
    here, namely that :math:`\\tau_{\\mathrm{corr}}` varies strongly along the
    arc.  Using a real propagation would test the propagator, not the grid.
    """
    t = np.linspace(0.0, DURATION, n)
    phase = 2.0 * np.pi * t / DURATION
    radius = 1.8e6 + 1.8e6 * (1.0 - np.cos(phase))     # 1.8e6 .. 5.4e6
    speed = np.sqrt(4.9028e12 * (2.0 / radius - 1.0 / 3.6e6))
    return t, radius, speed


# ---------------------------------------------------------------------------
# Correlation time
# ---------------------------------------------------------------------------


def test_correlation_time_matches_the_definition() -> None:
    r = np.array([1.8374e6])
    v = np.array([1633.0])
    assert correlation_time(r, v, 300)[0] == pytest.approx(
        np.pi * 1.8374e6 / (300 * 1633.0))


def test_correlation_time_is_inverse_in_degree() -> None:
    """Doubling the degree halves the scale the signed integral must resolve."""
    r, v = np.array([1.8e6]), np.array([1600.0])
    assert correlation_time(r, v, 600)[0] == pytest.approx(
        0.5 * correlation_time(r, v, 300)[0])


@pytest.mark.parametrize(("radius", "speed", "degree", "message"), [
    (np.array([1.8e6]), np.array([1600.0]), 0, "degree must be positive"),
    (np.array([0.0]), np.array([1600.0]), 300, "strictly positive"),
    (np.array([1.8e6]), np.array([0.0]), 300, "strictly positive"),
    (np.array([1.8e6, 1.0]), np.array([1600.0]), 300, "must match"),
])
def test_correlation_time_rejects_degenerate_input(radius, speed, degree,
                                                   message) -> None:
    with pytest.raises(ValueError, match=message):
        correlation_time(radius, speed, degree)


# ---------------------------------------------------------------------------
# Quadrature weights
# ---------------------------------------------------------------------------


def test_trapezoid_weights_sum_to_the_span() -> None:
    nodes = np.array([0.0, 1.0, 3.0, 7.0, 7.5])
    assert trapezoid_weights(nodes).sum() == pytest.approx(7.5)


def test_trapezoid_weights_are_non_negative() -> None:
    """Non-negativity is what makes ``Q`` positive semidefinite."""
    rng = np.random.default_rng(11)
    nodes = np.concatenate([[0.0], np.sort(rng.uniform(0, 100, 40)), [100.0]])
    assert np.all(trapezoid_weights(nodes) >= 0.0)


def test_trapezoid_weights_integrate_a_linear_function_exactly() -> None:
    nodes = np.array([0.0, 0.3, 1.7, 2.0, 5.0])
    w = trapezoid_weights(nodes)
    assert float(w @ (3.0 * nodes + 1.0)) == pytest.approx(
        3.0 * 5.0**2 / 2.0 + 5.0)


@pytest.mark.parametrize(("nodes", "message"), [
    (np.array([1.0]), "at least two"),
    (np.array([0.0, 1.0, 0.5]), "strictly increasing"),
    (np.array([0.0, 0.0]), "strictly increasing"),
])
def test_trapezoid_weights_reject_a_bad_grid(nodes, message) -> None:
    with pytest.raises(ValueError, match=message):
        trapezoid_weights(nodes)


# ---------------------------------------------------------------------------
# Accumulation grid
# ---------------------------------------------------------------------------


def test_grid_spans_the_arc_exactly() -> None:
    """The inner weights are the budget's denominator; they must sum to T."""
    t, r, v = _eccentric_sample()
    grid = build_accumulation_grid(t, r, v, 300, GridConfig())

    assert grid.edges[0] == 0.0
    assert grid.edges[-1] == DURATION
    assert grid.widths.sum() == pytest.approx(DURATION, rel=1e-12)
    assert grid.outer_weights.sum() == pytest.approx(DURATION, rel=1e-12)
    assert np.all(grid.outer_weights >= 0.0)


def test_nodes_are_the_cell_midpoints() -> None:
    """The midpoint rule is what makes the inner quadrature second-order."""
    t, r, v = _eccentric_sample()
    grid = build_accumulation_grid(t, r, v, 300, GridConfig())
    assert np.allclose(grid.nodes, 0.5 * (grid.edges[:-1] + grid.edges[1:]),
                       rtol=0.0, atol=1e-9)
    assert np.all(grid.nodes > grid.edges[:-1])
    assert np.all(grid.nodes < grid.edges[1:])


def test_step_bounds_hold_to_one_part_in_m_cells() -> None:
    """Clamping the density, not the steps, is what makes the bounds hold.

    The tolerance is not a fudge: rescaling the density by
    ``ceil(S(T))/S(T)`` to land the last edge on ``T`` tightens every cell by
    at most one part in ``M``, so the lower bound can be undershot by exactly
    that and no more.  Asserting a tighter bound would be asserting something
    the construction does not provide; asserting a looser one would stop
    testing the clamp.
    """
    t, r, v = _eccentric_sample()
    cfg = GridConfig(samples_per_tau=2.0, dt_acc_min_s=5.0, dt_acc_max_s=60.0)
    grid = build_accumulation_grid(t, r, v, 300, cfg)

    slack = 1.0 + 2.0 / len(grid)
    assert grid.widths.min() >= cfg.dt_acc_min_s / slack
    assert grid.widths.max() <= cfg.dt_acc_max_s
    # The clamp must actually bind somewhere, or the test proves nothing.
    assert grid.widths.min() < 1.5 * cfg.dt_acc_min_s


def test_refinement_follows_the_correlation_time() -> None:
    """Cells must be short where ``tau_corr`` is short, not merely on average.

    The measurement is a rank correlation between cell width and the local
    correlation time; a uniform grid scores zero and would pass a test that
    only checked the total.
    """
    t, r, v = _eccentric_sample()
    cfg = GridConfig(samples_per_tau=2.0, dt_acc_min_s=0.5,
                     dt_acc_max_s=600.0)
    grid = build_accumulation_grid(t, r, v, 300, cfg)

    tau_at_nodes = np.interp(grid.nodes, t, correlation_time(r, v, 300))
    order_tau = np.argsort(np.argsort(tau_at_nodes))
    order_width = np.argsort(np.argsort(grid.widths))
    rank_corr = np.corrcoef(order_tau, order_width)[0, 1]
    assert rank_corr > 0.99

    # And the ratio the equidistribution targets is nearly constant.
    samples_per_tau = tau_at_nodes / grid.widths
    spread = samples_per_tau.max() / samples_per_tau.min()
    assert spread < 1.2


def test_halving_the_target_step_roughly_doubles_the_cells() -> None:
    t, r, v = _eccentric_sample()
    coarse = build_accumulation_grid(
        t, r, v, 300, GridConfig(samples_per_tau=2.0, dt_acc_min_s=0.5,
                                 dt_acc_max_s=600.0))
    fine = build_accumulation_grid(
        t, r, v, 300, GridConfig(samples_per_tau=4.0, dt_acc_min_s=0.5,
                                 dt_acc_max_s=600.0))
    assert len(fine) == pytest.approx(2 * len(coarse), rel=0.02)


def test_required_epochs_interleave_edges_and_nodes() -> None:
    """``t_eval`` must contain every epoch a rule samples, in order."""
    t, r, v = _eccentric_sample()
    grid = build_accumulation_grid(t, r, v, 300, GridConfig())
    epochs = grid.required_epochs()

    assert epochs.size == 2 * len(grid) + 1
    assert np.all(np.diff(epochs) > 0.0)
    assert epochs[0] == 0.0
    assert np.array_equal(epochs[grid.edge_indices()], grid.edges)
    assert np.array_equal(epochs[grid.node_indices()], grid.nodes)


@pytest.mark.parametrize(("kwargs", "message"), [
    ({"samples_per_tau": 0.0}, "samples_per_tau must be positive"),
    ({"dt_acc_min_s": 100.0, "dt_acc_max_s": 10.0}, "dt_acc_min_s <"),
])
def test_accumulation_grid_rejects_degenerate_settings(kwargs,
                                                       message) -> None:
    t, r, v = _eccentric_sample()
    with pytest.raises(ValueError, match=message):
        build_accumulation_grid(t, r, v, 300, GridConfig(**kwargs))


@pytest.mark.parametrize(("sample_t", "message"), [
    (np.array([0.0]), "at least two"),
    (np.array([1.0, 2.0]), "start at 0.0"),
    (np.array([0.0, 2.0, 1.0]), "strictly increasing"),
])
def test_accumulation_grid_rejects_a_bad_sample(sample_t, message) -> None:
    r = np.full(sample_t.size, 1.8e6)
    v = np.full(sample_t.size, 1600.0)
    with pytest.raises(ValueError, match=message):
        build_accumulation_grid(sample_t, r, v, 300, GridConfig())


# ---------------------------------------------------------------------------
# Decision grid
# ---------------------------------------------------------------------------


def test_decision_weights_partition_the_arc() -> None:
    """``sum W_q`` is ``T``: every cell is charged to exactly one interval."""
    t, r, v = _eccentric_sample()
    grid = build_accumulation_grid(t, r, v, 300, GridConfig())
    dec = build_decision_grid(grid, GridConfig())

    assert dec.time_weight.sum() == pytest.approx(DURATION, rel=1e-12)
    assert np.all(dec.time_weight > 0.0)
    assert len(dec) == int(np.ceil(DURATION / GridConfig().dt_dec_s))


def test_decision_weight_is_time_not_a_cell_count() -> None:
    """The distinction the refined grid makes expensive to get wrong.

    Where cells are dense the count per interval is large but the time weight
    is not, and charging the budget by count would bill perilune for most of
    the arc.
    """
    t, r, v = _eccentric_sample()
    cfg = GridConfig(samples_per_tau=2.0, dt_acc_min_s=0.5,
                     dt_acc_max_s=600.0, dt_dec_s=120.0)
    grid = build_accumulation_grid(t, r, v, 300, cfg)
    dec = build_decision_grid(grid, cfg)

    counts = np.bincount(dec.interval_of, minlength=len(dec))
    assert counts.max() > 3 * counts.min(), "sample is not eccentric enough"

    # A cell is charged whole to the interval its midpoint falls in, so an
    # interval's weight is its nominal length give or take one cell at each
    # end.  That is the exact bound; a tighter one would not be true, and the
    # point is that it holds while the *count* varies several-fold.
    widest = float(grid.widths.max())
    assert np.all(np.abs(dec.time_weight - cfg.dt_dec_s) <= widest + 1e-9)

    weight_spread = dec.time_weight.max() / dec.time_weight.min()
    count_spread = counts.max() / counts.min()
    assert weight_spread < count_spread / 2.0


def test_every_cell_belongs_to_exactly_one_interval() -> None:
    t, r, v = _eccentric_sample()
    grid = build_accumulation_grid(t, r, v, 300, GridConfig())
    dec = build_decision_grid(grid, GridConfig())

    assert dec.interval_of.shape == grid.widths.shape
    assert dec.interval_of.min() >= 0
    assert dec.interval_of.max() == len(dec) - 1
    recovered = np.concatenate([dec.cells_in(q) for q in range(len(dec))])
    assert np.array_equal(np.sort(recovered), np.arange(len(grid)))


def test_intervals_are_contiguous_in_time() -> None:
    """``g(i)`` must be non-decreasing; a shuffled map would still sum right."""
    t, r, v = _eccentric_sample()
    grid = build_accumulation_grid(t, r, v, 300, GridConfig())
    dec = build_decision_grid(grid, GridConfig())
    assert np.all(np.diff(dec.interval_of) >= 0)


def test_decision_grid_rejects_an_interval_with_no_cell() -> None:
    """A degree chosen for a stretch that contributes nothing is a bug."""
    t, r, v = _eccentric_sample()
    coarse = GridConfig(samples_per_tau=2.0, dt_acc_min_s=400.0,
                        dt_acc_max_s=600.0, dt_dec_s=10.0)
    grid = build_accumulation_grid(t, r, v, 300, coarse)
    with pytest.raises(ValueError, match="no accumulation cell"):
        build_decision_grid(grid, coarse)


@pytest.mark.parametrize(("dt_dec_s", "message"), [
    (0.0, "dt_dec_s must be positive"),
    (2.0 * DURATION, "exceeds the arc length"),
])
def test_decision_grid_rejects_degenerate_settings(dt_dec_s, message) -> None:
    t, r, v = _eccentric_sample()
    grid = build_accumulation_grid(t, r, v, 300, GridConfig())
    with pytest.raises(ValueError, match=message):
        build_decision_grid(grid, GridConfig(dt_dec_s=dt_dec_s))
