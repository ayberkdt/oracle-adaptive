"""Verification of the accumulation and decision grids.

Three properties carry the campaign and are asserted rather than assumed:
the inner weights sum to the arc length (so the budget constraint means what
it says), the refinement equidistributes the correlation time (so ``n_s`` is
the quantity it is swept as), and **every decision boundary is an accumulation
edge** (so no cell straddles a degree switch).
"""

from __future__ import annotations

import numpy as np
import pytest

from tda.config import GridConfig
from tda.grids import (
    build_accumulation_grid,
    build_decision_edges,
    build_decision_grid,
    correlation_time,
    trapezoid_weights,
)

DURATION = 7200.0


def _eccentric_sample(n: int = 4001):
    """A sampled arc whose radius swings by a factor of three.

    Not a real orbit: an analytic stand-in with the property that matters
    here, namely that the correlation time varies strongly along the arc.
    Using a real propagation would test the propagator, not the grid.
    """
    t = np.linspace(0.0, DURATION, n)
    phase = 2.0 * np.pi * t / DURATION
    radius = 1.8e6 + 1.8e6 * (1.0 - np.cos(phase))     # 1.8e6 .. 5.4e6
    speed = np.sqrt(4.9028e12 * (2.0 / radius - 1.0 / 3.6e6))
    return t, radius, speed


def _build(cfg: GridConfig, degree: int = 300):
    """Build both grids against one set of decision edges."""
    t, r, v = _eccentric_sample()
    edges = build_decision_edges(DURATION, cfg)
    grid = build_accumulation_grid(t, r, v, degree, cfg, edges)
    return grid, build_decision_grid(grid, edges), edges


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
# Decision edges
# ---------------------------------------------------------------------------


def test_decision_edges_are_uniform_and_span_the_arc() -> None:
    edges = build_decision_edges(DURATION, GridConfig(dt_dec_s=120.0))
    assert edges[0] == 0.0
    assert edges[-1] == DURATION
    assert np.allclose(np.diff(edges), 120.0, rtol=0.0, atol=1e-9)


@pytest.mark.parametrize(("duration", "cfg", "message"), [
    (DURATION, GridConfig(dt_dec_s=0.0), "dt_dec_s must be positive"),
    (DURATION, GridConfig(dt_dec_s=2 * DURATION), "exceeds the arc length"),
    (0.0, GridConfig(), "duration must be positive"),
    (DURATION, GridConfig(dt_dec_s=1.0, dt_acc_min_s=2.0),
     "below dt_acc_min_s"),
])
def test_decision_edges_reject_degenerate_settings(duration, cfg,
                                                   message) -> None:
    with pytest.raises(ValueError, match=message):
        build_decision_edges(duration, cfg)


# ---------------------------------------------------------------------------
# Accumulation grid
# ---------------------------------------------------------------------------


def test_grid_spans_the_arc_exactly() -> None:
    """The inner weights are the budget's denominator; they must sum to T."""
    grid, _, _ = _build(GridConfig())

    assert grid.edges[0] == 0.0
    assert grid.edges[-1] == DURATION
    assert grid.widths.sum() == pytest.approx(DURATION, rel=1e-12)
    assert grid.outer_weights.sum() == pytest.approx(DURATION, rel=1e-12)
    assert np.all(grid.outer_weights >= 0.0)


def test_every_decision_boundary_is_an_accumulation_edge() -> None:
    """The property that keeps a cell from straddling a degree switch.

    Without it a cell spanning a boundary is charged whole to one degree while
    the policy switches inside it, which puts a boundary error of up to one
    cell width into a signed accumulation.
    """
    grid, _, decision_edges = _build(GridConfig())
    matched = np.isin(np.round(decision_edges, 9), np.round(grid.edges, 9))
    assert matched.all(), f"{(~matched).sum()} boundaries are not edges"


def test_no_cell_crosses_a_decision_boundary() -> None:
    """The same property stated on the cells rather than on the edges."""
    grid, _, decision_edges = _build(GridConfig())
    interior = decision_edges[1:-1]
    for lo, hi in zip(grid.edges[:-1], grid.edges[1:], strict=True):
        crossing = np.flatnonzero((interior > lo + 1e-9)
                                  & (interior < hi - 1e-9))
        assert crossing.size == 0, f"cell [{lo}, {hi}] crosses a boundary"


def test_nodes_are_the_cell_midpoints() -> None:
    """The midpoint rule is what makes the inner quadrature second-order."""
    grid, _, _ = _build(GridConfig())
    assert np.allclose(grid.nodes, 0.5 * (grid.edges[:-1] + grid.edges[1:]),
                       rtol=0.0, atol=1e-9)
    assert np.all(grid.nodes > grid.edges[:-1])
    assert np.all(grid.nodes < grid.edges[1:])


def test_step_bounds_follow_from_the_density_clamp() -> None:
    """Clipping the cell count carries the bounds into the widths -- inexactly.

    An interval demanding ``S`` cells is given ``round(S)``, so the widths are
    scaled by ``S/m`` and can miss a bound by ``1/(2m)``.  The assertion is
    that bound and not a tighter one: at twenty-four cells the undershoot is
    two per cent, and it was a real test failure that established the number
    rather than an allowance made in advance.
    """
    cfg = GridConfig(samples_per_tau=2.0, dt_acc_min_s=5.0, dt_acc_max_s=60.0)
    grid, _, _ = _build(cfg)

    # The exact bound: a cell carrying density mass alpha over a region where
    # the clamped density lies in [1/dt_max, 1/dt_min] has its width between
    # alpha*dt_min and alpha*dt_max.  No rounding argument enters.
    alpha = grid.density_mass
    assert np.all(grid.widths >= alpha * cfg.dt_acc_min_s - 1e-9)
    assert np.all(grid.widths <= alpha * cfg.dt_acc_max_s + 1e-9)

    # The clamp must actually bind somewhere, or the test proves nothing.
    assert grid.widths.min() < 1.5 * cfg.dt_acc_min_s
    # And the *requested* bound is missed only by the integer rounding, which
    # is 1/(2m) per interval -- a couple of per cent, not a factor.
    assert grid.widths.min() > 0.95 * cfg.dt_acc_min_s


def test_refinement_follows_the_correlation_time() -> None:
    """Cells must be short where the correlation time is short.

    Measured as a rank correlation between cell width and the local
    correlation time; a uniform grid scores zero and would pass a test that
    only checked the total.
    """
    cfg = GridConfig(samples_per_tau=2.0, dt_acc_min_s=0.5,
                     dt_acc_max_s=600.0)
    grid, _, _ = _build(cfg)

    t, r, v = _eccentric_sample()
    tau_at_nodes = np.interp(grid.nodes, t, correlation_time(r, v, 300))
    order_tau = np.argsort(np.argsort(tau_at_nodes))
    order_width = np.argsort(np.argsort(grid.widths))
    assert np.corrcoef(order_tau, order_width)[0, 1] > 0.99


def test_halving_the_target_step_roughly_doubles_the_cells() -> None:
    coarse, _, _ = _build(GridConfig(samples_per_tau=2.0, dt_acc_min_s=0.5,
                                     dt_acc_max_s=600.0))
    fine, _, _ = _build(GridConfig(samples_per_tau=4.0, dt_acc_min_s=0.5,
                                   dt_acc_max_s=600.0))
    assert len(fine) == pytest.approx(2 * len(coarse), rel=0.05)


def test_required_epochs_interleave_edges_and_nodes() -> None:
    """``t_eval`` must contain every epoch a rule samples, in order."""
    grid, _, _ = _build(GridConfig())
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
    edges = build_decision_edges(DURATION, GridConfig())
    with pytest.raises(ValueError, match=message):
        build_accumulation_grid(t, r, v, 300, GridConfig(**kwargs), edges)


@pytest.mark.parametrize(("sample_t", "message"), [
    (np.array([0.0]), "at least two"),
    (np.array([1.0, 2.0]), "start at 0.0"),
    (np.array([0.0, 2.0, 1.0]), "strictly increasing"),
    (np.array([0.0, 0.5 * DURATION]), "does not cover the arc"),
])
def test_accumulation_grid_rejects_a_bad_sample(sample_t, message) -> None:
    r = np.full(sample_t.size, 1.8e6)
    v = np.full(sample_t.size, 1600.0)
    edges = build_decision_edges(DURATION, GridConfig())
    with pytest.raises(ValueError, match=message):
        build_accumulation_grid(sample_t, r, v, 300, GridConfig(), edges)


# ---------------------------------------------------------------------------
# Decision grid
# ---------------------------------------------------------------------------


def test_decision_weight_is_exactly_the_interval_length() -> None:
    """Boundary alignment makes ``W_q`` exact rather than approximate.

    Before the accumulation grid was refined inside intervals, an interval's
    weight was its length give or take one cell -- up to 35 per cent at a
    120 s interval with 42 s cells.
    """
    cfg = GridConfig(samples_per_tau=2.0, dt_acc_min_s=0.5,
                     dt_acc_max_s=600.0, dt_dec_s=120.0)
    grid, dec, _ = _build(cfg)

    assert np.allclose(dec.time_weight, cfg.dt_dec_s, rtol=0.0, atol=1e-9)
    assert dec.time_weight.sum() == pytest.approx(DURATION, rel=1e-12)
    # The cell *count* still varies several-fold; only the weight is uniform.
    counts = np.bincount(dec.interval_of, minlength=len(dec))
    assert counts.max() > 3 * counts.min()


def test_every_cell_belongs_to_exactly_one_interval() -> None:
    grid, dec, _ = _build(GridConfig())

    assert dec.interval_of.shape == grid.widths.shape
    assert dec.interval_of.min() == 0
    assert dec.interval_of.max() == len(dec) - 1
    recovered = np.concatenate([dec.cells_in(q) for q in range(len(dec))])
    assert np.array_equal(np.sort(recovered), np.arange(len(grid)))


def test_intervals_are_contiguous_in_time() -> None:
    """``g(i)`` must be non-decreasing; a shuffled map would still sum right."""
    _, dec, _ = _build(GridConfig())
    assert np.all(np.diff(dec.interval_of) >= 0)


@pytest.mark.parametrize("other_dt_dec", [60.0, 90.0, 100.0, 250.0])
def test_misaligned_decision_edges_are_detected(other_dt_dec) -> None:
    """Boundaries that are not accumulation edges must not pair silently.

    The empty-interval check alone is not enough.  It catches edges *finer*
    than the cells, but a grid of small cells populates every interval of a
    coarser misaligned grid, so that pairing used to succeed and produce a
    ``W_q`` that was wrong by up to a cell at each end.  The invariant the
    construction guarantees is that every boundary is an edge, and that is
    what is now checked.
    """
    grid, _, _ = _build(GridConfig(dt_dec_s=120.0))
    other = build_decision_edges(DURATION, GridConfig(dt_dec_s=other_dt_dec))
    with pytest.raises(ValueError, match="not accumulation"):
        build_decision_grid(grid, other)


@pytest.mark.parametrize("coarser_dt_dec", [240.0, 360.0, 1200.0])
def test_a_coarser_aligned_grid_is_accepted(coarser_dt_dec) -> None:
    """Subsets of the boundaries are legitimate, and the guard must allow them.

    A 240 s grid built on top of 120 s boundaries still has every boundary on
    an accumulation edge, so no cell straddles a switch and ``W_q`` stays
    exact.  Rejecting it would make the check a coincidence rather than the
    invariant; the ablation that coarsens the decision grid depends on this.
    """
    grid, _, _ = _build(GridConfig(dt_dec_s=120.0))
    coarser = build_decision_edges(DURATION,
                                   GridConfig(dt_dec_s=coarser_dt_dec))
    dec = build_decision_grid(grid, coarser)
    assert np.allclose(dec.time_weight, coarser_dt_dec, rtol=0.0, atol=1e-9)


def test_matched_grids_pair_without_complaint() -> None:
    """The guard must not fire on the pairing the construction produces."""
    grid, dec, edges = _build(GridConfig())
    assert np.allclose(dec.edges, edges, rtol=0.0, atol=0.0)


def test_density_mass_is_uniform_within_an_interval() -> None:
    """Equidistribution splits an interval's demand into equal shares.

    The property the exact step bound rests on: if the shares were unequal,
    ``alpha_i`` would not characterise the cell and the bound would have to
    fall back on the rounding argument.
    """
    grid, dec, _ = _build(GridConfig())
    for q in range(len(dec)):
        share = grid.density_mass[dec.cells_in(q)]
        assert np.allclose(share, share[0], rtol=1e-12, atol=0.0), q


def test_density_mass_totals_the_integrated_density() -> None:
    """The shares must add up to what the density actually demanded."""
    cfg = GridConfig(samples_per_tau=2.0, dt_acc_min_s=0.5, dt_acc_max_s=600.0)
    grid, _, _ = _build(cfg)

    t, r, v = _eccentric_sample()
    density = np.clip(cfg.samples_per_tau / correlation_time(r, v, 300),
                      1.0 / cfg.dt_acc_max_s, 1.0 / cfg.dt_acc_min_s)
    total = float(np.trapezoid(density, t))
    # Not a rounding tolerance: an interval's shares are ``demand/count``
    # repeated ``count`` times, so they sum to the demand exactly and the
    # totals agree to interpolation error.  A loose bound here would let a
    # genuinely wrong split pass.
    assert grid.density_mass.sum() == pytest.approx(total, rel=1e-9)
