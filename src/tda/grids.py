"""Accumulation and decision grids.

Two grids, never conflated (manuscript §3.5).  The *accumulation* grid carries
the integrand and is refined where the omitted acceleration decorrelates
quickly; the *decision* grid carries the piecewise-constant degree and is
coarse enough to be flown.  This module builds both and nothing else: it takes
a sampled trajectory and returns index structures, so it is pure, fast and
testable without a gravity model.

Which quadrature rule the prefix structure allows
-------------------------------------------------
The objective's tractability rests on
:math:`\\mathbf S_j=\\sum_{i\\le j}\\mathbf u_i` being a plain prefix sum,
which requires the inner quadrature weight of epoch :math:`i` to be
independent of the outer epoch :math:`j`.  That rules out the trapezoid rule
for the inner integral: its endpoint coefficient is
:math:`(t_j-t_{j-1})/2`, which depends on :math:`j`, and carrying it would
turn :math:`\\mathbf Q_{ik}=\\mathbf A_{\\max(i,k)}` into that plus a
:math:`j`-dependent correction.

A rectangle rule keeps the structure but is only first-order, which at the
campaign's default of two samples per correlation time is not a detail.  The
resolution used here costs nothing and keeps both properties: the
accumulation **nodes are cell midpoints** and the objective is sampled at
**cell edges**.  The inner weight is then the full cell width -- independent
of :math:`j`, so the prefix structure survives verbatim -- while the rule
itself is the midpoint rule and second-order.  The prefix sum through cell
:math:`j` is exactly the inner integral up to edge :math:`j{+}1`, so the two
grids interlock rather than approximate one another.

The outer integral has no such constraint (there is no prefix sum over
:math:`j`), so it uses the trapezoid rule on the edges.  Both weight sets are
non-negative, which is what makes :math:`\\mathbf Q` positive semidefinite and
the Frank--Wolfe relaxation convex; higher-order Newton--Cotes rules would not
be, and are not offered.

Refinement
----------
Node placement is by equidistribution [deBoor1973]_, the standard device of
adaptive mesh generation [Huang2011]_: a node density
:math:`\\rho(t)=n_s/\\tau_{\\mathrm{corr}}(t)` is integrated to a monotone
:math:`S(t)`, and the edges are placed where :math:`S` crosses the integers.
Clamping the *density* rather than the resulting steps is what makes the step
bounds hold cell by cell: each cell spans unit density, so its width is the
reciprocal of its own mean density and therefore inherits the clamp.  A
post-hoc clip of the widths would move the edges and leave the neighbours
outside the bound it was applied for.

One residual: the node count must be an integer, so the density is rescaled by
:math:`\\lceil S(T)\\rceil/S(T)` to land the last edge exactly on :math:`T`.
That tightens every cell by at most one part in :math:`M`, so the lower bound
can be undershot by that much and no more.  Snapping the last edge instead
would leave a single cell up to a full step outside the bound, which is the
worse trade.

References
----------
.. [deBoor1973] C. de Boor, "Good approximation by splines with variable
   knots II", in *Conference on the Numerical Solution of Differential
   Equations*, Lecture Notes in Mathematics 363, Springer, 1973 -- the
   equidistribution principle.
.. [Huang2011] W. Huang and R. D. Russell, *Adaptive Moving Mesh Methods*,
   Springer, 2011, Ch. 2 -- equidistribution in one dimension.
.. [Davis1984] P. J. Davis and P. Rabinowitz, *Methods of Numerical
   Integration*, 2nd ed., Academic Press, 1984, §2.1 -- midpoint and
   trapezoid rules and their orders.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from tda.config import GridConfig

__all__ = [
    "AccumulationGrid",
    "DecisionGrid",
    "build_accumulation_grid",
    "build_decision_grid",
    "correlation_time",
    "trapezoid_weights",
]

Arr = NDArray[np.float64]
IntArr = NDArray[np.int_]


# ---------------------------------------------------------------------------
# Correlation time
# ---------------------------------------------------------------------------


def correlation_time(radius: Arr, speed: Arr, degree: int) -> Arr:
    """Along-track decorrelation time of the omitted acceleration.

    .. math::
        \\tau_{\\mathrm{corr}} \\;=\\; \\frac{\\pi r}{N\\,v} ,

    the time to traverse a half-wavelength of the leading omitted harmonic at
    the spacecraft's own radius.  This is the scale on which the *direction*
    of the defect changes, and therefore the scale a signed integral has to
    resolve -- a magnitude statistic can be converged on a much coarser grid,
    which is why the previous campaign's convergence check does not carry
    over.

    Which degree
    ------------
    ``degree`` should be the **largest** degree the campaign will evaluate on
    this arc, not the policy degree.  The grid is shared by every candidate,
    and the finest texture belongs to the highest truncation; refining on a
    lower degree would alias the very comparisons the campaign is built to
    make.  Passing the reference degree is the conservative choice.

    Parameters
    ----------
    radius:
        Selenocentric radius, m, shape ``(K,)``.
    speed:
        Inertial speed, m s^-1, shape ``(K,)``.
    degree:
        Truncation degree :math:`N`, positive.

    Returns
    -------
    ndarray, shape (K,)
        Correlation time in seconds.

    Raises
    ------
    ValueError
        If ``degree`` is not positive, or any radius or speed is not positive
        -- either would make the expression meaningless rather than large.

    Examples
    --------
    A 100 km circular lunar orbit at degree 300 decorrelates in a few seconds.

    >>> import numpy as np
    >>> tau = correlation_time(np.array([1.8374e6]), np.array([1633.0]), 300)
    >>> bool(1.0 < tau[0] < 20.0)
    True
    """
    if degree <= 0:
        raise ValueError(f"degree must be positive, got {degree}")
    radius = np.asarray(radius, dtype=float)
    speed = np.asarray(speed, dtype=float)
    if radius.shape != speed.shape:
        raise ValueError(
            f"radius {radius.shape} and speed {speed.shape} must match")
    if np.any(radius <= 0.0) or np.any(speed <= 0.0):
        raise ValueError("radius and speed must be strictly positive")
    return math.pi * radius / (degree * speed)


# ---------------------------------------------------------------------------
# Quadrature weights
# ---------------------------------------------------------------------------


def trapezoid_weights(nodes: Arr) -> Arr:
    """Non-negative trapezoid weights for a possibly non-uniform grid.

    For nodes :math:`x_0<\\dots<x_K` the composite trapezoid rule is
    :math:`\\sum_k w_k f(x_k)` with
    :math:`w_0=(x_1-x_0)/2`, :math:`w_K=(x_K-x_{K-1})/2` and
    :math:`w_k=(x_{k+1}-x_{k-1})/2` in between [Davis1984]_.  The weights sum
    to :math:`x_K-x_0` and are non-negative for any increasing grid, which is
    the property :math:`\\mathbf Q\\succeq0` depends on.

    Parameters
    ----------
    nodes:
        Strictly increasing, shape ``(K+1,)`` with ``K >= 1``.

    Returns
    -------
    ndarray, shape (K+1,)

    Raises
    ------
    ValueError
        If the grid has fewer than two nodes or is not strictly increasing.

    Examples
    --------
    >>> import numpy as np
    >>> trapezoid_weights(np.array([0.0, 1.0, 3.0]))
    array([0.5, 1.5, 1. ])
    """
    nodes = np.asarray(nodes, dtype=float)
    if nodes.ndim != 1 or nodes.size < 2:
        raise ValueError("need at least two nodes")
    steps = np.diff(nodes)
    if np.any(steps <= 0.0):
        raise ValueError("nodes must be strictly increasing")
    weights = np.empty_like(nodes)
    weights[0] = 0.5 * steps[0]
    weights[-1] = 0.5 * steps[-1]
    weights[1:-1] = 0.5 * (steps[:-1] + steps[1:])
    return weights


# ---------------------------------------------------------------------------
# Accumulation grid
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class AccumulationGrid:
    """Cells of the accumulation grid, plus the epochs each rule samples.

    Attributes
    ----------
    edges:
        Cell boundaries, shape ``(M+1,)``, ``edges[0] == 0`` and
        ``edges[-1] == T``.  The objective is sampled here.
    nodes:
        Cell midpoints, shape ``(M,)``.  The defect and the transported
        contributions :math:`\\mathbf u_i` are evaluated here.
    widths:
        Cell widths :math:`\\Delta t_i`, shape ``(M,)``.  These are the inner
        quadrature weights; they sum to ``T`` and are independent of the outer
        epoch, which is what preserves the prefix structure.
    outer_weights:
        :math:`\\omega_j\\ge0` on ``edges``, shape ``(M+1,)``, summing to
        ``T``.

    Notes
    -----
    Every epoch a propagation must visit is
    :meth:`required_epochs`; passing anything less would force the state at
    the missing points to be interpolated, and the campaign does not
    interpolate a trajectory it can integrate.
    """

    edges: Arr
    nodes: Arr
    widths: Arr
    outer_weights: Arr

    def __len__(self) -> int:
        """Number of cells, :math:`M`."""
        return int(self.widths.size)

    @property
    def duration(self) -> float:
        """Arc length :math:`T`, seconds."""
        return float(self.edges[-1])

    def required_epochs(self) -> Arr:
        """Sorted union of edges and midpoints, for use as ``t_eval``.

        Returns
        -------
        ndarray, shape (2M+1,)
            Strictly increasing, starting at zero.
        """
        merged = np.empty(self.edges.size + self.nodes.size, dtype=float)
        merged[0::2] = self.edges
        merged[1::2] = self.nodes
        return merged

    def edge_indices(self) -> IntArr:
        """Positions of the edges within :meth:`required_epochs`."""
        return np.arange(0, 2 * len(self) + 1, 2)

    def node_indices(self) -> IntArr:
        """Positions of the midpoints within :meth:`required_epochs`."""
        return np.arange(1, 2 * len(self) + 1, 2)


def build_accumulation_grid(sample_t: Arr, sample_radius: Arr,
                            sample_speed: Arr, degree: int,
                            cfg: GridConfig) -> AccumulationGrid:
    """Refine an accumulation grid on the correlation time.

    The target step is :math:`\\Delta t=\\tau_{\\mathrm{corr}}/n_s`, which
    varies by more than an order of magnitude along an eccentric lunar orbit:
    a uniform grid fine enough for perilune spends most of its samples where
    nothing happens, and one coarse enough for apolune aliases the passage
    that carries the signal.

    Method
    ------
    Equidistribution [deBoor1973]_.  The node density
    :math:`\\rho=n_s/\\tau_{\\mathrm{corr}}` is clamped to
    :math:`[1/\\Delta t_{\\max},\\,1/\\Delta t_{\\min}]`, integrated to a
    strictly increasing :math:`S(t)`, and the edges are placed where
    :math:`S` crosses the integers.  Clamping the density rather than the
    resulting steps is what makes the bounds exact: a post-hoc clip of the
    steps would move the nodes and break the equidistribution it was applied
    to.

    The number of cells is :math:`M=\\lceil S(T)\\rceil`, and the density is
    then rescaled by :math:`M/S(T)` so that the last edge lands exactly on
    :math:`T` rather than being snapped there.  That rescaling tightens every
    cell by at most one part in :math:`M`; the step bounds therefore hold to
    that tolerance rather than exactly, and the module docstring says why the
    alternative is worse.

    Parameters
    ----------
    sample_t:
        Epochs at which the trajectory has been sampled, shape ``(K,)``,
        strictly increasing and starting at zero.  These come from a cheap
        first-pass propagation: the grid must follow the *flown* radius and
        speed, and a Keplerian approximation drifts in phase by minutes over a
        seven-day arc, which would place the refinement beside the perilune it
        was meant for rather than on it.
    sample_radius, sample_speed:
        Radius and speed at those epochs, SI, shape ``(K,)``.
    degree:
        Degree to refine on; see :func:`correlation_time`.
    cfg:
        Grid settings.

    Returns
    -------
    AccumulationGrid

    Raises
    ------
    ValueError
        If the sample grid is malformed, or the settings are degenerate
        (non-positive ``samples_per_tau``, or ``dt_acc_min_s`` not below
        ``dt_acc_max_s``).

    Notes
    -----
    Cost.  ``M`` is decided by the physics, not by the caller, and on an
    eccentric arc it can be large.  The caller is expected to check
    ``len(grid)`` against the memory budget of the defect table rather than
    discovering it downstream.
    """
    sample_t = np.asarray(sample_t, dtype=float)
    if sample_t.ndim != 1 or sample_t.size < 2:
        raise ValueError("sample_t must hold at least two epochs")
    if sample_t[0] != 0.0:
        raise ValueError(f"sample_t must start at 0.0, got {sample_t[0]}")
    if np.any(np.diff(sample_t) <= 0.0):
        raise ValueError("sample_t must be strictly increasing")
    if cfg.samples_per_tau <= 0.0:
        raise ValueError(
            f"samples_per_tau must be positive, got {cfg.samples_per_tau}")
    if not 0.0 < cfg.dt_acc_min_s < cfg.dt_acc_max_s:
        raise ValueError(
            f"need 0 < dt_acc_min_s < dt_acc_max_s, got "
            f"{cfg.dt_acc_min_s} and {cfg.dt_acc_max_s}")

    tau = correlation_time(sample_radius, sample_speed, degree)
    density = np.clip(cfg.samples_per_tau / tau,
                      1.0 / cfg.dt_acc_max_s, 1.0 / cfg.dt_acc_min_s)

    # Monotone node count S(t); strictly increasing because density > 0.
    cumulative = np.concatenate(
        [[0.0], np.cumsum(0.5 * (density[:-1] + density[1:]) * np.diff(sample_t))]
    )
    n_cells = max(1, math.ceil(cumulative[-1]))
    # Rescale so the last edge is exactly T instead of being snapped to it.
    cumulative *= n_cells / cumulative[-1]

    edges = np.interp(np.arange(n_cells + 1, dtype=float),
                      cumulative, sample_t)
    edges[0] = 0.0
    edges[-1] = float(sample_t[-1])

    widths = np.diff(edges)
    if np.any(widths <= 0.0):
        raise ValueError(
            "equidistribution produced a non-increasing edge sequence; "
            "the sampled trajectory is probably too coarse for the requested "
            f"samples_per_tau={cfg.samples_per_tau}")

    return AccumulationGrid(
        edges=edges,
        nodes=0.5 * (edges[:-1] + edges[1:]),
        widths=widths,
        outer_weights=trapezoid_weights(edges),
    )


# ---------------------------------------------------------------------------
# Decision grid
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class DecisionGrid:
    """Intervals on which the degree is held constant.

    Attributes
    ----------
    edges:
        Interval boundaries, shape ``(K+1,)``.
    interval_of:
        :math:`g(i)`, the interval each accumulation cell belongs to, shape
        ``(M,)``.
    time_weight:
        :math:`W_q=\\sum_{i\\in I_q}\\Delta t_i`, shape ``(K,)``, the weight
        the budget constraint uses.

    Notes
    -----
    ``time_weight`` is a sum of cell *widths*, not a count of cells.  On a
    refined grid the two differ by a large factor, and using the count would
    charge the budget as though perilune occupied most of the arc merely
    because it is sampled most densely.
    """

    edges: Arr
    interval_of: IntArr
    time_weight: Arr

    def __len__(self) -> int:
        """Number of decision intervals, :math:`K_{\\mathrm{dec}}`."""
        return int(self.time_weight.size)

    def cells_in(self, interval: int) -> IntArr:
        """Indices of the accumulation cells inside one decision interval.

        Returns
        -------
        ndarray
            The fibre :math:`I_q`.  Contiguous by construction, but returned
            as explicit indices so that callers need not rely on that.
        """
        return np.flatnonzero(self.interval_of == interval)


def build_decision_grid(grid: AccumulationGrid,
                        cfg: GridConfig) -> DecisionGrid:
    """Map accumulation cells onto uniform decision intervals.

    The decision grid is uniform even though the accumulation grid is not.
    That is deliberate and is the reason the comparison is about information
    rather than resolution: every policy in the campaign -- the benchmark that
    sees the whole arc and the controller that sees one interval -- switches
    on the *same* boundaries, so a margin cannot come from one of them being
    allowed to switch more often.

    A cell is assigned to an interval by its **midpoint**, so each cell
    belongs to exactly one interval and the widths partition the arc exactly.
    Assigning by overlap would split cells across boundaries and make
    :math:`W_q` depend on how the split is apportioned.

    Parameters
    ----------
    grid:
        The accumulation grid to partition.
    cfg:
        Grid settings; ``dt_dec_s`` is used.

    Returns
    -------
    DecisionGrid

    Raises
    ------
    ValueError
        If ``dt_dec_s`` is not positive, exceeds the arc, or is finer than the
        coarsest accumulation cell -- the last of which would leave decision
        intervals containing no cell at all, so that a degree would be chosen
        for a stretch of arc that contributes nothing to the objective.
    """
    if cfg.dt_dec_s <= 0.0:
        raise ValueError(f"dt_dec_s must be positive, got {cfg.dt_dec_s}")
    duration = grid.duration
    if cfg.dt_dec_s > duration:
        raise ValueError(
            f"dt_dec_s={cfg.dt_dec_s} exceeds the arc length {duration}")

    n_intervals = math.ceil(duration / cfg.dt_dec_s)
    edges = np.linspace(0.0, duration, n_intervals + 1)

    interval_of = np.clip(
        np.searchsorted(edges, grid.nodes, side="right") - 1,
        0, n_intervals - 1,
    ).astype(np.int64)

    time_weight = np.bincount(interval_of, weights=grid.widths,
                              minlength=n_intervals)
    empty = np.flatnonzero(time_weight <= 0.0)
    if empty.size:
        raise ValueError(
            f"{empty.size} decision interval(s) contain no accumulation cell "
            f"(first is index {int(empty[0])}); dt_dec_s={cfg.dt_dec_s} is "
            f"finer than the coarsest cell {float(grid.widths.max()):.3g} s")

    return DecisionGrid(edges=edges, interval_of=interval_of,
                        time_weight=time_weight)
