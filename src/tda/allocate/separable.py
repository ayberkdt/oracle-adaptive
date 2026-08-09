"""The two separable rungs: ``A-force`` and ``A-sens``.

Both minimise a sum of per-interval terms under one multiplier, so both
decouple and a single bisection solves them exactly.  They are the ladder's
lower rungs and the denominator of RQ2: whatever ``A-sign`` gains beyond
``A-sens`` is what the *signed cross-epoch* coupling buys, as opposed to what
correcting the weighting buys.

.. math::
    \\textsc{a-force}:\\;\\; \\argmin_{N}\\Bigl[\\textstyle\\sum_{i\\in I_q}
      \\lVert\\Delta\\mathbf a_i(N)\\rVert^2\\Delta t_i+\\lambda W_qN^2\\Bigr],
    \\qquad
    \\textsc{a-sens}:\\;\\; \\argmin_{N}\\Bigl[\\textstyle\\sum_{i\\in I_q}
      \\Delta t_i\\,\\Delta\\mathbf a_i(N)^{\\top}\\mathbf K_i\\,
      \\Delta\\mathbf a_i(N)+\\lambda W_qN^2\\Bigr].

The weight is not borrowed
--------------------------
``A-sens`` uses the objective's own kernel evaluated on its diagonal,
:math:`\\mathbf K_i`, and not a scalar such as
:math:`\\lVert\\Phi(T,t)\\mathbf B\\rVert^2` lifted from a terminal proxy.
Nor does it use the raw discrete diagonal :math:`\\mathbf u_i^{\\top}
\\mathbf Q_{ii}\\mathbf u_i`, which carries :math:`\\Delta t_i^2` and shrinks
as the grid is refined --- on an adaptive grid that would also manufacture a
within-arc trend, since the cell width varies along the arc
(``DECISIONS.md`` D110).

Why the bisection is exact here, and where it still stops short
---------------------------------------------------------------
For a fixed multiplier each interval is minimised independently over a finite
candidate set, so the solve is a global minimum, and for a global minimiser
the attained work is non-increasing in the multiplier.  The bisection of
:func:`tda.allocate.budget.solve_to_budget` is therefore sound without the
monotonicity check that the coordinate descent needs.

What it does not deliver is a schedule that spends the ceiling.  Over a
discrete candidate set the attained work is a *step* function of the
multiplier, so between the largest feasible multiplier and the smallest
infeasible one the work jumps, and the bisection returns the point below the
jump.  Measured on a pilot arc that left ``A-force`` at forty-five per cent of
its ceiling and, worse, returned the *same* schedule at two ceilings a factor
of 1.44 apart --- a comparator that does not respond to its budget is not a
comparator (``DECISIONS.md`` D177).  This is exactly the duality gap
[Shoham1988]_ describes: a multiplier sweep recovers only the supported
efficient solutions, and the unsupported ones are unreachable from any
multiplier.

:func:`greedy_fill` closes it by marginal analysis, the classical completion
of an Everett sweep: repeatedly take the single interval change that buys the
most criterion per unit of extra work and still fits.  It is not a
"spend the ceiling" rule.  It only accepts moves that *reduce* the criterion,
which matters because the criterion is **not** monotone in degree ---
:math:`\\Delta\\mathbf a(N)=-\\sum_{n>N}\\mathbf a_n` is the norm of a partial
sum, and dropping a band can leave a larger residual when the bands cancel.
A rule that raised degrees until the budget was gone would sometimes make the
comparator worse, and would do it silently.

References
----------
.. [Shoham1988] Y. Shoham and A. Gersho, "Efficient bit allocation for an
   arbitrary set of quantizers", *IEEE Trans. ASSP* 36(9), 1988.
.. [Everett1963] H. Everett III, "Generalized Lagrange multiplier method for
   solving problems of optimum allocation of resources", *Operations
   Research* 11(3), 1963 -- and the marginal-analysis completion of the sweep.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from tda.allocate.budget import ScheduleSolution, solve_to_budget

__all__ = ["force_values", "greedy_fill", "sensitivity_values",
           "solve_separable"]

Arr = NDArray[np.float64]
IntArr = NDArray[np.int_]


def _accumulate(per_cell: Arr, interval_of: IntArr, n_intervals: int) -> Arr:
    """Sum a ``(M, P)`` per-cell array into ``(K, P)`` per-interval."""
    out = np.zeros((n_intervals, per_cell.shape[1]), dtype=float)
    np.add.at(out, interval_of, per_cell)
    return out


def force_values(defect: Arr, widths: Arr, interval_of: IntArr,
                 n_intervals: int) -> Arr:
    """Per-interval value of the force-level criterion, shape ``(K, P)``.

    :math:`\\sum_{i\\in I_q}\\lVert\\Delta\\mathbf a_i(N)\\rVert^2\\Delta t_i`
    --- the previous paper's benchmark, which never sees :math:`\\Phi` at all.

    Parameters
    ----------
    defect:
        :math:`\\Delta\\mathbf a(m_i,N)`, shape ``(M, P, 3)``.
    widths:
        Cell widths, shape ``(M,)``.
    interval_of:
        :math:`g(i)`, shape ``(M,)``.
    n_intervals:
        :math:`K_{\\mathrm{dec}}`.

    Returns
    -------
    ndarray, shape (K, P)
    """
    per_cell = np.einsum("ipa,ipa->ip", defect, defect) * widths[:, None]
    return _accumulate(per_cell, interval_of, n_intervals)


def sensitivity_values(defect: Arr, widths: Arr, local_kernel: Arr,
                       interval_of: IntArr, n_intervals: int) -> Arr:
    """Per-interval value of the sensitivity-weighted criterion, ``(K, P)``.

    :math:`\\sum_{i\\in I_q}\\Delta t_i\\,\\Delta\\mathbf a_i(N)^{\\top}
    \\mathbf K_i\\,\\Delta\\mathbf a_i(N)`, with :math:`\\mathbf K_i` from
    :meth:`tda.kernel.CouplingKernel.local_kernel`.

    Parameters
    ----------
    defect:
        Shape ``(M, P, 3)``.
    widths:
        Shape ``(M,)``.
    local_kernel:
        :math:`\\mathbf K_i`, shape ``(M, 3, 3)``, symmetric.
    interval_of, n_intervals:
        The decision-grid map and its size.

    Returns
    -------
    ndarray, shape (K, P)

    Raises
    ------
    ValueError
        If the shapes disagree.
    """
    if local_kernel.shape != (defect.shape[0], 3, 3):
        raise ValueError(
            f"local_kernel must have shape ({defect.shape[0]}, 3, 3), got "
            f"{local_kernel.shape}")
    per_cell = np.einsum("ipa,iab,ipb->ip", defect, local_kernel, defect)
    return _accumulate(per_cell * widths[:, None], interval_of, n_intervals)


def greedy_fill(values: Arr, cost: Arr, picked: IntArr, budget: float,
                max_moves: int | None = None) -> tuple[IntArr, int]:
    """Improve a schedule by marginal analysis until nothing fits.

    Parameters
    ----------
    values:
        Per-interval, per-candidate criterion, shape ``(K, P)``.
    cost:
        :math:`W_qN^2` for the same grid, shape ``(K, P)``.
    picked:
        Candidate column per interval, shape ``(K,)``; the bisection's answer.
    budget:
        The work ceiling.
    max_moves:
        Cap on accepted moves.  Defaults to ``K * P``.  Each move strictly
        decreases the criterion so the loop cannot cycle, but a cap keeps a
        pathological instance from consuming a campaign.

    Returns
    -------
    tuple of (ndarray, int)
        The improved columns and how many moves were accepted.

    Raises
    ------
    ValueError
        If the shapes disagree or the incoming schedule already overspends ---
        the fill improves a feasible point and cannot rescue an infeasible
        one.

    Notes
    -----
    Two kinds of move are accepted and the order matters.  A move that lowers
    the criterion *and* costs no more is taken unconditionally: it dominates,
    and taking it first releases work for the moves that need it.  Otherwise
    the move maximising the criterion reduction per unit of extra work is
    taken, which is the marginal-analysis rule and the same quantity the
    multiplier prices.

    A move is only ever accepted if it lowers the criterion, so an instance
    where the incoming schedule is already the best reachable point returns
    unchanged with zero moves.  That is the case the campaign must be able to
    distinguish from a fill that was never attempted, which is why the count
    is returned rather than inferred.
    """
    values = np.asarray(values, dtype=float)
    cost = np.asarray(cost, dtype=float)
    if values.shape != cost.shape or values.ndim != 2:
        raise ValueError(
            f"values {values.shape} and cost {cost.shape} must be equal 2-D")
    picked = np.array(picked, dtype=np.int64)
    if picked.shape != (values.shape[0],):
        raise ValueError(
            f"picked must have shape ({values.shape[0]},), got {picked.shape}")

    rows = np.arange(values.shape[0])
    work = float(cost[rows, picked].sum())
    if work > budget * (1.0 + 1.0e-12):
        raise ValueError(
            f"the incoming schedule spends {work:.6g} against a ceiling of "
            f"{budget:.6g}; the fill improves a feasible point")

    limit = values.size if max_moves is None else max_moves
    moves = 0
    while moves < limit:
        gain = values - values[rows, picked][:, None]          # < 0 is better
        extra = cost - cost[rows, picked][:, None]
        room = budget - work

        improves = gain < 0.0
        free = improves & (extra <= 0.0)
        if free.any():
            # Dominating moves first: better and no more expensive.
            score = np.where(free, gain, 0.0)
        else:
            affordable = improves & (extra <= room * (1.0 + 1.0e-12))
            if not affordable.any():
                break
            with np.errstate(divide="ignore", invalid="ignore"):
                rate = np.where(affordable, gain / extra, 0.0)
            score = rate
        flat = int(np.argmin(score))
        if score.ravel()[flat] >= 0.0:
            break
        interval, column = divmod(flat, values.shape[1])
        work += float(cost[interval, column] - cost[interval, picked[interval]])
        picked[interval] = column
        moves += 1
    return picked, moves


def solve_separable(values: Arr, time_weight: Arr,
                    candidate_degrees: tuple[int, ...],
                    budget: float, fill: bool = True) -> ScheduleSolution:
    """Minimise a separable criterion under the work ceiling.

    Parameters
    ----------
    values:
        Per-interval, per-candidate criterion, shape ``(K, P)``, from
        :func:`force_values` or :func:`sensitivity_values`.
    time_weight:
        :math:`W_q`, shape ``(K,)``.
    candidate_degrees:
        The candidate axis, shape ``(P,)``, matching ``values``.
    budget:
        The work ceiling.
    fill:
        Whether to run :func:`greedy_fill` on the bisection's answer.  On by
        default because a comparator that leaves half its ceiling unspent is a
        straw man; ``False`` is the control that shows what the fill was
        worth, and is how the pilot measurement of D177 was made.

    Returns
    -------
    ScheduleSolution
        ``objective`` is the separable criterion's own value, **not** the
        trajectory objective :math:`J`.  The two are different functionals and
        the campaign never quotes one for the other; evaluate the returned
        schedule through :class:`tda.kernel.CouplingKernel` to get :math:`J`.
        ``diagnostics["fill_moves"]`` records how many marginal moves the fill
        accepted; zero means the bisection's answer was already the best
        reachable point, which is a different statement from the fill not
        having run.

    Raises
    ------
    ValueError
        If the shapes disagree.

    Notes
    -----
    Ties are broken toward the **lower** degree.  A tie means two candidates
    buy the same criterion value at different cost, and the cheaper one is the
    one a budgeted allocation should take; ``argmin`` on the penalized array
    already does this because the penalty separates them, but exact ties
    survive at :math:`\\lambda=0` and the convention is fixed here rather than
    left to array order.
    """
    values = np.asarray(values, dtype=float)
    time_weight = np.asarray(time_weight, dtype=float)
    degrees = np.asarray(candidate_degrees, dtype=float)
    if values.ndim != 2 or values.shape[1] != degrees.size:
        raise ValueError(
            f"values {values.shape} must be (K, {degrees.size})")
    if time_weight.shape != (values.shape[0],):
        raise ValueError(
            f"time_weight {time_weight.shape} must be ({values.shape[0]},)")

    cost = time_weight[:, None] * degrees[None, :] ** 2

    def solve_at(multiplier: float) -> tuple[IntArr, float]:
        penalized = values + multiplier * cost
        # Lexicographic argmin: value first, then cost, so a tie takes the
        # cheaper degree rather than whichever index came first.
        picked = np.lexsort((cost, penalized), axis=1)[:, 0]
        chosen = degrees[picked].astype(np.int64)
        return chosen, float(values[np.arange(values.shape[0]), picked].sum())

    solution = solve_to_budget(solve_at, time_weight, budget)
    if not fill:
        return solution

    column = {int(d): p for p, d in enumerate(candidate_degrees)}
    picked = np.array([column[int(n)] for n in solution.degrees],
                      dtype=np.int64)
    picked, moves = greedy_fill(values, cost, picked, budget)
    rows = np.arange(values.shape[0])
    diagnostics = dict(solution.diagnostics)
    diagnostics["fill_moves"] = float(moves)
    return ScheduleSolution(
        degrees=degrees[picked].astype(np.int64),
        objective=float(values[rows, picked].sum()),
        work=float(cost[rows, picked].sum()),
        budget=budget,
        multiplier=solution.multiplier,
        feasible=True,
        spread=solution.spread,
        diagnostics=diagnostics,
    )
