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

Why the bisection is exact here
-------------------------------
For a fixed multiplier each interval is minimised independently over a finite
candidate set, so the solve is a global minimum, and for a global minimiser
the attained work is non-increasing in the multiplier.  The bisection of
:func:`tda.allocate.budget.solve_to_budget` is therefore sound without the
monotonicity check that the coordinate descent needs.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from tda.allocate.budget import ScheduleSolution, solve_to_budget

__all__ = ["force_values", "sensitivity_values", "solve_separable"]

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


def solve_separable(values: Arr, time_weight: Arr,
                    candidate_degrees: tuple[int, ...],
                    budget: float) -> ScheduleSolution:
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

    Returns
    -------
    ScheduleSolution
        ``objective`` is the separable criterion's own value, **not** the
        trajectory objective :math:`J`.  The two are different functionals and
        the campaign never quotes one for the other; evaluate the returned
        schedule through :class:`tda.kernel.CouplingKernel` to get :math:`J`.

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

    return solve_to_budget(solve_at, time_weight, budget)
