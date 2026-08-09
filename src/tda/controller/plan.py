"""The offline half: forward--backward allocation on a pilot arc.

The cross term :math:`\\mathbf c_i` depends on the degrees chosen at *every
other* epoch, so a single backward adjoint sweep cannot supply it --- there is
no schedule yet to sweep against.  Algorithm 2 of the manuscript therefore
iterates: accumulate the forward prefix from the current schedule, form the
suffix blocks and the cross terms, re-select every interval, bisect the
multiplier to the budget, repeat.

That iteration is *exactly* the budgeted coordinate descent of
:mod:`tda.allocate.descent`
--------------------------------------------------------------------------
Sweeping in temporal order already carries the prefix forward and reads the
suffix to the right of the block being updated; the forward and backward passes
of Algorithm 2 are the two halves of one sweep, and :math:`J` is the sweep
count.  This module therefore does not re-implement the solver, and that is a
design decision rather than an economy: the benchmark and the plan must differ
in the *information* they are given and in nothing else, or the capture
fraction of RQ3 measures the two solvers against each other instead of
measuring what the reference field was worth (``DECISIONS.md`` D169).

What differs is the contribution tensor.  The benchmark fills
:class:`~tda.allocate.descent.DescentProblem` from the reference field's
truncation defect; the plan fills it from the band probe's estimate
:math:`\\hat{\\mathbf v}` on a pilot arc the controller computed itself.  The
caller builds it, so this module cannot enforce which was used --- it records
the answer on the plan instead, where a serialized run states it.

The plan is small, and the reason it is small matters
-----------------------------------------------------
The online decision never needs :math:`\\mathbf c_i` in the state space.  It
contracts it immediately against the transport, so what the vehicle carries is

.. math::
    \\mathbf z_i \\;=\\; \\mathbf S_i^{\\top}\\mathbf c_i^{(-q)},\\qquad
    \\mathbf S_i=\\Delta t_i\\,\\Phi(t_0,m_i)\\mathbf B ,

a three-vector per cell: the direction in *acceleration* space that the
already-accumulated displacement asks this cell to cancel.  It is the same
:math:`\\mathbf z` the retained-fraction diagnostic is written in
(:func:`tda.probe.retained_fraction`), which is not a coincidence --- that
diagnostic asks how much of this contraction the probe recovers.

The exclusion is not cosmetic
-----------------------------
:math:`\\mathbf c_i^{(-q)}` **excludes** the interval being decided.  The online
score computes the within-block quadratic explicitly, so a target that still
contained the block's own planned contribution would count it twice, with a
sign that depends on the schedule --- an error no norm and no shape check can
see.  The exclusion is pre-computed here, once, from the same decomposition the
sweep uses (D162).
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

from tda.allocate.budget import ScheduleSolution
from tda.allocate.descent import DescentProblem, solve_descent, sweep_once
from tda.kernel import CouplingKernel

__all__ = ["PLAN_SOURCES", "OfflinePlan", "build_plan", "cancellation_target"]

Arr = NDArray[np.float64]
IntArr = NDArray[np.int_]

PLAN_SOURCES = ("pilot-probe", "reference", "pilot-truth")
"""What filled the contribution tensor, declared by the caller.

``pilot-probe``
    The deployable case: band-probe estimates on a pilot arc.
``reference``
    The benchmark's own tensor.  A plan built from it is an oracle and is only
    ever used as the upper anchor of the capture fraction.
``pilot-truth``
    The reference field evaluated on the pilot arc --- the ``abl-probe``
    control, which separates the cost of the pilot arc from the cost of the
    probe.
"""


def cancellation_target(kernel: CouplingKernel, transport: Arr, u: Arr,
                        cell_slices: IntArr, interval_of: IntArr) -> Arr:
    """:math:`\\mathbf z_i=\\mathbf S_i^{\\top}\\mathbf c_i^{(-q)}`, ``(M, 3)``.

    Parameters
    ----------
    kernel:
        Supplies :math:`\\mathbf A_i` and :math:`T`.
    transport:
        :math:`\\mathbf S_i=\\Delta t_i\\,\\Phi(t_0,m_i)\\mathbf B`, shape
        ``(M, 6, 3)`` --- the cell width folded in, so that the online score
        never has to reapply it.
    u:
        The planned contributions :math:`\\mathbf u_i`, shape ``(M, 6)``.
    cell_slices:
        Each interval's ``[start, stop)`` in cell order, shape ``(K, 2)``.
    interval_of:
        :math:`g(i)`, shape ``(M,)``.

    Returns
    -------
    ndarray, shape (M, 3)

    Raises
    ------
    ValueError
        If the shapes disagree.

    Notes
    -----
    Built from the two-vector collapse of D162:
    :math:`\\mathbf c_i^{(-q)}=T^{-1}[\\mathbf A_i\\mathbf P_{\\text{left}}
    +\\mathbf R_{\\text{right}}]` with
    :math:`\\mathbf P_{\\text{left}}=\\sum_{k<\\min I_q}\\mathbf u_k` and
    :math:`\\mathbf R_{\\text{right}}=\\sum_{k>\\max I_q}\\mathbf A_k\\mathbf
    u_k`, so all :math:`M` of them cost one prefix scan and one suffix scan
    rather than a pass over the interval structure per cell.
    """
    n_cells = len(kernel)
    transport = np.asarray(transport, dtype=float)
    if transport.shape != (n_cells, 6, 3):
        raise ValueError(
            f"transport must have shape ({n_cells}, 6, 3), got "
            f"{transport.shape}")
    u = np.asarray(u, dtype=float)
    if u.shape != (n_cells, 6):
        raise ValueError(f"u must have shape ({n_cells}, 6), got {u.shape}")

    starts, stops = cell_slices[:, 0], cell_slices[:, 1]
    # P_left[q] = sum of u over every cell strictly before interval q.
    running = np.zeros((n_cells + 1, 6), dtype=float)
    running[1:] = np.cumsum(u, axis=0)
    p_left = running[starts]

    # R_right[q] = sum_{k >= stop_q} A_k u_k.
    au = np.einsum("iab,ib->ia", kernel.suffix, u)
    tail = np.zeros((n_cells + 1, 6), dtype=float)
    tail[:-1] = np.cumsum(au[::-1], axis=0)[::-1]
    r_right = tail[stops]

    outer = (np.einsum("iab,ib->ia", kernel.suffix, p_left[interval_of])
             + r_right[interval_of]) / kernel.duration
    return np.einsum("iba,ib->ia", transport, outer)


@dataclass(frozen=True, slots=True)
class OfflinePlan:
    """Everything the vehicle carries, indexed by orbital phase.

    Attributes
    ----------
    coordinate:
        ``revolution + phase`` of each accumulation cell's midpoint, shape
        ``(M,)``, strictly increasing.
    edge_coordinate:
        The same for the decision-interval boundaries, shape ``(K+1,)``.
    interval_of:
        :math:`g(i)`, shape ``(M,)``.
    cell_slices:
        Each interval's ``[start, stop)``, shape ``(K, 2)``.
    transport:
        :math:`\\mathbf S_i`, shape ``(M, 6, 3)``.
    kernel:
        The pilot arc's coupling structure; supplies :math:`\\mathbf A_i` for
        the within-block quadratic and :math:`T`.
    target:
        :math:`\\mathbf z_i`, shape ``(M, 3)``, with the decided interval
        already excluded.
    time_weight:
        :math:`W_q`, shape ``(K,)``.
    nominal_degrees:
        :math:`N_{\\mathrm{plan}}`, shape ``(K,)``; the centre of the online
        search window.
    candidate_degrees:
        The tabulated degree axis.
    multiplier:
        :math:`\\lambda_0`.  Exactly zero is legitimate and means the ceiling
        did not bind on the pilot arc (D142); the online feedback may still
        raise it.
    objective:
        :math:`J` of the plan *on the pilot arc*, in the plan's own
        information.  Not comparable with the benchmark's :math:`J`, and never
        quoted against it.
    work, budget:
        Realized nominal work and the ceiling it was calibrated to.
    source:
        One of :data:`PLAN_SOURCES`.
    iterations:
        :math:`J`, the sweep cap the plan was built with.  ``abl-J`` sets it
        to one.
    converged:
        Whether one further sweep at the winning multiplier would move
        nothing.  Measured, not assumed: Q4 asks how many passes the iteration
        needs and a plan that cannot answer it does not close the question.
    n_starts:
        How many starting schedules were used.  The declared algorithm uses
        one, the sensitivity-weighted schedule; more is admissible and is
        recorded because it makes the plan a different algorithm.
    spread:
        Objective across those starts, relative to the best.
    """

    coordinate: Arr
    edge_coordinate: Arr
    interval_of: IntArr
    cell_slices: IntArr
    transport: Arr
    kernel: CouplingKernel
    target: Arr
    time_weight: Arr
    nominal_degrees: IntArr
    candidate_degrees: tuple[int, ...]
    multiplier: float
    objective: float
    work: float
    budget: float
    source: str
    iterations: int
    converged: bool
    n_starts: int = 1
    spread: tuple[float, ...] = ()
    diagnostics: dict[str, float] = field(default_factory=dict)

    def __len__(self) -> int:
        """Number of decision intervals, :math:`K_{\\mathrm{dec}}`."""
        return int(self.time_weight.size)

    @property
    def n_cells(self) -> int:
        """:math:`M`."""
        return int(self.coordinate.size)

    @property
    def nbytes(self) -> int:
        """Size of the arrays the vehicle would actually carry.

        Counts the phase-indexed tables, not the pilot arc they came from.
        The suffix blocks dominate at thirty-six numbers per cell; a flight
        build would drop them and reconstruct :math:`\\mathbf A_i` from the
        :math:`3\\times6` edge rows the kernel already holds, by the same
        reverse accumulation that built them.  That reduction is available and
        is not claimed here, because it has not been implemented or measured.
        """
        return int(self.transport.nbytes + self.kernel.suffix.nbytes
                   + self.target.nbytes + self.coordinate.nbytes
                   + self.time_weight.nbytes + self.nominal_degrees.nbytes)

    def locate(self, coordinate: float) -> int:
        """Decision interval containing a phase coordinate.

        Parameters
        ----------
        coordinate:
            ``revolution + phase`` on the *flown* arc.

        Returns
        -------
        int

        Raises
        ------
        IndexError
            If the coordinate is outside the plan's coverage.  Clamping would
            apply the last interval's cancellation target for the rest of the
            arc, which is a plausible-looking way to fly the wrong plan; the
            campaign counts these instead.
        """
        edges = self.edge_coordinate
        if not edges[0] <= coordinate <= edges[-1]:
            raise IndexError(
                f"phase coordinate {coordinate:.6f} is outside the plan's "
                f"coverage [{edges[0]:.6f}, {edges[-1]:.6f}]; the flown arc "
                "has drifted past the pilot arc it was planned on")
        return int(min(np.searchsorted(edges, coordinate, side="right") - 1,
                       len(self) - 1))

    def cell_fractions(self, interval: int) -> Arr:
        """Where each of an interval's cells sits within it, in :math:`[0,1]`.

        The fine half of the alignment.  The phase index says *which* planned
        interval the vehicle is in; this says where inside it each planned cell
        falls, so a probe taken at a fraction of the flown interval can be
        matched to the plan cells at the same fraction.  Comparing fractions
        rather than absolute phase is what keeps the match insensitive to the
        two arcs' intervals having slightly different phase widths.
        """
        lo, hi = self.cell_slices[interval]
        left = self.edge_coordinate[interval]
        width = self.edge_coordinate[interval + 1] - left
        return (self.coordinate[lo:hi] - left) / width


def _verify_transport(contributions: Arr, transport: Arr,
                      tolerance: float = 1.0e-8) -> None:
    """Check the transport against the tensor it is supposed to have built.

    Every :math:`\\mathbf u_i(N)` is :math:`\\mathbf S_i` applied to a
    three-vector, so it must lie in the column space of :math:`\\mathbf S_i`.
    Passing the wrong transport --- a different arc's, or one built at
    different epochs --- leaves every array the right shape and every number
    finite, and produces an online score that is a quadratic form in the wrong
    subspace.  The residual of the projection catches it in one batched QR.

    It does **not** catch a transport that is merely mis-scaled, since scaling
    does not move the column space.  That is why :func:`build_plan` folds the
    cell widths itself instead of accepting them already folded.

    Raises
    ------
    ValueError
        If any contribution lies outside its transport's range.
    """
    basis, _ = np.linalg.qr(transport)
    residual = contributions - np.einsum(
        "iab,ipb->ipa", basis, np.einsum("iba,ipb->ipa", basis, contributions))
    size = np.linalg.norm(contributions, axis=-1)
    error = np.linalg.norm(residual, axis=-1)
    live = size > 0.0
    worst = float(np.max(error[live] / size[live])) if live.any() else 0.0
    if worst > tolerance:
        raise ValueError(
            f"the contribution tensor does not lie in the transport's column "
            f"space (worst relative residual {worst:.3e} against a tolerance "
            f"of {tolerance:.3e}); the transport and the problem were not "
            "built from the same arc")


def build_plan(problem: DescentProblem, node_transport: Arr, widths: Arr,
               coordinate: Arr, edge_coordinate: Arr, budget: float,
               starts: Sequence[IntArr], *, source: str = "pilot-probe",
               iterations: int = 2) -> OfflinePlan:
    """Run Algorithm 2 and reduce its output to what the vehicle carries.

    Parameters
    ----------
    problem:
        Assembled from the **pilot** arc.  Whether its contribution tensor
        came from the probe or from the reference field is the whole content
        of RQ3 and is declared through ``source``.
    node_transport:
        :math:`\\Phi(t_0,m_i)\\mathbf B`, shape ``(M, 6, 3)``, *without* the
        cell width.  The width is folded here rather than by the caller
        because a transport that is right except for a missing
        :math:`\\Delta t_i` passes every structural check there is --- the
        column space is unchanged --- and silently rescales the online score
        against the plan it is supposed to agree with.
    widths:
        :math:`\\Delta t_i`, shape ``(M,)``, the same array the contribution
        tensor was built with.
    coordinate, edge_coordinate:
        Phase coordinates of the cell midpoints and the decision edges, from
        :meth:`tda.controller.phase.RevolutionIndex.coordinate`.
    budget:
        The work ceiling.
    starts:
        Starting schedules.  Algorithm 2 declares one --- the separable
        sensitivity-weighted schedule --- and passing more is recorded rather
        than refused.
    source:
        One of :data:`PLAN_SOURCES`.
    iterations:
        :math:`J`, the sweep cap.  ``abl-J`` sets it to one.

    Returns
    -------
    OfflinePlan

    Raises
    ------
    ValueError
        If ``source`` is not declared, ``iterations`` is not positive, the
        coordinate arrays do not match the grid or are not increasing --- a
        non-monotone phase coordinate would make the interval lookup silently
        return the wrong interval rather than fail --- or the transport does
        not match the contribution tensor.

    Notes
    -----
    Convergence is measured by running one further sweep at the winning
    multiplier and asking whether anything moves.  That costs one sweep and
    turns "``J=2`` was enough" from an assumption into a recorded fact.
    """
    if source not in PLAN_SOURCES:
        raise ValueError(
            f"source must be one of {PLAN_SOURCES}, got {source!r}")
    if iterations <= 0:
        raise ValueError(f"iterations must be positive, got {iterations}")
    coordinate = np.asarray(coordinate, dtype=float)
    edge_coordinate = np.asarray(edge_coordinate, dtype=float)
    n_cells = problem.contributions.shape[0]
    if coordinate.shape != (n_cells,):
        raise ValueError(
            f"coordinate must have shape ({n_cells},), got {coordinate.shape}")
    if edge_coordinate.shape != (problem.n_intervals + 1,):
        raise ValueError(
            f"edge_coordinate must have shape ({problem.n_intervals + 1},), "
            f"got {edge_coordinate.shape}")
    for name, values in (("coordinate", coordinate),
                         ("edge_coordinate", edge_coordinate)):
        if np.any(np.diff(values) <= 0.0):
            raise ValueError(
                f"{name} must be strictly increasing; a non-monotone phase "
                "coordinate makes the interval lookup wrong rather than loud")

    node_transport = np.asarray(node_transport, dtype=float)
    if node_transport.shape != (n_cells, 6, 3):
        raise ValueError(
            f"node_transport must have shape ({n_cells}, 6, 3), got "
            f"{node_transport.shape}")
    widths = np.asarray(widths, dtype=float)
    if widths.shape != (n_cells,):
        raise ValueError(
            f"widths must have shape ({n_cells},), got {widths.shape}")
    transport = node_transport * widths[:, None, None]
    _verify_transport(problem.contributions, transport)

    solution: ScheduleSolution = solve_descent(
        problem, budget, starts, max_sweeps=iterations)

    index = {d: p for p, d in enumerate(problem.candidate_degrees)}
    columns = np.array([index[int(n)] for n in solution.degrees],
                       dtype=np.int64)
    _, changed = sweep_once(problem, columns.copy(), solution.multiplier)

    u = problem.gather(solution.degrees)
    target = cancellation_target(problem.kernel, transport, u,
                                 problem.cell_slices, problem.interval_of)

    return OfflinePlan(
        coordinate=coordinate,
        edge_coordinate=edge_coordinate,
        interval_of=np.asarray(problem.interval_of),
        cell_slices=np.asarray(problem.cell_slices),
        transport=transport,
        kernel=problem.kernel,
        target=target,
        time_weight=np.asarray(problem.time_weight, dtype=float),
        nominal_degrees=np.asarray(solution.degrees),
        candidate_degrees=problem.candidate_degrees,
        multiplier=solution.multiplier,
        objective=solution.objective,
        work=solution.work,
        budget=solution.budget,
        source=source,
        iterations=iterations,
        converged=not changed,
        n_starts=len(starts),
        spread=solution.spread,
        diagnostics=dict(solution.diagnostics),
    )
