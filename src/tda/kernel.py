"""The objective's coupling structure.

The mathematical core: the suffix sequence that generates :math:`\\mathbf Q`,
the objective value, the cross terms a coordinate method needs, and the local
sensitivity kernel.  Pure array arithmetic --- no gravity model, no
integrator, no configuration --- so every claim here is checkable against a
brute-force evaluation, and :mod:`tests.test_kernel` does exactly that.

Index convention, which is load-bearing
---------------------------------------
Two index sets, deliberately distinct:

* **cells** :math:`i=0,\\dots,M-1`, spanning :math:`[e_i,e_{i+1}]`, with the
  defect sampled at the midpoint and width :math:`\\Delta t_i`;
* **edges** :math:`j=0,\\dots,M`, where the objective is sampled.

Cell :math:`i` is felt at edge :math:`j` exactly when :math:`j\\ge i+1` --- at
its own left edge the impulse has not happened yet.  Hence

.. math::
    \\mathbf S_j=\\sum_{i\\le j-1}\\mathbf u_i , \\qquad
    \\mathbf A_i=\\sum_{j=i+1}^{M}\\omega_j\\mathbf M_j , \\qquad
    \\mathbf Q_{ik}=\\tfrac1T\\mathbf A_{\\max(i,k)} .

**The suffix starts at** :math:`i+1`, **not at** :math:`i`.  Writing
:math:`\\sum_{j\\ge i}` adds :math:`\\omega_i\\mathbf M_i` to every block,
which credits a cell with displacing the state before its own impulse.  That
is not a small perturbation: on a random five-cell instance it misstates
:math:`J` by 26 per cent and breaks the identity between :math:`\\mathbf K_i`
and the diagonal of the objective kernel.  The single-grid form in earlier
drafts was self-consistent; separating cells from edges is what introduced the
offset (``DECISIONS.md`` D153).

What is never formed
--------------------
:math:`\\mathbf Q` is :math:`M\\times M` in :math:`6\\times6` blocks and is
never assembled: at fifty thousand cells that is :math:`10^{11}` numbers.
Everything is expressed through the suffix sequence, which is
:math:`M\\times6\\times6`, and through prefix and suffix scans over it.
:math:`\\mathbf M_j` is not formed either --- only its factor
:math:`\\mathbf H_r\\Phi(e_j,t_0)`, which is :math:`3\\times6`.

References
----------
.. [Neumaier1974] A. Neumaier, "Rundungsfehleranalyse einiger Verfahren zur
   Summation endlicher Summen", *ZAMM* 54, 1974 -- compensated summation, used
   for the prefix scans because their cancellation is the quantity of
   interest.
.. [Higham2002] N. J. Higham, *Accuracy and Stability of Numerical
   Algorithms*, 2nd ed., SIAM, 2002, Ch. 4 -- error bounds for summation and
   the condition number of a sum.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

__all__ = ["CouplingKernel", "compensated_prefix_sum"]

Arr = NDArray[np.float64]


def compensated_prefix_sum(terms: Arr) -> Arr:
    """Running sum of ``terms`` along axis 0, in Neumaier arithmetic.

    The prefix sums of :math:`\\mathbf u` are sums of order :math:`10^5`
    signed vectors whose *cancellation is the quantity the method exploits*,
    so the summation is ill-conditioned by construction: the condition number
    of a sum is :math:`\\sum|x_i|/|\\sum x_i|` [Higham2002]_, which here is
    exactly the ratio the campaign reports as a diagnostic.  Compensated
    summation [Neumaier1974]_ keeps the error at roundoff of the *result*
    rather than of the intermediate partial sums.

    Parameters
    ----------
    terms:
        Shape ``(M, ...)``.  Summed along the leading axis.

    Returns
    -------
    ndarray
        Same shape; entry ``i`` is the sum of ``terms[:i+1]``.

    Examples
    --------
    Catastrophic cancellation that a naive running sum gets wrong.

    >>> import numpy as np
    >>> x = np.array([1.0, 1e16, 1.0, -1e16])
    >>> float(compensated_prefix_sum(x)[-1])
    2.0
    >>> float(np.cumsum(x)[-1])
    0.0
    """
    terms = np.asarray(terms, dtype=float)
    out = np.empty_like(terms)
    total = np.zeros(terms.shape[1:], dtype=float)
    compensation = np.zeros_like(total)
    for i, term in enumerate(terms):
        candidate = total + term
        # Neumaier: the lost low-order part depends on which operand is larger.
        bigger = np.abs(total) >= np.abs(term)
        compensation += np.where(bigger,
                                 (total - candidate) + term,
                                 (term - candidate) + total)
        total = candidate
        out[i] = total + compensation
    return out


@dataclass(frozen=True, slots=True)
class CouplingKernel:
    """The suffix sequence :math:`\\mathbf A_i` and what it generates.

    Attributes
    ----------
    suffix:
        :math:`\\mathbf A_i=\\sum_{j\\ge i+1}\\omega_j\\mathbf M_j`, shape
        ``(M, 6, 6)``, symmetric positive semidefinite and decreasing in the
        Loewner order.
    duration:
        Arc length :math:`T`, the normalisation in
        :math:`\\mathbf Q_{ik}=\\mathbf A_{\\max(i,k)}/T`.
    edge_rows:
        :math:`\\mathbf H_r\\Phi(e_j,t_0)`, shape ``(M+1, 3, 6)``, retained so
        that the objective can be evaluated in :math:`\\mathcal O(M)` without
        going through :math:`\\mathbf Q`.
    outer_weights:
        :math:`\\omega_j\\ge0` on the edges, shape ``(M+1,)``.
    """

    suffix: Arr
    duration: float
    edge_rows: Arr
    outer_weights: Arr

    def __len__(self) -> int:
        """Number of cells, :math:`M`."""
        return int(self.suffix.shape[0])

    # -- construction -------------------------------------------------------

    @classmethod
    def from_arc(cls, edge_rows: Arr, outer_weights: Arr,
                 duration: float) -> CouplingKernel:
        """Build the suffix sequence from the position rows at the edges.

        Parameters
        ----------
        edge_rows:
            :math:`\\mathbf H_r\\Phi(e_j,t_0)` for every edge, shape
            ``(M+1, 3, 6)`` --- the output of
            :func:`tda.stm.position_rows` on the edge state-transition
            matrices.  Passing the :math:`3\\times6` factor rather than the
            :math:`6\\times6` product :math:`\\mathbf M_j` halves the memory
            and makes the positive semidefiniteness structural.
        outer_weights:
            :math:`\\omega_j\\ge0`, shape ``(M+1,)``, summing to ``duration``.
        duration:
            Arc length :math:`T`.

        Returns
        -------
        CouplingKernel

        Raises
        ------
        ValueError
            If the shapes disagree, a weight is negative -- which would break
            :math:`\\mathbf Q\\succeq0` and with it the convex relaxation --
            or the duration is not positive.
        """
        edge_rows = np.asarray(edge_rows, dtype=float)
        outer_weights = np.asarray(outer_weights, dtype=float)
        if edge_rows.ndim != 3 or edge_rows.shape[1:] != (3, 6):
            raise ValueError(
                f"edge_rows must have shape (M+1, 3, 6), got {edge_rows.shape}")
        if outer_weights.shape != (edge_rows.shape[0],):
            raise ValueError(
                f"outer_weights {outer_weights.shape} must match the "
                f"{edge_rows.shape[0]} edges")
        if np.any(outer_weights < 0.0):
            raise ValueError(
                "outer weights must be non-negative; a negative weight breaks "
                "the positive semidefiniteness the relaxation rests on")
        if duration <= 0.0:
            raise ValueError(f"duration must be positive, got {duration}")

        # weighted[j] = omega_j * M_j, formed from the 3x6 factor.
        weighted = outer_weights[:, None, None] * np.einsum(
            "jai,jak->jik", edge_rows, edge_rows)
        # A_i = sum over j >= i+1, so reverse-accumulate and drop edge 0.
        reverse = np.cumsum(weighted[::-1], axis=0)[::-1]
        return cls(suffix=np.ascontiguousarray(reverse[1:]),
                   duration=float(duration),
                   edge_rows=edge_rows,
                   outer_weights=outer_weights)

    # -- evaluations --------------------------------------------------------

    def prefix_states(self, u: Arr) -> Arr:
        """:math:`\\mathbf S_j` at every edge, shape ``(M+1, 6)``.

        ``S[0]`` is zero by construction: no cell has acted yet.
        """
        u = self._check_u(u)
        out = np.zeros((len(self) + 1, 6), dtype=float)
        out[1:] = compensated_prefix_sum(u)
        return out

    def objective(self, u: Arr) -> float:
        """:math:`J(\\mathbf u)=\\mathbf u^{\\top}\\mathbf Q\\mathbf u`.

        Evaluated as the quadrature it is,
        :math:`T^{-1}\\sum_j\\omega_j\\lVert\\mathbf H_r\\Phi_j
        \\mathbf S_j\\rVert^2`, which is :math:`\\mathcal O(M)` and never
        touches :math:`\\mathbf Q`.
        """
        states = self.prefix_states(u)
        displacement = np.einsum("jai,ji->ja", self.edge_rows, states)
        return float(self.outer_weights
                     @ np.einsum("ja,ja->j", displacement, displacement)
                     / self.duration)

    def cross_terms(self, u: Arr) -> Arr:
        """:math:`\\mathbf c_i=\\sum_k\\mathbf Q_{ik}\\mathbf u_k`.

        One prefix scan and one suffix scan:

        .. math::
            \\mathbf c_i=\\frac1T\\Bigl[\\mathbf A_i\\sum_{k\\le i}\\mathbf u_k
            +\\sum_{k>i}\\mathbf A_k\\mathbf u_k\\Bigr],

        so all :math:`M` of them cost :math:`\\mathcal O(M)` rather than the
        :math:`\\mathcal O(M^2)` a dense :math:`\\mathbf Q` would.  The
        gradient is :math:`\\nabla J=2\\mathbf c`.

        Returns
        -------
        ndarray, shape (M, 6)
        """
        u = self._check_u(u)
        prefix = compensated_prefix_sum(u)
        head = np.einsum("iab,ib->ia", self.suffix, prefix)
        au = np.einsum("iab,ib->ia", self.suffix, u)
        tail = np.zeros_like(au)
        tail[:-1] = np.cumsum(au[::-1], axis=0)[::-1][1:]
        return (head + tail) / self.duration

    def gradient(self, u: Arr) -> Arr:
        """:math:`\\nabla J=2\\mathbf Q\\mathbf u`, shape ``(M, 6)``."""
        return 2.0 * self.cross_terms(u)

    def local_kernel(self, node_transport: Arr) -> Arr:
        """:math:`\\mathbf K_i`, the objective kernel on its diagonal.

        .. math::
            \\mathbf K_i=\\frac1T\\,\\mathbf B^{\\top}\\Phi(t_0,m_i)^{\\top}
            \\mathbf A_i\\,\\Phi(t_0,m_i)\\,\\mathbf B ,

        the weight of the separable sensitivity criterion.  Read off the same
        suffix sequence, so it costs nothing beyond a small product.

        Parameters
        ----------
        node_transport:
            :math:`\\Phi(t_0,m_i)\\mathbf B`, shape ``(M, 6, 3)`` --- the
            output of :func:`tda.stm.velocity_transport` at the cell
            midpoints.

        Returns
        -------
        ndarray, shape (M, 3, 3)
            Symmetric positive semidefinite.
        """
        node_transport = np.asarray(node_transport, dtype=float)
        if node_transport.shape != (len(self), 6, 3):
            raise ValueError(
                f"node_transport must have shape ({len(self)}, 6, 3), got "
                f"{node_transport.shape}")
        return np.einsum("iba,ibc,icd->iad", node_transport, self.suffix,
                         node_transport) / self.duration

    def conditioning(self, u: Arr) -> float:
        """Condition number of the final prefix sum.

        :math:`\\sum_i\\lVert\\mathbf u_i\\rVert / \\lVert\\mathbf S_M\\rVert`
        [Higham2002]_.  Large is expected and is not a fault: it *is* the
        cancellation the method exploits.  It is reported per orbit so that a
        result can be read against how much cancellation produced it.
        """
        u = self._check_u(u)
        total = np.linalg.norm(u.sum(axis=0))
        if total == 0.0:
            return float("inf")
        return float(np.linalg.norm(u, axis=1).sum() / total)

    # -- internals ----------------------------------------------------------

    def _check_u(self, u: Arr) -> Arr:
        """Coerce and shape-check a stack of transported cell contributions."""
        u = np.asarray(u, dtype=float)
        if u.shape != (len(self), 6):
            raise ValueError(
                f"u must have shape ({len(self)}, 6), got {u.shape}")
        return u
