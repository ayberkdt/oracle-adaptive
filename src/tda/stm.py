"""Algebra of the state-transition matrix.

Pure linear algebra on :math:`6\\times6` blocks: no gravity model, no
integrator, no configuration.  It is a separate module because it has no
dependency on either and can therefore be verified against closed-form
identities alone, which is what makes a failure here unambiguous.  A test that
exercised the inverse *through* a propagation could not say whether a
discrepancy came from the algebra or from the integration.

Why the symplectic route
------------------------
The allocation needs the backward transport :math:`\\Phi(t_0,t_i)` at every
accumulation epoch, and an integration produces the forward one.  Inverting
each matrix numerically would work, but the flow of a time-dependent
Hamiltonian is symplectic [Hairer2006]_, so the inverse is a rearrangement of
transposed blocks [Montenbruck2000]_ -- no factorisation in the inner loop.

**What that does and does not buy.**  The *exact* flow is symplectic; DOP853
is not a symplectic integrator, so the computed
:math:`\\Phi_h` satisfies :math:`\\Phi_h^\\top\\mathbf J\\Phi_h=\\mathbf J+
\\varepsilon` and :math:`-\\mathbf J\\Phi_h^\\top\\mathbf J` is therefore
*not* its exact inverse.  It is expected to be close, but that is a claim to
measure rather than assert.  :func:`inverse_residual` compares the identity
against a direct solve and :func:`condition_number` records how much room
there is between them; the campaign reports both per arc rather than saying
the conditioning question went away.

References
----------
.. [Montenbruck2000] O. Montenbruck and E. Gill, *Satellite Orbits*, Springer,
   2000, §7.1 -- state-transition matrices and the block-transpose inverse of
   a symplectic one.
.. [Hairer2006] E. Hairer, C. Lubich, G. Wanner, *Geometric Numerical
   Integration*, 2nd ed., Springer, 2006, Ch. VI -- symplecticity of
   Hamiltonian flows and its use as a numerical diagnostic.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

__all__ = [
    "STATE_DIM",
    "SYMPLECTIC_J",
    "condition_number",
    "inverse_residual",
    "position_rows",
    "symplectic_defect",
    "symplectic_inverse",
    "velocity_transport",
]

Arr = NDArray[np.float64]

STATE_DIM = 6
"""Dimension of the Cartesian state, ordered ``(r, v)``."""

SYMPLECTIC_J: Arr = np.block([
    [np.zeros((3, 3)), np.eye(3)],
    [-np.eye(3), np.zeros((3, 3))],
])
"""The standard symplectic form :math:`\\mathbf J`, ordering ``(r, v)``."""


def _check_square(phi: Arr) -> Arr:
    """Coerce to float and require a trailing ``(6, 6)``."""
    phi = np.asarray(phi, dtype=float)
    if phi.shape[-2:] != (STATE_DIM, STATE_DIM):
        raise ValueError(f"expected (...,6,6), got {phi.shape}")
    return phi


def _nondimensionalise(phi: Arr, length_scale: float,
                       time_scale: float) -> Arr:
    """Return S^-1 phi S with S = diag(L,L,L,L/T,L/T,L/T).

    The raw state-transition matrix mixes units across its blocks, so any
    norm or condition number taken on it is meaningless.  This scaling leaves
    the symplectic condition invariant because S scales J by a single
    factor.
    """
    if length_scale <= 0.0 or time_scale <= 0.0:
        raise ValueError(
            f"scales must be positive, got L={length_scale}, T={time_scale}")
    scale = np.concatenate([
        np.full(3, length_scale),
        np.full(3, length_scale / time_scale),
    ])
    return phi / scale[:, None] * scale[None, :]


def symplectic_inverse(phi: Arr) -> Arr:
    """Invert a symplectic state-transition matrix by block transposition.

    For :math:`\\Phi=\\bigl(\\begin{smallmatrix}A&B\\\\C&D\\end{smallmatrix}
    \\bigr)` satisfying :math:`\\Phi^\\top\\mathbf J\\Phi=\\mathbf J`,

    .. math::
        \\Phi^{-1} \\;=\\; -\\mathbf J\\,\\Phi^\\top\\mathbf J
        \\;=\\; \\begin{pmatrix} D^\\top & -B^\\top \\\\
                                -C^\\top & A^\\top \\end{pmatrix},

    a standard identity [Montenbruck2000]_.

    Parameters
    ----------
    phi:
        Shape ``(6, 6)`` or ``(M, 6, 6)``.  A batch is transposed in one pass
        rather than looped.

    Returns
    -------
    ndarray
        Same shape as ``phi``.

    Raises
    ------
    ValueError
        If the trailing shape is not ``(6, 6)``.

    Warnings
    --------
    Exact only for a symplectic matrix.  This is **not** a general inverse and
    returns a plausible wrong answer for one, silently.  Check
    :func:`symplectic_defect` before relying on it; the campaign does so once
    per arc and records the result.
    """
    phi = _check_square(phi)

    def transposed(block: Arr) -> Arr:
        """Transpose the trailing two axes, leaving any batch axes alone."""
        return np.swapaxes(block, -1, -2)

    out = np.empty_like(phi)
    out[..., 0:3, 0:3] = transposed(phi[..., 3:6, 3:6])   #  D^T
    out[..., 0:3, 3:6] = -transposed(phi[..., 0:3, 3:6])  # -B^T
    out[..., 3:6, 0:3] = -transposed(phi[..., 3:6, 0:3])  # -C^T
    out[..., 3:6, 3:6] = transposed(phi[..., 0:3, 0:3])   #  A^T
    return out


def symplectic_defect(phi: Arr, length_scale: float,
                      time_scale: float) -> Arr:
    """Dimensionless residual of :math:`\\Phi^\\top\\mathbf J\\Phi=\\mathbf J`.

    The exact flow preserves the symplectic form [Hairer2006]_, so any
    residual is integration error -- provided the gravity gradient that
    generated the matrix was symmetrised, which
    :func:`tda.field.gravity_gradient` does, or the residual would read out
    differencing asymmetry instead.

    Scaling
    -------
    The raw matrix is dimensionally inhomogeneous -- its position--velocity
    blocks carry seconds and inverse seconds -- so a norm taken on it mixes
    units and means nothing.  It is first non-dimensionalised by
    :math:`\\mathbf S=\\operatorname{diag}(L,L,L,L/T,L/T,L/T)`,
    :math:`\\tilde\\Phi=\\mathbf S^{-1}\\Phi\\mathbf S`, which leaves the
    symplectic condition invariant because :math:`\\mathbf S` scales
    :math:`\\mathbf J` by a single factor.

    Parameters
    ----------
    phi:
        Shape ``(6, 6)`` or ``(M, 6, 6)``.
    length_scale, time_scale:
        Characteristic scales, e.g. the initial radius and radius over speed.
        Any consistent pair works.

    Returns
    -------
    ndarray
        Scalar for a single matrix, shape ``(M,)`` for a batch: the maximum
        absolute entry of :math:`\\tilde\\Phi^\\top\\mathbf J\\tilde\\Phi
        -\\mathbf J`.

    Raises
    ------
    ValueError
        If the trailing shape is not ``(6, 6)``, or a scale is not positive.

    Examples
    --------
    >>> import numpy as np
    >>> float(symplectic_defect(np.eye(6), 1.0, 1.0))
    0.0
    """
    tilde = _nondimensionalise(_check_square(phi), length_scale, time_scale)
    residual = (np.swapaxes(tilde, -1, -2) @ SYMPLECTIC_J @ tilde
                - SYMPLECTIC_J)
    return np.max(np.abs(residual), axis=(-2, -1))


def velocity_transport(phi_inverse: Arr) -> Arr:
    """Extract :math:`\\Phi(t_0,t_i)\\mathbf B`, the columns acceleration hits.

    :math:`\\mathbf B` maps an acceleration into the velocity block of the
    state, so :math:`\\Phi\\mathbf B` is the last three columns of
    :math:`\\Phi`.  Building :math:`\\mathbf B` and multiplying would allocate
    a 6x3 matrix and do a 6x6 by 6x3 product per epoch to reproduce a slice.

    Parameters
    ----------
    phi_inverse:
        :math:`\\Phi(t_0,t_i)`, shape ``(...,6,6)`` -- typically the output of
        :func:`symplectic_inverse`.

    Returns
    -------
    ndarray, shape ``(...,6,3)``
        Multiply by :math:`\\Delta\\mathbf a\\,\\Delta t` to obtain
        :math:`\\mathbf u_i`.
    """
    return _check_square(phi_inverse)[..., :, 3:6]


def position_rows(phi: Arr) -> Arr:
    """Extract :math:`\\mathbf H_r\\Phi`, the position rows.

    The objective selects position before taking a norm, and
    :math:`\\mathbf H_r=[\\,\\mathbf I_3\\;\\;\\mathbf 0\\,]` is a row slice.
    The per-epoch block of the objective is
    :math:`\\mathbf M_j=(\\mathbf H_r\\Phi_j)^\\top(\\mathbf H_r\\Phi_j)`,
    which :mod:`tda.kernel` forms from this.

    Parameters
    ----------
    phi:
        Shape ``(...,6,6)``.

    Returns
    -------
    ndarray, shape ``(...,3,6)``
    """
    return _check_square(phi)[..., 0:3, :]


def inverse_residual(phi: Arr, length_scale: float,
                     time_scale: float) -> Arr:
    """How far the symplectic identity is from the true inverse.

    :math:`\\lVert\\tilde\\Phi^{-1}_{\\mathrm{symp}}\\tilde\\Phi-
    \\mathbf I\\rVert_\\infty`, computed on the non-dimensionalised matrix so
    that the norm does not mix units.  Zero for an exactly symplectic matrix;
    for one produced by a non-symplectic integrator it is the price of using
    the identity instead of a factorisation, and the campaign records it
    rather than assuming it away.

    Parameters
    ----------
    phi:
        Shape ``(6, 6)`` or ``(M, 6, 6)``.
    length_scale, time_scale:
        As for :func:`symplectic_defect`.

    Returns
    -------
    ndarray
        Scalar or shape ``(M,)``.

    See Also
    --------
    condition_number : how much room a direct solve would have had.

    Examples
    --------
    >>> import numpy as np
    >>> float(inverse_residual(np.eye(6), 1.0, 1.0))
    0.0
    """
    tilde = _nondimensionalise(_check_square(phi), length_scale, time_scale)
    product = symplectic_inverse(tilde) @ tilde
    return np.max(np.abs(product - np.eye(STATE_DIM)), axis=(-2, -1))


def condition_number(phi: Arr, length_scale: float,
                     time_scale: float) -> Arr:
    """Spectral condition number of the non-dimensionalised matrix.

    Recorded alongside :func:`inverse_residual` so that a small residual can
    be read as "the identity agrees with a solve" rather than as "the matrix
    was too well conditioned for the comparison to mean anything".  Over a
    seven-day low lunar arc the along-track block grows secularly, so this is
    expected to be large and is worth watching.

    Parameters
    ----------
    phi:
        Shape ``(6, 6)`` or ``(M, 6, 6)``.
    length_scale, time_scale:
        As for :func:`symplectic_defect`.  Non-dimensionalising first is
        essential: the condition number of the raw matrix is dominated by the
        unit mismatch between its blocks and says nothing about the operator.

    Returns
    -------
    ndarray
        Scalar or shape ``(M,)``.
    """
    tilde = _nondimensionalise(_check_square(phi), length_scale, time_scale)
    return np.linalg.cond(tilde)
