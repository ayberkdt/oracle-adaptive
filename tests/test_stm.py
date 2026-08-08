"""Verification of the state-transition-matrix algebra.

Pure linear algebra, so the assertions are closed-form identities rather than
tolerances: a failure here cannot be blamed on an integrator or a field.

The pair :func:`test_free_drift_shear_is_symplectic` and
:func:`test_symplectic_defect_detects_a_broken_matrix` is the point of this
file.  Perturbing a matrix "obviously wrongly" is not enough -- most
single-entry perturbations of the identity are still symplectic -- so the
broken cases are chosen against the criterion instead of by eye.
"""

from __future__ import annotations

import numpy as np
import pytest

from tda.stm import (
    SYMPLECTIC_J,
    condition_number,
    inverse_residual,
    position_rows,
    symplectic_defect,
    symplectic_inverse,
    velocity_transport,
)


def _random_symplectic(rng: np.random.Generator) -> np.ndarray:
    """Build a symplectic matrix as the exponential of a Hamiltonian one.

    ``H = J S`` with ``S`` symmetric is Hamiltonian, and ``expm(H)`` is
    symplectic.  A small scaling keeps the exponential well conditioned so
    the test measures the identity rather than the matrix exponential.
    """
    from scipy.linalg import expm

    s = rng.normal(size=(6, 6))
    s = 0.5 * (s + s.T)
    return expm(0.05 * SYMPLECTIC_J @ s)


def test_symplectic_inverse_is_an_inverse() -> None:
    rng = np.random.default_rng(20260808)
    for _ in range(8):
        phi = _random_symplectic(rng)
        assert np.allclose(symplectic_inverse(phi) @ phi, np.eye(6),
                           rtol=0.0, atol=1e-12)


def test_symplectic_inverse_is_batched() -> None:
    rng = np.random.default_rng(1)
    batch = np.stack([_random_symplectic(rng) for _ in range(5)])
    inv = symplectic_inverse(batch)
    assert inv.shape == batch.shape
    assert np.allclose(inv @ batch, np.eye(6)[None], rtol=0.0, atol=1e-12)


def test_symplectic_defect_vanishes_on_identity() -> None:
    assert symplectic_defect(np.eye(6), 1.0, 1.0) == 0.0


def test_free_drift_shear_is_symplectic() -> None:
    """Not every off-diagonal term breaks the form, and it matters which.

    ``I + tau*e_0 e_3^T`` is one component of free drift, and free drift is a
    Hamiltonian flow, so the residual vanishes identically -- at first order
    *and* at second, since ``E^T J E = 0`` for ``E = e_0 e_3^T``.  Recorded as
    a test because it is the trap a hand-written "obviously wrong" matrix
    falls into: perturbing a position row towards its *own* velocity column
    produces a matrix that is still symplectic.
    """
    shear = np.eye(6)
    shear[0, 3] = 1.0e-3
    assert symplectic_defect(shear, 1.0, 1.0) == pytest.approx(0.0, abs=1e-15)


@pytest.mark.parametrize(("index", "value", "why"), [
    ((0, 0), 1.1, "position scaled with no compensating velocity contraction"),
    ((0, 4), 1.0e-3, "x coupled to v_y; J e_0 = -e_3, so e_0 e_4^T is not "
                     "Hamiltonian"),
    ((1, 3), 1.0e-3, "y coupled to v_x; J e_1 = -e_4, so e_1 e_3^T is not "
                     "Hamiltonian"),
])
def test_symplectic_defect_detects_a_broken_matrix(index, value, why) -> None:
    """A matrix that is not symplectic must not pass silently.

    Guards the warning in :func:`symplectic_inverse`: the block-transpose
    identity returns a wrong answer for a general matrix, so the defect is the
    only thing standing between a bad integration and a plausible number.

    The perturbations are chosen against the criterion rather than by eye.
    For :math:`\\Phi=\\mathbf I+\\varepsilon\\,\\mathbf e_a\\mathbf e_b^\\top`
    the first-order residual vanishes exactly when
    :math:`\\mathbf J\\mathbf e_a\\propto\\mathbf e_b`, and the second-order
    term :math:`\\mathbf E^\\top\\mathbf J\\mathbf E` vanishes for every
    rank-one :math:`\\mathbf E`.  Picking a pair that satisfies the
    proportionality would therefore give a *symplectic* matrix and a test that
    silently proves nothing -- which is what
    :func:`test_free_drift_shear_is_symplectic` records.
    """
    broken = np.eye(6)
    broken[index] = value
    assert symplectic_defect(broken, 1.0, 1.0) > 1.0e-6, why


def test_symplectic_defect_is_scale_invariant() -> None:
    """The dimensionless residual must not depend on the chosen units."""
    rng = np.random.default_rng(7)
    phi = _random_symplectic(rng)
    # Re-express the same operator in units where length and time differ.
    scale = np.concatenate([np.full(3, 1.0e6), np.full(3, 1.0e6 / 5.0e3)])
    dimensional = phi * scale[:, None] / scale[None, :]
    a = symplectic_defect(phi, 1.0, 1.0)
    b = symplectic_defect(dimensional, 1.0e6, 5.0e3)
    # Both are roundoff-level; what is asserted is that re-expressing the same
    # operator in other units does not manufacture a defect. Without the
    # non-dimensionalisation, ``b`` would be of order the scale factors.
    assert b == pytest.approx(a, rel=1e-6, abs=1e-12)
    assert b < 1.0e-12


# ---------------------------------------------------------------------------
# Transport helpers
# ---------------------------------------------------------------------------


def test_velocity_transport_selects_the_acceleration_columns() -> None:
    rng = np.random.default_rng(3)
    phi = rng.normal(size=(4, 6, 6))
    b = np.vstack([np.zeros((3, 3)), np.eye(3)])
    assert np.array_equal(velocity_transport(phi), phi @ b)


def test_position_rows_selects_the_position_block() -> None:
    rng = np.random.default_rng(4)
    phi = rng.normal(size=(4, 6, 6))
    h_r = np.hstack([np.eye(3), np.zeros((3, 3))])
    assert np.array_equal(position_rows(phi), h_r @ phi)


# ---------------------------------------------------------------------------
# The identity against a direct solve
# ---------------------------------------------------------------------------


def test_inverse_residual_vanishes_on_an_exactly_symplectic_matrix() -> None:
    rng = np.random.default_rng(31)
    for _ in range(5):
        phi = _random_symplectic(rng)
        assert inverse_residual(phi, 1.0, 1.0) < 1e-12


def test_inverse_residual_is_nonzero_when_the_matrix_is_not_symplectic() -> None:
    """DOP853 is not symplectic, so this quantity is not identically zero.

    The module claims the block-transpose identity avoids a factorisation, not
    that it is exact for a numerically integrated matrix.  This pins that the
    diagnostic can actually see the difference.
    """
    broken = np.eye(6)
    broken[0, 0] = 1.1
    assert inverse_residual(broken, 1.0, 1.0) > 1e-3


def test_symplectic_inverse_agrees_with_a_direct_solve() -> None:
    """The comparison the campaign records per arc, in miniature."""
    rng = np.random.default_rng(32)
    phi = _random_symplectic(rng)
    direct = np.linalg.solve(phi, np.eye(6))
    assert np.allclose(symplectic_inverse(phi), direct, rtol=0.0, atol=1e-11)


def test_condition_number_is_reported_on_the_scaled_matrix() -> None:
    """Unscaled, the condition number reads out the unit mismatch.

    A pure drift over 1800 s has a modest condition number in natural units
    and an enormous one in SI, and only the first says anything about the
    operator.
    """
    drift = np.eye(6)
    drift[0:3, 3:6] = 1800.0 * np.eye(3)
    scaled = condition_number(drift, 1.8e6, 1.0e3)
    unscaled = float(np.linalg.cond(drift))
    assert scaled < unscaled / 100.0


def test_condition_number_of_the_identity_is_one() -> None:
    assert float(condition_number(np.eye(6), 1.0, 1.0)) == pytest.approx(1.0)


@pytest.mark.parametrize("fn", [inverse_residual, condition_number,
                                symplectic_defect])
def test_diagnostics_reject_a_non_positive_scale(fn) -> None:
    with pytest.raises(ValueError, match="scales must be positive"):
        fn(np.eye(6), 0.0, 1.0)
