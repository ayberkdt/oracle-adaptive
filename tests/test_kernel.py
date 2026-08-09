"""Verification of the coupling structure against brute force.

The load-bearing tests here are the ones that build the objective twice: once
from its definition as a quadrature over edges, and once through the suffix
sequence, and require the two to agree.  That comparison is what caught the
index error of D153, in which the suffix started at the cell rather than at
the edge after it -- a 26 per cent misstatement of ``J`` that no structural
check would have noticed, because the wrong form is a perfectly well-formed
positive semidefinite quadratic.

Nothing here needs a gravity model or an integrator: the state-transition
matrices are random, which is stronger for this purpose than physical ones,
since a physical arc might accidentally satisfy an identity that the algebra
does not.
"""

from __future__ import annotations

import numpy as np
import pytest

from tda.kernel import CouplingKernel, compensated_prefix_sum

M_CELLS = 7
H_R = np.hstack([np.eye(3), np.zeros((3, 3))])
B_ACC = np.vstack([np.zeros((3, 3)), np.eye(3)])


@pytest.fixture
def instance():
    """A small random instance: edge rows, weights, duration, contributions."""
    rng = np.random.default_rng(20260809)
    phi = np.stack([np.eye(6)]
                   + [rng.normal(size=(6, 6)) for _ in range(M_CELLS)])
    edge_rows = np.einsum("ab,jbc->jac", H_R, phi)
    weights = rng.uniform(0.5, 2.0, M_CELLS + 1)
    u = rng.normal(size=(M_CELLS, 6))
    return edge_rows, weights, float(weights.sum()), u, phi


def _objective_from_definition(edge_rows, weights, duration, u) -> float:
    """J built directly from the quadrature, cell by cell and edge by edge."""
    total = 0.0
    for j in range(len(weights)):
        state = u[:j].sum(axis=0) if j > 0 else np.zeros(6)
        displacement = edge_rows[j] @ state
        total += weights[j] * displacement @ displacement
    return total / duration


# ---------------------------------------------------------------------------
# The index convention
# ---------------------------------------------------------------------------


def test_objective_matches_the_definition(instance) -> None:
    """The whole derivation, checked end to end against the quadrature."""
    edge_rows, weights, duration, u, _ = instance
    kernel = CouplingKernel.from_arc(edge_rows, weights, duration)
    assert kernel.objective(u) == pytest.approx(
        _objective_from_definition(edge_rows, weights, duration, u), rel=1e-12)


def test_dense_q_from_the_suffix_reproduces_the_objective(instance) -> None:
    """``u' Q u`` with ``Q_ik = A_max(i,k)/T`` must give the same number.

    Assembling ``Q`` is exactly what the implementation avoids, which is why
    it is done here: the cheap path is only worth having if it agrees with the
    expensive one.
    """
    edge_rows, weights, duration, u, _ = instance
    kernel = CouplingKernel.from_arc(edge_rows, weights, duration)

    total = 0.0
    for i in range(M_CELLS):
        for k in range(M_CELLS):
            total += u[i] @ kernel.suffix[max(i, k)] @ u[k]
    assert total / duration == pytest.approx(kernel.objective(u), rel=1e-12)


def test_the_off_by_one_suffix_is_detectably_wrong(instance) -> None:
    """``A_i = sum_(j>=i)`` is not a small perturbation of the right form.

    Pinned as a test because the wrong form is well-posed -- symmetric,
    positive semidefinite, structurally identical -- so nothing but a
    comparison against the definition rejects it.  It credits a cell with
    displacing the state at its own left edge, before its impulse.
    """
    edge_rows, weights, duration, u, _ = instance
    kernel = CouplingKernel.from_arc(edge_rows, weights, duration)
    truth = _objective_from_definition(edge_rows, weights, duration, u)

    weighted = weights[:, None, None] * np.einsum(
        "jai,jak->jik", edge_rows, edge_rows)
    wrong_suffix = np.cumsum(weighted[::-1], axis=0)[::-1][:-1]   # j >= i
    wrong = sum(u[i] @ wrong_suffix[max(i, k)] @ u[k]
                for i in range(M_CELLS) for k in range(M_CELLS)) / duration

    assert kernel.objective(u) == pytest.approx(truth, rel=1e-12)
    assert abs(wrong / truth - 1.0) > 0.05, "the two forms must differ visibly"


def test_first_cell_is_not_felt_at_the_first_edge(instance) -> None:
    """Causality, stated on the smallest case that can express it.

    Perturbing cell 0 must leave the displacement at edge 0 untouched.
    """
    edge_rows, weights, duration, u, _ = instance
    kernel = CouplingKernel.from_arc(edge_rows, weights, duration)
    states = kernel.prefix_states(u)
    assert np.array_equal(states[0], np.zeros(6))

    bumped = u.copy()
    bumped[0] += 1.0
    assert np.array_equal(kernel.prefix_states(bumped)[0], np.zeros(6))
    assert not np.allclose(kernel.prefix_states(bumped)[1], states[1])


# ---------------------------------------------------------------------------
# Cross terms and gradient
# ---------------------------------------------------------------------------


def test_cross_terms_match_the_dense_product(instance) -> None:
    edge_rows, weights, duration, u, _ = instance
    kernel = CouplingKernel.from_arc(edge_rows, weights, duration)

    dense = np.stack([
        sum(kernel.suffix[max(i, k)] @ u[k] for k in range(M_CELLS)) / duration
        for i in range(M_CELLS)
    ])
    assert np.allclose(kernel.cross_terms(u), dense, rtol=1e-11, atol=0.0)


def test_gradient_matches_a_finite_difference(instance) -> None:
    """``grad J = 2 Q u``, verified as the derivative it claims to be."""
    edge_rows, weights, duration, u, _ = instance
    kernel = CouplingKernel.from_arc(edge_rows, weights, duration)

    analytic = kernel.gradient(u)
    step = 1e-6
    numeric = np.empty_like(u)
    for i in range(M_CELLS):
        for a in range(6):
            plus, minus = u.copy(), u.copy()
            plus[i, a] += step
            minus[i, a] -= step
            numeric[i, a] = ((kernel.objective(plus) - kernel.objective(minus))
                             / (2.0 * step))
    assert np.allclose(analytic, numeric, rtol=1e-6,
                       atol=1e-8 * np.abs(numeric).max())


def test_objective_is_the_quadratic_form_of_its_own_gradient(instance) -> None:
    """Euler's identity for a quadratic form: ``u . grad J = 2 J``."""
    edge_rows, weights, duration, u, _ = instance
    kernel = CouplingKernel.from_arc(edge_rows, weights, duration)
    assert float(np.sum(u * kernel.gradient(u))) == pytest.approx(
        2.0 * kernel.objective(u), rel=1e-11)


# ---------------------------------------------------------------------------
# Positive semidefiniteness
# ---------------------------------------------------------------------------


def test_objective_is_non_negative_for_random_inputs(instance) -> None:
    """``Q >= 0`` is what the convex relaxation rests on."""
    edge_rows, weights, duration, _, _ = instance
    kernel = CouplingKernel.from_arc(edge_rows, weights, duration)
    rng = np.random.default_rng(3)
    for _ in range(50):
        assert kernel.objective(rng.normal(size=(M_CELLS, 6))) >= 0.0


def test_suffix_blocks_are_symmetric_and_decreasing(instance) -> None:
    """Each ``A_i`` is a sum of fewer PSD terms than ``A_{i-1}``."""
    edge_rows, weights, duration, _, _ = instance
    kernel = CouplingKernel.from_arc(edge_rows, weights, duration)
    for i in range(M_CELLS):
        block = kernel.suffix[i]
        assert np.allclose(block, block.T, rtol=1e-12, atol=0.0)
        assert np.min(np.linalg.eigvalsh(block)) > -1e-9
        if i:
            drop = kernel.suffix[i - 1] - block
            assert np.min(np.linalg.eigvalsh(drop)) > -1e-9


def test_last_cell_sees_only_the_final_edge(instance) -> None:
    """``A_{M-1}`` is a single term, which pins the top of the recursion."""
    edge_rows, weights, duration, _, _ = instance
    kernel = CouplingKernel.from_arc(edge_rows, weights, duration)
    last = weights[-1] * edge_rows[-1].T @ edge_rows[-1]
    assert np.allclose(kernel.suffix[-1], last, rtol=1e-12, atol=0.0)


# ---------------------------------------------------------------------------
# Local sensitivity kernel
# ---------------------------------------------------------------------------


def test_local_kernel_is_the_diagonal_of_the_objective_kernel(instance) -> None:
    """``K_i`` must equal the direct sum over the edges the cell is felt at."""
    edge_rows, weights, duration, _, phi = instance
    kernel = CouplingKernel.from_arc(edge_rows, weights, duration)

    rng = np.random.default_rng(5)
    transport = rng.normal(size=(M_CELLS, 6, 3))
    got = kernel.local_kernel(transport)

    for i in (0, 3, M_CELLS - 1):
        direct = sum(
            weights[j] * (edge_rows[j] @ transport[i]).T
            @ (edge_rows[j] @ transport[i])
            for j in range(i + 1, M_CELLS + 1)
        ) / duration
        assert np.allclose(got[i], direct, rtol=1e-11, atol=0.0), i


def test_local_kernel_is_symmetric_psd(instance) -> None:
    edge_rows, weights, duration, _, _ = instance
    kernel = CouplingKernel.from_arc(edge_rows, weights, duration)
    transport = np.random.default_rng(6).normal(size=(M_CELLS, 6, 3))
    for block in kernel.local_kernel(transport):
        assert np.allclose(block, block.T, rtol=1e-12, atol=0.0)
        assert np.min(np.linalg.eigvalsh(block)) > -1e-9


# ---------------------------------------------------------------------------
# Terminal special case
# ---------------------------------------------------------------------------


def test_terminal_objective_is_the_collapsed_quadrature(instance) -> None:
    """With all weight on the last edge, ``J`` is the terminal displacement.

    And every suffix block collapses to the same matrix, which is the cheapest
    available check that the recursion is aligned.
    """
    edge_rows, _, _, u, _ = instance
    duration = 3.0
    weights = np.zeros(M_CELLS + 1)
    weights[-1] = duration
    kernel = CouplingKernel.from_arc(edge_rows, weights, duration)

    expected = edge_rows[-1] @ u.sum(axis=0)
    assert kernel.objective(u) == pytest.approx(float(expected @ expected),
                                                rel=1e-12)
    for i in range(M_CELLS):
        assert np.allclose(kernel.suffix[i], kernel.suffix[0], rtol=1e-12,
                           atol=0.0)


# ---------------------------------------------------------------------------
# Summation
# ---------------------------------------------------------------------------


def test_compensated_sum_survives_cancellation() -> None:
    """The reason the prefix scan is not ``np.cumsum``."""
    x = np.array([1.0, 1e16, 1.0, -1e16])
    assert float(compensated_prefix_sum(x)[-1]) == 2.0
    assert float(np.cumsum(x)[-1]) == 0.0


def test_compensated_sum_matches_cumsum_when_benign() -> None:
    rng = np.random.default_rng(9)
    x = rng.uniform(0.5, 1.5, (200, 6))
    assert np.allclose(compensated_prefix_sum(x), np.cumsum(x, axis=0),
                       rtol=1e-14, atol=0.0)


def test_conditioning_reports_the_cancellation(instance) -> None:
    """Large is expected: it is the cancellation the method exploits."""
    edge_rows, weights, duration, _, _ = instance
    kernel = CouplingKernel.from_arc(edge_rows, weights, duration)

    aligned = np.tile(np.array([1.0, 0, 0, 0, 0, 0]), (M_CELLS, 1))
    assert kernel.conditioning(aligned) == pytest.approx(1.0)

    cancelling = aligned.copy()
    cancelling[1::2] *= -1.0
    assert kernel.conditioning(cancelling) > 5.0


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(("kwargs", "message"), [
    ({"edge_rows": np.zeros((5, 6, 6))}, r"shape \(M\+1, 3, 6\)"),
    ({"outer_weights": np.zeros(3)}, "must match"),
    ({"outer_weights": -np.ones(M_CELLS + 1)}, "non-negative"),
    ({"duration": 0.0}, "duration must be positive"),
])
def test_construction_rejects_bad_input(instance, kwargs, message) -> None:
    edge_rows, weights, duration, _, _ = instance
    args = {"edge_rows": edge_rows, "outer_weights": weights,
            "duration": duration}
    args.update(kwargs)
    with pytest.raises(ValueError, match=message):
        CouplingKernel.from_arc(**args)


def test_wrong_shaped_u_is_rejected(instance) -> None:
    edge_rows, weights, duration, _, _ = instance
    kernel = CouplingKernel.from_arc(edge_rows, weights, duration)
    with pytest.raises(ValueError, match=r"shape \(7, 6\)"):
        kernel.objective(np.zeros((M_CELLS, 3)))
