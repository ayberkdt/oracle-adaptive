"""The certificate: a convex-hull relaxation and a Frank--Wolfe lower bound.

Calling the descent result a benchmark requires knowing how far it can be from
the best possible schedule.  Relax each decision interval to the convex hull
of its candidate contributions,

.. math::
    \\mathbf u_i=\\sum_{N}\\theta_{g(i)N}\\,\\mathbf u_i(N),\\qquad
    \\theta_{qN}\\ge0,\\quad\\sum_N\\theta_{qN}=1,\\quad
    \\sum_qW_q\\sum_N\\theta_{qN}N^2\\le B .

Every integer schedule is feasible here with :math:`\\theta` an indicator and
attains the same objective, so the relaxed optimum is a lower bound on the
integer one, and Frank--Wolfe's duality gap supplies it directly:

.. math::
    J(\\mathbf u)+\\langle\\nabla J,\\mathbf s-\\mathbf u\\rangle
    \\;\\le\\;J^{*}_{\\mathrm{rel}}\\;\\le\\;J^{*}\\;\\le\\;
    J(\\mathbf u_{\\mathrm{descent}}) .

The subproblem is solved as the linear programme it is
------------------------------------------------------
The bound is valid only if :math:`\\mathbf s` attains the minimum of the
linearisation over the **relaxed** feasible set --- and over that set, not
over its boundary.  Two ways to get this wrong, both tempting:

*Forcing the budget to bind.*  The objective is not monotone in degree, and
neither is :math:`\\langle\\nabla_qJ,\\mathbf u_q(N)\\rangle`, so the
programme's optimum may leave slack.  Pushing it onto the budget surface
returns a feasible point rather than a minimiser and the bound stops being one
(``DECISIONS.md`` D142).

*Mixing the wrong pair.*  At a critical multiplier a fractional solution mixes
neighbouring vertices of the **lower convex envelope** of the points
:math:`(N^2,\\langle\\nabla_qJ,\\mathbf u_q(N)\\rangle)`, which need not be
neighbours in the candidate set: a candidate above the envelope is never
selected at any multiplier, so with candidates :math:`\\{40,50,60\\}` the exact
solution may mix 40 with 60.  Several intervals can also tie at the same
multiplier (D132).

Rather than implement that reasoning, the subproblem is handed to
``scipy.optimize.linprog``.  It is exact by construction rather than by an
argument about envelopes, and the certificate is one of the paper's stronger
claims: a dull verified solver is worth more here than a clever fragile one.

Why the away step is not an optimisation
----------------------------------------
Classical Frank--Wolfe converges as :math:`\\mathcal O(1/t)` and stalls badly
when the optimum lies on a face of the polytope: the iterate zig-zags between
vertices it cannot leave, because the only direction available points *toward*
a new vertex.  On this objective that is the normal case rather than a corner
one --- the whole content of the signed formulation is near-cancellation, so
the quadratic is severely ill-conditioned and the relaxed optimum sits on a
low-dimensional face.  A pilot run at four hundred iterations returned a bound
that never rose above zero, which said nothing about the schedule and
everything about the step.

The away step [GuelatMarcotte1986]_ adds the missing direction: move *away*
from the worst vertex currently carrying weight, which lets the iterate shed
an atom instead of merely diluting it.  On a polytope with a strongly convex
objective this restores a linear rate [LacosteJulien2015]_.  The bound itself
is unchanged --- it is still the linearisation at the current iterate --- so
away steps buy convergence, not validity, and a bound obtained with them is
the same kind of object as one obtained without.

Distinguishing a slow certificate from a useless one
----------------------------------------------------
A vacuous bound has two possible causes and they call for opposite responses.
Either the relaxation has not been solved yet, or it has been solved and its
optimum really is near zero --- which happens when fractional mixing can
cancel the accumulated displacement far better than any integer schedule, and
means no amount of further computation will produce a certificate.
:attr:`Certificate.relaxed_gap` separates them: it is the Frank--Wolfe duality
gap at termination, which bounds how far the *relaxation's* own iterate still
is from the relaxed optimum.  Small gap with a zero bound means the relaxation
is loose and the certificate has to be built differently;
:attr:`Certificate.structurally_vacuous` reports that case by name instead of
letting it read as "run it longer".

References
----------
.. [Frank1956] M. Frank and P. Wolfe, "An algorithm for quadratic
   programming", *Naval Research Logistics Quarterly* 3, 1956.
.. [Jaggi2013] M. Jaggi, "Revisiting Frank--Wolfe: projection-free sparse
   convex optimization", *ICML*, 2013 -- the duality gap as a certificate.
.. [GuelatMarcotte1986] J. Guelat and P. Marcotte, "Some comments on
   Wolfe's away step", *Mathematical Programming* 35, 1986.
.. [LacosteJulien2015] S. Lacoste-Julien and M. Jaggi, "On the global linear
   convergence of Frank--Wolfe optimization variants", *NeurIPS*, 2015 --
   the away-step and pairwise variants and their linear rates.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy import sparse
from scipy.optimize import linprog

from tda.allocate.descent import DescentProblem

__all__ = ["Certificate", "certify", "linear_minimisation"]

Arr = NDArray[np.float64]

_WEIGHT_FLOOR = 1.0e-12
"""Below this an active vertex is treated as gone rather than as tiny."""


@dataclass(frozen=True, slots=True)
class Certificate:
    """The certified distance between the descent schedule and the optimum.

    Attributes
    ----------
    lower_bound:
        :math:`L_{\\mathrm{FW}}`, the best (largest) bound the iteration
        attained, clipped at zero since :math:`J\\ge0` always.
    descent_objective:
        :math:`J_{\\mathrm{desc}}`.
    gap_objective:
        :math:`g_J=(J_{\\mathrm{desc}}-L_{\\mathrm{FW}})/J_{\\mathrm{desc}}`.
    gap_error:
        :math:`g_E=1-\\sqrt{L_{\\mathrm{FW}}/J_{\\mathrm{desc}}}`.  The
        quantity the paper's claims are in: the descent schedule's position
        error is within this of the best any schedule of the class could
        achieve.  A ten per cent error gap is a nineteen per cent objective
        gap, which is why the two are never quoted for one another.
    vacuous:
        Whether the bound stayed at zero for the whole iteration.  Reported in
        its own column rather than dropped, because an orbit with no
        certificate is a different statement from an orbit with a wide one.
    iterations:
        Frank--Wolfe steps taken.
    relaxed_objective:
        :math:`J` at the relaxation's final iterate.  Not a bound on anything;
        carried so that :attr:`relaxed_gap` can be read against a scale.
    relaxed_gap:
        The Frank--Wolfe duality gap at termination,
        :math:`\\langle\\nabla J,\\mathbf u-\\mathbf s\\rangle`, which bounds
        :math:`J(\\mathbf u)-J^{*}_{\\mathrm{rel}}`.  This is a statement about
        how well the *relaxation* was solved, not about the schedule.
    away_steps, drop_steps:
        How many of the steps moved away from an active vertex, and how many
        of those removed it entirely.
    active_atoms:
        Size of the active vertex set at termination.
    """

    lower_bound: float
    descent_objective: float
    gap_objective: float
    gap_error: float
    vacuous: bool
    iterations: int
    relaxed_objective: float = 0.0
    relaxed_gap: float = float("inf")
    away_steps: int = 0
    drop_steps: int = 0
    active_atoms: int = 0

    @property
    def structurally_vacuous(self) -> bool:
        """Whether the relaxation, not the solver, is what produced no bound.

        True when the bound is zero *and* the relaxation's own remaining gap
        is small next to the quantity being certified.  In that case the
        relaxed optimum really is near zero: fractional mixing cancels the
        accumulated displacement in a way no integer schedule can, and further
        iterations cannot help.  The response is a different certificate ---
        a tighter relaxation, or a bound that does not go through one --- and
        not a longer run.
        """
        return (self.vacuous and self.descent_objective > 0.0
                and self.relaxed_gap < 0.01 * self.descent_objective)

    def earns_the_name(self, threshold: float = 0.10) -> bool:
        """Whether the gap clears the naming rule's threshold.

        The rule is fixed before the runs and bound to :attr:`gap_error`, not
        to the result: below the threshold the quantity may be called an
        oracle, otherwise it is reported as a linearised trajectory-aware
        allocation benchmark and every ratio against it read as conservative.
        """
        return (not self.vacuous) and self.gap_error < threshold


def _feasible_set(problem: DescentProblem, budget: float):
    """Build the relaxation's constraint matrices.

    Returns the simplex equalities, the single knapsack row and the cost
    vector, all in the flattened ``(q, N)`` variable order.  Sparse because
    the equality block is :math:`K\\times KP` with exactly :math:`KP`
    non-zeros -- dense it would be a hundred gigabytes at campaign scale.
    """
    n_intervals = problem.n_intervals
    degrees = np.asarray(problem.candidate_degrees, dtype=float)
    n_candidates = degrees.size
    n_vars = n_intervals * n_candidates

    rows = np.repeat(np.arange(n_intervals), n_candidates)
    cols = np.arange(n_vars)
    a_eq = sparse.csr_matrix(
        (np.ones(n_vars), (rows, cols)), shape=(n_intervals, n_vars))
    b_eq = np.ones(n_intervals)

    work = (problem.time_weight[:, None] * degrees[None, :] ** 2).ravel()
    a_ub = sparse.csr_matrix(work.reshape(1, -1))
    return a_eq, b_eq, a_ub, np.array([budget]), work


def objective_coefficients(problem: DescentProblem, gradient: Arr) -> Arr:
    """:math:`\\langle\\nabla_qJ,\\mathbf u_q(N)\\rangle`, shape ``(K, P)``.

    The linearised cost of putting interval :math:`q` entirely on candidate
    :math:`N`.  Exposed rather than buried in the subproblem because the
    away-step search needs the same numbers: scoring a stored vertex is then
    :math:`\\langle\\mathbf c,\\theta\\rangle`, an inner product over ``K*P``
    entries, instead of a contraction over all :math:`M` cells.  With a few
    hundred active vertices that difference is what keeps the away step from
    costing more than the iteration it accelerates.
    """
    per_cell = np.einsum("ipa,ia->ip", problem.contributions, gradient)
    out = np.zeros((problem.n_intervals, per_cell.shape[1]))
    np.add.at(out, problem.interval_of, per_cell)
    return out


def linear_minimisation(problem: DescentProblem, gradient: Arr,
                        budget: float) -> Arr:
    """Solve the Frank--Wolfe subproblem exactly, as a linear programme.

    Parameters
    ----------
    problem:
        The relaxation's data.
    gradient:
        :math:`\\nabla J` at the current iterate, shape ``(M, 6)``.
    budget:
        The work ceiling.

    Returns
    -------
    ndarray, shape (K, P)
        The vertex :math:`\\theta` attaining the minimum.

    Raises
    ------
    RuntimeError
        If the solver does not certify optimality.  A run in which it does not
        produces no bound: the whole value of the certificate is that the
        subproblem was solved rather than approximated.
    """
    a_eq, b_eq, a_ub, b_ub, _ = _feasible_set(problem, budget)
    coefficients = objective_coefficients(problem, gradient)

    result = linprog(coefficients.ravel(), A_ub=a_ub, b_ub=b_ub,
                     A_eq=a_eq, b_eq=b_eq, bounds=(0.0, 1.0),
                     method="highs")
    if not result.success:
        raise RuntimeError(
            f"the linear minimisation did not solve: {result.message}. "
            "No bound is reported for this orbit rather than a bound from an "
            "unsolved subproblem.")
    return result.x.reshape(coefficients.shape)


def _contract(problem: DescentProblem, theta: Arr) -> Arr:
    """:math:`\\mathbf u(\\theta)`, shape ``(M, 6)``."""
    return np.einsum("ipa,ip->ia", problem.contributions,
                     theta[problem.interval_of])


def certify(problem: DescentProblem, descent_objective: float, budget: float,
            iterations: int = 60, away_steps: bool = True,
            tolerance: float = 1.0e-9) -> Certificate:
    """Run Frank--Wolfe on the relaxation and return the certified gap.

    Parameters
    ----------
    problem:
        The relaxation's data.
    descent_objective:
        :math:`J` of the schedule being certified.
    budget:
        The work ceiling.
    iterations:
        Step budget.  The bound is valid at every iteration; more only makes
        it tighter.
    away_steps:
        Whether to allow away steps.  ``False`` is classical Frank--Wolfe and
        is retained as the control: the two produce the same *kind* of bound,
        so comparing them measures the solver and nothing else.
    tolerance:
        Stop once the relaxation's own duality gap falls below this fraction
        of the scale being certified.  With away steps the gap reaches zero
        rather than merely shrinking, and once it does no further step can
        raise the bound; continuing would only spend the iteration budget of
        the next orbit.

    Returns
    -------
    Certificate

    Notes
    -----
    The step is the exact line search, which for a quadratic is closed form:
    :math:`\\gamma^{\\star}=-\\langle\\nabla J,\\mathbf d\\rangle/(2J(\\mathbf
    d))`, clipped to :math:`[0,\\gamma_{\\max}]`.  The usual :math:`2/(t+2)`
    schedule is a fallback for objectives whose curvature is unknown; here it
    is known exactly and guessing would only slow the bound down.

    For a Frank--Wolfe step :math:`\\gamma_{\\max}=1`; for an away step it is
    :math:`\\alpha_a/(1-\\alpha_a)`, the largest move that keeps the away
    vertex's weight non-negative.  Reaching it is a *drop step* and removes
    the vertex from the active set --- the mechanism that classical
    Frank--Wolfe lacks and the reason it stalls on a face.

    The bound can be negative early in the iteration.  Since :math:`J\\ge0`
    always, it is clipped at zero, and a bound that never rises above zero is
    reported as vacuous rather than as a gap of one hundred per cent.

    Active vertices are stored sparsely.  A vertex of this polytope assigns
    each interval one candidate except at most one interval split by the
    knapsack, so it has at most :math:`K+1` non-zeros out of :math:`KP`;
    holding them densely would cost a megabyte apiece at campaign scale.
    """
    kernel = problem.kernel
    n_intervals = problem.n_intervals
    n_candidates = len(problem.candidate_degrees)
    shape = (n_intervals, n_candidates)

    def as_atom(theta: Arr):
        """Store a vertex as a sparse row."""
        return sparse.csr_matrix(np.asarray(theta).reshape(1, -1))

    def contract_atom(atom) -> Arr:
        """:math:`\\mathbf u` of a stored vertex, shape ``(M, 6)``."""
        return _contract(problem, np.asarray(atom.todense()).reshape(shape))

    # Start from the cheapest feasible vertex: all mass on the lowest degree.
    start = np.zeros(shape)
    start[:, 0] = 1.0
    atoms = [as_atom(start)]
    weights = np.array([1.0])
    u = _contract(problem, start)

    best_bound = -np.inf
    taken = away_taken = drops = 0
    relaxed_gap = float("inf")

    for step_index in range(1, iterations + 1):
        taken = step_index
        gradient = kernel.gradient(u)
        coefficients = objective_coefficients(problem, gradient)
        s = linear_minimisation(problem, gradient, budget)
        u_s = _contract(problem, s)

        forward = u_s - u
        forward_gain = -float(np.sum(gradient * forward))
        relaxed_gap = max(forward_gain, 0.0)
        current = kernel.objective(u)
        best_bound = max(best_bound, current - forward_gain)
        if relaxed_gap <= tolerance * max(abs(descent_objective), current):
            break                      # the relaxation is solved

        direction, ceiling, is_away, worst = forward, 1.0, False, -1
        if away_steps and len(atoms) > 1:
            flat = coefficients.ravel()
            scores = np.array([float(atom.dot(flat)[0]) for atom in atoms])
            worst = int(np.argmax(scores))
            alpha = float(weights[worst])
            u_away = contract_atom(atoms[worst])
            candidate = u - u_away
            # A vertex carrying no weight cannot be moved away from: its
            # ceiling is zero, the step is zero, and the iteration reads that
            # as convergence and stops. Guarded here as well as by pruning,
            # because "the step was zero" and "there is nowhere to go" have
            # to stay different statements.
            if (alpha > _WEIGHT_FLOOR
                    and -float(np.sum(gradient * candidate)) > forward_gain):
                direction, is_away = candidate, True
                ceiling = (np.inf if alpha >= 1.0
                           else alpha / (1.0 - alpha))

        curvature = kernel.objective(direction)
        if curvature <= 0.0:
            break                      # the iterate is already a vertex
        step = -float(np.sum(gradient * direction)) / (2.0 * curvature)
        step = min(max(step, 0.0), ceiling)
        if step == 0.0:
            break

        u = u + step * direction
        if is_away:
            away_taken += 1
            weights = weights * (1.0 + step)
            weights[worst] -= step
            if weights[worst] <= _WEIGHT_FLOOR:
                drops += 1
        else:
            weights = weights * (1.0 - step)
            atoms.append(as_atom(s))
            weights = np.append(weights, step)
        weights = np.clip(weights, 0.0, None)
        weights /= weights.sum()

        # Prune emptied vertices. A full Frank--Wolfe step leaves every
        # previous atom at zero weight, and an away search that then picks one
        # of them gets a zero ceiling and terminates the iteration three steps
        # in --- which is how this looked the first time: a converged-looking
        # certificate that had barely started.
        if float(weights.min()) <= _WEIGHT_FLOOR:
            keep = np.flatnonzero(weights > _WEIGHT_FLOOR)
            atoms = [atoms[i] for i in keep]
            weights = weights[keep]
            weights /= weights.sum()

    lower = max(best_bound, 0.0)
    vacuous = lower <= 0.0
    if descent_objective <= 0.0:
        gap_objective = gap_error = 0.0
    else:
        gap_objective = (descent_objective - lower) / descent_objective
        gap_error = 1.0 - float(np.sqrt(lower / descent_objective))
    return Certificate(
        lower_bound=lower,
        descent_objective=descent_objective,
        gap_objective=gap_objective,
        gap_error=gap_error,
        vacuous=vacuous,
        iterations=taken,
        relaxed_objective=kernel.objective(u),
        relaxed_gap=relaxed_gap,
        away_steps=away_taken,
        drop_steps=drops,
        active_atoms=len(atoms),
    )
