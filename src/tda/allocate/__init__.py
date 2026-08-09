"""Solvers for the budgeted allocation.

Four rungs of one ladder, sharing one budget contract:

``separable``
    ``A-force`` and ``A-sens``.  Both decouple under a single multiplier, so
    both are solved exactly; they are the denominator of RQ2.
``descent``
    ``A-sign``.  The only solver that sees the off-diagonal coupling, and the
    only one that is not exact --- it returns an attainable schedule.
``frankwolfe``
    How far that schedule can be from the optimum, certified on a convex-hull
    relaxation whose linear subproblem is solved as a linear programme.
``rounding``
    ``S-round``, a second solver out of the same relaxation.  If it beats the
    descent it becomes the reported schedule and the descent the control.

The shared contract is in ``budget``: the constraint is a **ceiling**, the
multiplier is tested at zero before anything is bracketed, and selection among
starts is on :math:`J+\\lambda W` rather than on :math:`J`.  All three are
places where the obvious implementation is wrong for this objective, because
:math:`J` is not monotone in degree.
"""

from tda.allocate.budget import (
    ScheduleSolution,
    expand_to_cells,
    schedule_work,
    solve_to_budget,
)
from tda.allocate.descent import (
    DescentProblem,
    monotonicity_report,
    solve_descent,
)
from tda.allocate.frankwolfe import Certificate, certify, linear_minimisation
from tda.allocate.rounding import (
    argmax_rounding,
    round_and_polish,
    sum_up_rounding,
)
from tda.allocate.separable import (
    force_values,
    sensitivity_values,
    solve_separable,
)

__all__ = [
    "Certificate",
    "DescentProblem",
    "ScheduleSolution",
    "argmax_rounding",
    "certify",
    "expand_to_cells",
    "force_values",
    "linear_minimisation",
    "monotonicity_report",
    "round_and_polish",
    "schedule_work",
    "sensitivity_values",
    "solve_descent",
    "solve_separable",
    "solve_to_budget",
    "sum_up_rounding",
]
