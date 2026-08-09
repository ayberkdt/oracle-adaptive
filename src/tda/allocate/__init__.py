"""Solvers for the budgeted allocation.

Four rungs of one ladder, sharing one budget contract:

``separable``
    ``A-force`` and ``A-sens``.  Both decouple under a single multiplier, so
    both are solved exactly; they are the denominator of RQ2.
``descent``
    ``A-sign``.  The only solver that sees the off-diagonal coupling, and the
    only one that is not exact --- it returns an attainable schedule.
``exhaustive``
    How far that schedule is from the optimum, by enumerating every feasible
    schedule on instances small enough to allow it.  This is the primary
    verification: it compares against truth rather than against a relaxation,
    and on the pilot the descent was exactly optimal in every enumerable
    instance.
``frankwolfe``
    The secondary diagnostic, for the full-size problem where enumeration is
    impossible: a lower bound on a convex-hull relaxation whose linear
    subproblem is solved as a linear programme.  Demoted from the primary role
    once the relaxation's own integrality gap was measured and found not to
    shrink with the interval count (``DECISIONS.md`` D178).
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
from tda.allocate.exhaustive import (
    VerificationRecord,
    panel_summary,
    schedule_count,
    solve_exhaustive,
    verify_schedule,
)
from tda.allocate.frankwolfe import (
    Certificate,
    certify,
    linear_minimisation,
    objective_coefficients,
)
from tda.allocate.polish import PolishReport, block_deltas, polish_to_budget
from tda.allocate.rounding import (
    argmax_rounding,
    round_and_polish,
    sum_up_rounding,
)
from tda.allocate.separable import (
    force_values,
    greedy_fill,
    sensitivity_values,
    solve_separable,
)

__all__ = [
    "Certificate",
    "DescentProblem",
    "PolishReport",
    "ScheduleSolution",
    "VerificationRecord",
    "argmax_rounding",
    "block_deltas",
    "certify",
    "expand_to_cells",
    "force_values",
    "greedy_fill",
    "linear_minimisation",
    "monotonicity_report",
    "objective_coefficients",
    "panel_summary",
    "polish_to_budget",
    "round_and_polish",
    "schedule_count",
    "schedule_work",
    "sensitivity_values",
    "solve_descent",
    "solve_exhaustive",
    "solve_separable",
    "solve_to_budget",
    "sum_up_rounding",
    "verify_schedule",
]
