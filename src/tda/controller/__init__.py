"""The deployable controllers: what is left when the reference field is gone.

The benchmark of :mod:`tda.allocate` evaluates the reference field at every
epoch of a reference trajectory it already knows.  This package removes that
access and measures what each removal costs --- RQ3.  The distinction is not
between having and not having the gravity model, which a propagator carries by
definition, but between *global, offline model metadata* and
*trajectory-specific evaluations of it*: the controller may know the degree
variances, and may not be told what the omitted acceleration is at the epochs
it will fly through.

Four modules, and the split between them is the architecture
------------------------------------------------------------
``phase``
    Indexing by ``revolution + phase`` rather than by absolute time.  A loose
    pilot arc drifts along track by a timing offset of order a perilune
    passage, and the plan's features *are* the perilune passages.
``plan``
    Algorithm 2, offline, on a pilot arc: the cancellation target, the
    multiplier and the nominal degrees.  It reuses the benchmark's solver
    unchanged, so the plan and the benchmark differ in information and in
    nothing else.
``online``
    Algorithm 3, one decision per boundary: predict, probe ahead, score the
    window, hold.  ``C-plan`` and ``C-lite`` are two settings of one code path.
``feedback``
    ``C-fb``: the multiplier closed on *realized* work, in the accounting the
    integrator actually spends in.

Which half computes what is forced by a measurement, not by taste.  The local
omitted direction is field texture and decorrelates over :math:`\\pi r/N`, so a
value computed on the pilot arc is worthless by the time it is flown --- it has
to be probed online.  The cancellation target is an *accumulated* displacement,
dominated by a sustained trend rather than by local texture, so it survives the
transfer and can be planned.  Neither half reads the reference field or the
reference arc.

What this package does not claim
--------------------------------
That this is the best deployable controller.  It is the first interpretable,
auditable, low-complexity attempt, and the campaign's branching stage measures
which controller family the problem actually needs before the expensive
propagations are spent (``DECISIONS.md`` D116).  The claim that a
trajectory-aware allocation is worth having is separate, rests on the benchmark,
and does not depend on this package being the right answer.
"""

from tda.controller.feedback import (
    BudgetFeedback,
    ReferenceProfile,
    WorkTracker,
    overspend_ratio,
)
from tda.controller.online import (
    Decision,
    OnlineController,
    OnlineSettings,
    assign_probe_points,
    colocated_cost_fraction,
    probe_fractions,
    score_candidates,
    select_candidate,
)
from tda.controller.phase import (
    PhaseIndexError,
    RevolutionIndex,
    apsis_epochs,
    osculating_period,
    radial_rate,
)
from tda.controller.plan import OfflinePlan, build_plan, cancellation_target

__all__ = [
    "BudgetFeedback",
    "Decision",
    "OfflinePlan",
    "OnlineController",
    "OnlineSettings",
    "PhaseIndexError",
    "ReferenceProfile",
    "RevolutionIndex",
    "WorkTracker",
    "apsis_epochs",
    "assign_probe_points",
    "build_plan",
    "cancellation_target",
    "colocated_cost_fraction",
    "osculating_period",
    "overspend_ratio",
    "probe_fractions",
    "radial_rate",
    "score_candidates",
    "select_candidate",
]
