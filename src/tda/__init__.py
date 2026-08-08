"""Trajectory-aware truncation-degree allocation.

The package implements the allocation problem of the manuscript: choose a
piecewise-constant spherical-harmonic truncation degree :math:`N(t)` over an
arc so as to minimise the induced position error under a quadratic work
budget.  What distinguishes it from force-level degree selection is that the
objective takes its norm *after* a signed integral against the state-transition
matrix, so the epochs are coupled and the allocation is not separable.

Layering
--------
The modules form a strict dependency chain; nothing below reaches upward.

``config``
    Frozen run constants and their provenance hash.
``field``
    Adapter over the gravity kernel: acceleration, truncation defect,
    gravity gradient, and the omitted-band stack the online probe needs.
``dynamics``
    Reference trajectory and its state-transition matrix, from the
    variational equations.
``grids`` / ``tables``
    Accumulation and decision grids, and the tabulated defect vectors.
``kernel``
    The objective's coupling structure: the suffix sequence, the local
    sensitivity kernel, the cross terms and the objective value.
``allocate``
    Solvers: separable, coordinate descent, Frank--Wolfe certificate,
    relax-and-round.
``policies`` / ``probe`` / ``controller``
    Comparators, the band probe, and the deployable controllers.
``analysis``
    The architecture-selection tests T1--T7 of WP21.

Conventions
-----------
* SI throughout: metres, seconds, radians.  Degrees of the expansion are the
  only dimensionless integers.
* Vectors are ``numpy`` arrays with the component axis **last**, so an arc of
  positions has shape ``(M, 3)`` and an arc of state-transition matrices has
  shape ``(M, 6, 6)``.
* Every quantity that enters a reported number carries the run configuration
  hash of :func:`tda.config.RunConfig.digest`.

See ``../NOTATION.md`` for the symbol table this code follows and
``../DECISIONS.md`` for why each convention is what it is.
"""

from tda.config import RunConfig

__all__ = ["RunConfig"]
__version__ = "0.1.0"
