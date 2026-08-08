"""Frozen run configuration and its provenance digest.

Every constant that can change a reported number lives here, in one immutable
object, so that a result can be traced to the settings that produced it.  The
alternative -- module-level constants scattered across the code -- is what made
the previous campaign's manifest audit expensive, and it is the reason
``DECISIONS.md`` carries a manifest digest chain as an inherited convention.

Three properties are enforced rather than documented:

* the objects are frozen dataclasses, so nothing can mutate a setting after a
  run has started;
* :meth:`RunConfig.digest` hashes the *whole* configuration, so a table can
  record which settings produced it and a later reader can detect a mismatch;
* values inherited from the previous campaign are marked as such, because
  changing one of them breaks comparability with the archive rather than
  merely changing a number here.

References
----------
.. [Hairer1993] E. Hairer, S. P. Nørsett, G. Wanner, *Solving Ordinary
   Differential Equations I*, 2nd ed., Springer, 1993 -- integrator families
   and tolerance conventions.
.. [Press2007] W. H. Press et al., *Numerical Recipes*, 3rd ed., CUP, 2007,
   §5.7 -- optimal step size for numerical derivatives.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import asdict, dataclass, field
from typing import Literal

__all__ = [
    "BudgetConfig",
    "GradientConfig",
    "GridConfig",
    "IntegratorConfig",
    "RunConfig",
]

# --------------------------------------------------------------------------
# Inherited from the previous campaign.  Changing any of these breaks
# comparability with the read-only archive and requires a decision entry.
# --------------------------------------------------------------------------

ARC_DURATION_S: float = 7.0 * 86400.0
"""Seven-day arc.  Inherited; the archive's populations are calibrated on it."""

OUTPUT_STEP_S: float = 120.0
"""Reference output cadence of the archive, used for admissibility checks."""

SIDEREAL_MONTH_S: float = 27.321661 * 86400.0
"""Sidereal rotation period of the Moon, as the archive defines it."""

OMEGA_MOON_RAD_S: float = 2.0 * math.pi / SIDEREAL_MONTH_S
"""Uniform lunar rotation rate used by the archive's inertial/body transform.

A *simplification the archive adopted* -- a uniform rotation about the
inertial z-axis rather than the full librating orientation -- reproduced here
verbatim rather than improved, because the campaign compares against archived
numbers.

Written as the same expression the archive uses
(``rev3_common.OMEGA_MOON``) rather than as a rounded literal, so that the
admissibility check of WP0 compares two identical floats instead of two values
that agree to seven digits.  A rounded literal here would put a systematic
along-track phase error into every body-fixed evaluation, which is exactly the
class of error a signed method cannot tolerate and a magnitude statistic would
not notice.
"""


@dataclass(frozen=True, slots=True)
class IntegratorConfig:
    """Settings for the reference and variational integrations.

    The defaults reproduce the archive's forced-variational runs so that the
    admissibility check of WP0 compares like with like.

    Attributes
    ----------
    method:
        SciPy integrator name.  ``DOP853`` is the eighth-order explicit
        Runge--Kutta pair of Dormand and Prince [Hairer1993]_, which the
        previous campaign selected and whose behaviour under a discontinuous
        right-hand side (a degree switch) it characterised.
    rtol, atol:
        Relative and *vector* absolute tolerances.  ``atol`` is a per-component
        sequence rather than a scalar: position and velocity differ by six
        orders of magnitude in SI, and a scalar ``atol`` silently sets the
        velocity tolerance far too loose.  The previous campaign traced a
        noise-floor artefact to exactly that.
    max_step_s:
        Cap on the step size.  Present so that a degree switch cannot be
        stepped over.

        **This differs from the archive**, whose ``rev3_common.propagate``
        defaults to ``inf``.  A bit-exact reproduction of an archived value
        (WP0) must therefore pass the archive's own value; the default here is
        the campaign's, not the archive's, and the two are deliberately
        distinct rather than silently reconciled.
    """

    method: Literal["DOP853", "RK45", "Radau"] = "DOP853"
    rtol: float = 1.0e-11
    atol_position_m: float = 1.0e-9
    atol_velocity_m_s: float = 1.0e-12
    max_step_s: float = 60.0

    def atol_vector(self, n_extra: int = 0,
                    atol_extra: float = 1.0e-12) -> list[float]:
        """Return the per-component absolute tolerance for an augmented state.

        Parameters
        ----------
        n_extra:
            Number of trailing states beyond the six-dimensional Cartesian
            state -- 36 when the state-transition matrix is carried alongside.
        atol_extra:
            Tolerance for those trailing states.  The state-transition matrix
            is dimensionless in its position--position block but carries units
            elsewhere; a single tight value is the conservative choice and its
            cost is small because the same steps serve the physical state.

        Returns
        -------
        list of float
            Length ``6 + n_extra``.
        """
        base = [self.atol_position_m] * 3 + [self.atol_velocity_m_s] * 3
        return base + [atol_extra] * n_extra


@dataclass(frozen=True, slots=True)
class GradientConfig:
    """How the gravity gradient that generates the STM is obtained.

    The gradient :math:`\\mathbf G=\\partial\\mathbf a/\\partial\\mathbf r`
    is the only input to the variational equations, and its truncation degree
    is an open question of the plan (``DECISIONS.md`` Q13): the previous
    campaign measured that a degree-120 gradient is *not* adequate on
    31 km-perilune orbits, and this campaign carries a low-perilune
    population.  The degree is therefore an explicit, recorded parameter with
    no silent default, and ``None`` means "match the orbit's reference
    degree", which is the conservative choice.

    Attributes
    ----------
    degree:
        Truncation degree for the gradient, or ``None`` to match the
        trajectory's reference degree.
    step_m:
        Central-difference step.  The optimum for a central difference is
        near :math:`\\epsilon^{1/3}` times the scale over which the third
        derivative varies [Press2007]_; at lunar-orbit radii that is metres,
        and the archive's value of 1 m is retained.  A convergence sweep is
        part of WP1 rather than an assumption.
    symmetrise:
        Replace :math:`\\mathbf G` by :math:`(\\mathbf G+\\mathbf G^\\top)/2`.
        The exact Hessian of a potential is symmetric, so any asymmetry in a
        finite-difference estimate is noise.  Enforcing the symmetry is not
        cosmetic: it is what makes the flow's state-transition matrix
        symplectic, and hence what makes the symplectic defect of
        :func:`tda.dynamics.symplectic_defect` a diagnostic of the
        *integration* rather than a readout of differencing noise.
    """

    degree: int | None = None
    step_m: float = 1.0
    symmetrise: bool = True


@dataclass(frozen=True, slots=True)
class GridConfig:
    """Accumulation and decision grids.

    Two grids, never conflated (manuscript §3.5).  The accumulation grid
    carries the integrand and is refined on the correlation time
    :math:`\\tau_{\\mathrm{corr}}=\\pi r/(Nv)`; the decision grid carries the
    piecewise-constant degree and is coarse enough to be flyable.

    Attributes
    ----------
    samples_per_tau:
        :math:`n_s=\\tau_{\\mathrm{corr}}/\\Delta t_i`, the convergence
        parameter of the adaptive accumulation grid.  Swept in WP4; the value
        here is the campaign default, declared in advance.
    dt_acc_min_s, dt_acc_max_s:
        Clamps on the refined step, so that a near-zero :math:`\\tau` at a
        very low perilune cannot generate an unbounded number of samples.
    dt_dec_s:
        Decision interval.  Every policy compared -- benchmark and controller
        alike -- switches on this grid, so that the benchmark's margin is
        information rather than switching resolution.
    """

    samples_per_tau: float = 2.0
    dt_acc_min_s: float = 2.0
    dt_acc_max_s: float = 240.0
    dt_dec_s: float = 120.0


@dataclass(frozen=True, slots=True)
class BudgetConfig:
    """Work accounting.

    Four cost quantities that are never interchanged (manuscript §5.3).  Only
    the conventions live here; the measured values are run outputs.

    Attributes
    ----------
    beta:
        Allocation target :math:`\\beta=B_1/N_{\\mathrm{crit}}^2`, defined at
        the nominal level only.
    match_tolerance:
        Relative tolerance on :math:`|B_2+B_+-B_{\\mathrm{tot}}|
        /B_{\\mathrm{tot}}` for a candidate to count as calibrated to the
        anchor.
    max_calibration_passes:
        Cap on the propagate--measure--rescale loop, so a candidate that will
        not converge cannot consume the campaign.
    """

    beta: float = 1.0
    match_tolerance: float = 0.02
    max_calibration_passes: int = 4


@dataclass(frozen=True, slots=True)
class RunConfig:
    """The complete, immutable configuration of one run.

    Examples
    --------
    >>> cfg = RunConfig(label="wp1-smoke")
    >>> cfg.integrator.method
    'DOP853'
    >>> len(cfg.digest()) == 16
    True
    >>> cfg.integrator.rtol = 1e-9              # doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
    dataclasses.FrozenInstanceError: cannot assign to field 'rtol'
    """

    label: str
    integrator: IntegratorConfig = field(default_factory=IntegratorConfig)
    gradient: GradientConfig = field(default_factory=GradientConfig)
    grids: GridConfig = field(default_factory=GridConfig)
    budget: BudgetConfig = field(default_factory=BudgetConfig)
    arc_duration_s: float = ARC_DURATION_S

    def as_dict(self) -> dict:
        """Return a plain, JSON-serialisable view of the configuration."""
        return asdict(self)

    def digest(self, length: int = 16) -> str:
        """Return a stable short hash of the configuration.

        The hash covers every field, is independent of dictionary ordering,
        and is what a metrics file records so that a later reader can tell
        whether two tables were produced under the same settings.

        Parameters
        ----------
        length:
            Number of hex characters to keep.  Sixteen gives a 64-bit
            fingerprint, ample for distinguishing the configurations of one
            campaign and short enough to sit in a filename.
        """
        payload = json.dumps(self.as_dict(), sort_keys=True,
                             separators=(",", ":"))
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:length]
