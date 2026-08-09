"""``C-fb``: closing the budget on realized work rather than on nominal work.

The previous campaign's weakest measurement was a leak.  Its budget was
equalised on the nominal schedule, :math:`\\int N^2\\dd t`, but an adaptive
integrator does not distribute its calls uniformly in time: it takes short
steps exactly where a concentrated schedule runs its highest degrees, so the
realized :math:`\\sum_{\\text{calls}}N^2` overran the nominal one by a median
of a third, and by more than a factor of two at half the budget.

Two accounting bases, and mixing them is the original error
-----------------------------------------------------------
:math:`\\int N^2\\dd t` is a *rate* integrated over time; :math:`\\sum N^2` is a
sum over calls.  They agree only if the call density is uniform, which is the
very assumption that failed.  A feedback law that compared spent calls against
a nominal time profile would therefore be driving on the difference between two
accounting conventions and would inject the original category error into the
reference signal --- a subtler place than where it was found, and no easier to
see (``DECISIONS.md`` D174).

:class:`ReferenceProfile` accordingly carries its ``basis`` and refuses to be
anonymous.  Both constructors exist because the campaign measures the
difference: :meth:`ReferenceProfile.from_realized_calls` is the deployable
signal, built from the pilot arc's own call history, and
:meth:`ReferenceProfile.from_nominal_schedule` is the control that shows how
much the mismatch is worth.

The law
-------
.. math::
    \\lambda(\\varphi)=\\operatorname{clip}\\Bigl(\\lambda_0+g\\,
    \\lambda_{\\mathrm{ref}}\\bigl[\\tfrac{W_{\\mathrm{spent}}}{B}
    -\\rho(\\varphi)\\bigr],\\;0,\\;\\lambda_{\\max}\\Bigr) ,

with :math:`\\rho` the reference profile and :math:`\\lambda_{\\mathrm{ref}}=
J_{\\mathrm{plan}}/B` a scale that makes the gain dimensionless.

Proportional on an accumulated error is integral action on the rate, so the
law has no steady-state rate offset --- which is what "the declared budget is
met end to end rather than per call" requires.  Adding an explicit integral
term on top would integrate an integral and is not what the manuscript claims.

The lower clip at zero is the KKT condition, not a guard: a multiplier is a
price on a ceiling and a negative one would pay the schedule to spend.  The
upper clip is a guard, it is counted, and a run that saturates is reported as
saturated rather than as controlled.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

__all__ = [
    "PROFILE_BASES",
    "BudgetFeedback",
    "ReferenceProfile",
    "WorkTracker",
    "overspend_ratio",
]

Arr = NDArray[np.float64]

PROFILE_BASES = ("realized", "nominal")
"""How a reference profile counts work: over calls, or over time."""


@dataclass(slots=True)
class WorkTracker:
    """Realized gravity work, accumulated over the calls actually made.

    Attributes
    ----------
    spent:
        :math:`\\sum_{\\text{calls}}N^2`.
    calls:
        How many right-hand-side evaluations produced it.

    Notes
    -----
    Mutable, deliberately: it is the controller's one piece of state, and
    threading a new frozen copy through every right-hand-side call would make
    the cost of measuring the budget depend on how the integrator is written.
    """

    spent: float = 0.0
    calls: int = 0

    def charge(self, degree: int, n_calls: int = 1) -> None:
        """Record ``n_calls`` evaluations at ``degree``.

        Raises
        ------
        ValueError
            If the degree is negative or the call count is not positive.
        """
        if degree < 0:
            raise ValueError(f"degree must be non-negative, got {degree}")
        if n_calls <= 0:
            raise ValueError(f"n_calls must be positive, got {n_calls}")
        self.spent += float(n_calls) * float(degree) ** 2
        self.calls += int(n_calls)

    def fraction(self, ceiling: float) -> float:
        """Share of a realized-work ceiling consumed so far.

        Raises
        ------
        ValueError
            If the ceiling is not positive.
        """
        if ceiling <= 0.0:
            raise ValueError(f"ceiling must be positive, got {ceiling}")
        return self.spent / ceiling


@dataclass(frozen=True, slots=True)
class ReferenceProfile:
    """Where the budget was supposed to have been by a given phase.

    Attributes
    ----------
    coordinate:
        ``revolution + phase``, ascending, shape ``(n,)``.
    fraction:
        Cumulative share of the total work at those coordinates, shape
        ``(n,)``, non-decreasing, running from zero to one.
    basis:
        One of :data:`PROFILE_BASES`.  Carried through every construction and
        reported, because a profile whose basis is unknown cannot be read
        against a tracker.
    """

    coordinate: Arr
    fraction: Arr
    basis: str

    def __post_init__(self) -> None:
        """Reject a profile that cannot be interpolated or compared."""
        if self.basis not in PROFILE_BASES:
            raise ValueError(
                f"basis must be one of {PROFILE_BASES}, got {self.basis!r}")
        if self.coordinate.shape != self.fraction.shape:
            raise ValueError(
                f"coordinate {self.coordinate.shape} and fraction "
                f"{self.fraction.shape} must match")
        if self.coordinate.size < 2:
            raise ValueError("a profile needs at least two points")
        if np.any(np.diff(self.coordinate) <= 0.0):
            raise ValueError("coordinate must be strictly increasing")
        if np.any(np.diff(self.fraction) < 0.0):
            raise ValueError(
                "fraction must be non-decreasing; a profile that falls would "
                "ask the controller to un-spend work")

    def at(self, coordinate: float) -> float:
        """Planned cumulative fraction at a phase coordinate.

        Held at the end values outside the tabulated range.  A vehicle past
        the profile's last point should be at one hundred per cent of the
        plan, and that is what the clamp says.
        """
        return float(np.interp(coordinate, self.coordinate, self.fraction))

    @classmethod
    def from_realized_calls(cls, coordinate: Arr, degrees: Arr,
                            ) -> ReferenceProfile:
        """Build from a call history: one entry per right-hand-side call.

        The deployable signal.  The pilot arc is integrated by the same solver
        with the same step-size logic, so its call density is the best
        available prediction of the flown one, and it is measured rather than
        modelled.

        Parameters
        ----------
        coordinate:
            Phase coordinate of each call, non-decreasing.
        degrees:
            Degree evaluated at each call.

        Raises
        ------
        ValueError
            If the shapes disagree, the coordinates decrease, or no work was
            done at all.

        Notes
        -----
        Calls sharing a coordinate are collapsed onto one point carrying their
        combined work.  A multi-stage integrator can evaluate the right-hand
        side more than once at the same epoch, and requiring strict increase
        would reject a perfectly ordinary call history.  The retained value is
        the cumulative sum *after* every call at that coordinate, so the
        profile at a call's own coordinate includes it.
        """
        coordinate = np.asarray(coordinate, dtype=float)
        work = np.asarray(degrees, dtype=float) ** 2
        if work.shape != coordinate.shape:
            raise ValueError(
                f"degrees {work.shape} must match coordinate "
                f"{coordinate.shape}")
        if coordinate.ndim != 1:
            raise ValueError("the call history must be one-dimensional")
        if np.any(np.diff(coordinate) < 0.0):
            raise ValueError("the call history must be non-decreasing")
        total = float(work.sum())
        if total <= 0.0:
            raise ValueError("the call history carries no work")
        cumulative = np.cumsum(work) / total
        last = np.append(np.diff(coordinate) > 0.0, True)
        return cls(coordinate=coordinate[last], fraction=cumulative[last],
                   basis="realized")

    @classmethod
    def from_nominal_schedule(cls, edge_coordinate: Arr, time_weight: Arr,
                              degrees: Arr) -> ReferenceProfile:
        """Build from the nominal schedule, :math:`\\sum_qW_qN_q^2`.

        The **control**, not the default.  It assumes the calls are spread
        uniformly in time, which is exactly the assumption whose failure
        ``C-fb`` exists to absorb; using it as the reference reintroduces the
        mismatch in the set-point instead of in the plant.  Provided so the
        campaign can measure what that costs.

        Parameters
        ----------
        edge_coordinate:
            Decision-interval boundaries in phase, shape ``(K+1,)``.
        time_weight:
            :math:`W_q`, shape ``(K,)``.
        degrees:
            :math:`N_q`, shape ``(K,)``.
        """
        edge_coordinate = np.asarray(edge_coordinate, dtype=float)
        work = (np.asarray(time_weight, dtype=float)
                * np.asarray(degrees, dtype=float) ** 2)
        if edge_coordinate.shape != (work.size + 1,):
            raise ValueError(
                f"edge_coordinate must have shape ({work.size + 1},), got "
                f"{edge_coordinate.shape}")
        total = float(work.sum())
        if total <= 0.0:
            raise ValueError("the nominal schedule carries no work")
        fraction = np.zeros(work.size + 1, dtype=float)
        fraction[1:] = np.cumsum(work) / total
        return cls(coordinate=edge_coordinate, fraction=fraction,
                   basis="nominal")


@dataclass(slots=True)
class BudgetFeedback:
    """The ``C-fb`` loop: a tracker, a profile and the multiplier law.

    Attributes
    ----------
    ceiling:
        The realized-work ceiling :math:`B`, in the tracker's units.
    nominal_multiplier:
        :math:`\\lambda_0` from the plan.
    reference_multiplier:
        :math:`\\lambda_{\\mathrm{ref}}=J_{\\mathrm{plan}}/B`, the scale that
        makes ``gain`` dimensionless.
    profile:
        Where the spending was supposed to be.
    gain:
        :math:`g`.  Zero is the open-loop control and is how the campaign
        measures the leak the loop closes.
    max_factor:
        Upper clip, as a multiple of the larger of the two multipliers above.
    tracker:
        Realized work so far.
    saturations:
        How many decisions hit the upper clip.

    Notes
    -----
    The tracker's basis and the profile's must agree.  A ``nominal`` profile is
    accepted --- it is a declared control --- but is recorded on every
    decision, so a run cannot be read as deployable when it was not.
    """

    ceiling: float
    nominal_multiplier: float
    reference_multiplier: float
    profile: ReferenceProfile
    gain: float = 1.0
    max_factor: float = 1.0e3
    tracker: WorkTracker = field(default_factory=WorkTracker)
    saturations: int = 0

    def __post_init__(self) -> None:
        """Reject settings that would make the law meaningless."""
        if self.ceiling <= 0.0:
            raise ValueError(f"ceiling must be positive, got {self.ceiling}")
        if self.nominal_multiplier < 0.0:
            raise ValueError(
                f"lambda_0 must be non-negative, got "
                f"{self.nominal_multiplier}")
        if self.reference_multiplier <= 0.0:
            raise ValueError(
                "reference_multiplier must be positive; it is the scale the "
                "gain is measured in and a zero scale disables the loop "
                "silently")
        if self.gain < 0.0:
            raise ValueError(
                "gain must be non-negative; a negative gain would relax the "
                "price in response to an overspend")

    @property
    def upper_clip(self) -> float:
        """The largest multiplier the law may return."""
        return self.max_factor * max(self.reference_multiplier,
                                     self.nominal_multiplier)

    def charge(self, degree: int, n_calls: int = 1) -> None:
        """Record realized work; a thin pass-through to the tracker."""
        self.tracker.charge(degree, n_calls)

    def error(self, coordinate: float) -> float:
        """Overspend at a phase coordinate, as a fraction of the ceiling.

        Positive means ahead of plan.  Both terms are dimensionless shares of
        their own total, which is what makes them subtractable.
        """
        return (self.tracker.fraction(self.ceiling)
                - self.profile.at(coordinate))

    def multiplier_at(self, coordinate: float) -> float:
        """:math:`\\lambda` for the interval starting at this phase.

        Returns
        -------
        float
            Non-negative.  Saturation at the upper clip increments
            :attr:`saturations` and is reported; the value is still returned,
            because refusing to decide would leave the vehicle without a
            degree.
        """
        raw = (self.nominal_multiplier
               + self.gain * self.reference_multiplier * self.error(coordinate))
        upper = self.upper_clip
        if raw > upper:
            self.saturations += 1
            return upper
        return max(raw, 0.0)


def overspend_ratio(spent: float, calls: int, nominal_work: float,
                    duration: float) -> float:
    """Realized mean :math:`N^2` per call against the nominal time average.

    The quantity the previous campaign reported as :math:`B_2/B_1`, written so
    that the two sides are the same physical thing.  The numerator is the mean
    squared degree the integrator actually paid for, the denominator the mean
    squared degree the schedule declared; their ratio exceeds one exactly when
    the calls concentrated where the degree was high.

    Parameters
    ----------
    spent:
        :math:`\\sum_{\\text{calls}}N^2`.
    calls:
        Number of calls.
    nominal_work:
        :math:`B=\\sum_qW_qN_q^2`.
    duration:
        Arc length :math:`T`.

    Returns
    -------
    float

    Raises
    ------
    ValueError
        If any denominator is not positive.

    Examples
    --------
    A uniform call density on a constant schedule gives exactly one.

    >>> overspend_ratio(spent=100.0 * 4.0, calls=100, nominal_work=4.0 * 50.0,
    ...                 duration=50.0)
    1.0
    """
    if calls <= 0:
        raise ValueError(f"calls must be positive, got {calls}")
    if duration <= 0.0:
        raise ValueError(f"duration must be positive, got {duration}")
    if nominal_work <= 0.0:
        raise ValueError(f"nominal_work must be positive, got {nominal_work}")
    return (spent / calls) / (nominal_work / duration)
