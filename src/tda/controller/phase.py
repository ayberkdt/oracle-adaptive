"""Indexing the plan by orbital phase rather than by absolute time.

A pilot arc integrated at degree 40 with loose tolerances drifts along track
against the arc that is actually flown.  The drift is a *timing* offset ---
kilometres of along-track error at low-lunar-orbit speed is of order a hundred
seconds --- and a perilune passage lasts a few minutes.  Since the cancellation
target and the local kernel are dominated by exactly those passages, a plan
indexed by absolute time would apply the perilune correction before or after
the perilune it was computed for: an error of order the thing being corrected
(``DECISIONS.md`` D52).

The index is anchored at apolune
--------------------------------
Revolutions are counted apolune to apolune, not perilune to perilune, so that
the dominant feature sits near phase :math:`0.5` and is never split across the
wrap.  The alternative puts the one passage the plan cares about exactly at the
seam, where any interpolation has to cross a discontinuity in the index.

Where the index does not exist, it is not needed
------------------------------------------------
On a near-circular arc the apsides are ill-conditioned: the radius oscillates
under the harmonics rather than under the conic, so :math:`\\dot r` changes
sign many times per revolution and no crossing is the apolune.
:class:`RevolutionIndex` refuses to build there instead of returning an index
built from noise.  That refusal costs nothing, and the reason is worth stating:
a near-circular arc has no perilune passage for the plan to misalign, so the
failure mode the phase index exists to prevent is absent in precisely the
regime where the index cannot be constructed.  A time index is admissible
there, and the campaign records which arcs used which.

References
----------
.. [Vallado2013] D. A. Vallado, *Fundamentals of Astrodynamics and
   Applications*, 4th ed., Microcosm, 2013, §2.2 -- the osculating semi-major
   axis used for the independent period estimate.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

__all__ = [
    "PhaseIndexError",
    "RevolutionIndex",
    "apsis_epochs",
    "osculating_period",
    "radial_rate",
]

Arr = NDArray[np.float64]
IntArr = NDArray[np.int_]

MIN_RADIAL_VARIATION = 1.0e-3
"""Relative peak-to-peak radius below which the apsides are declared absent.

One part in a thousand is a kilometre at a lunar orbit radius, which is the
scale on which the harmonics move the radius; below it a sign change of
:math:`\\dot r` carries no information about the conic.
"""


class PhaseIndexError(RuntimeError):
    """Raised when an arc cannot carry an apsidal phase index.

    Not absorbed into a fallback.  A silently degraded index would look like a
    working one and would misalign the plan by an unknown amount, which is the
    failure the index exists to prevent.
    """


def radial_rate(states: Arr) -> Arr:
    """:math:`\\dot r=\\mathbf r\\cdot\\mathbf v/\\lVert\\mathbf r\\rVert`.

    Parameters
    ----------
    states:
        Cartesian states ``[r, v]``, shape ``(n, 6)``.

    Returns
    -------
    ndarray, shape (n,)

    Raises
    ------
    ValueError
        If the states are misshaped or one of them is at the origin.

    Examples
    --------
    >>> import numpy as np
    >>> radial_rate(np.array([[1.0, 0, 0, 0.5, 2.0, 0]]))
    array([0.5])
    """
    states = np.asarray(states, dtype=float)
    if states.ndim != 2 or states.shape[1] != 6:
        raise ValueError(f"states must have shape (n, 6), got {states.shape}")
    radius = np.linalg.norm(states[:, 0:3], axis=1)
    if np.any(radius == 0.0):
        raise ValueError("a state is at the origin; the radial rate is "
                         "undefined there")
    return np.einsum("ia,ia->i", states[:, 0:3], states[:, 3:6]) / radius


def osculating_period(state: Arr, mu: float) -> float:
    """Keplerian period of the osculating ellipse, seconds.

    Used only to set the minimum spacing two apsides of the same kind may
    have.  Taking it from the state rather than from the detected crossings is
    what keeps the filter from bootstrapping off the contamination it is
    meant to remove: a median over a list containing spurious apsides is
    itself biased by them.

    Parameters
    ----------
    state:
        ``[r, v]``, shape ``(6,)``.
    mu:
        Gravitational parameter.

    Returns
    -------
    float

    Raises
    ------
    PhaseIndexError
        If the state is not bound.  An unbound arc has one perilune and no
        revolutions, so there is nothing to index.
    """
    state = np.asarray(state, dtype=float)
    if state.shape != (6,):
        raise ValueError(f"state must have shape (6,), got {state.shape}")
    radius = float(np.linalg.norm(state[0:3]))
    alpha = 2.0 / radius - float(state[3:6] @ state[3:6]) / mu
    if alpha <= 0.0:
        raise PhaseIndexError(
            "the state is unbound (1/a <= 0); an escape arc has no "
            "revolutions to index")
    return 2.0 * math.pi * math.sqrt((1.0 / alpha) ** 3 / mu)


def apsis_epochs(times: Arr, rate: Arr, *, falling: bool) -> Arr:
    """Epochs where :math:`\\dot r` changes sign, by linear interpolation.

    Parameters
    ----------
    times:
        Sample epochs, strictly increasing, shape ``(n,)``.
    rate:
        :math:`\\dot r` at those epochs, shape ``(n,)``.
    falling:
        ``True`` selects downward crossings (apolune), ``False`` upward ones
        (perilune).

    Returns
    -------
    ndarray
        Crossing epochs, ascending.

    Raises
    ------
    ValueError
        If the shapes disagree or the epochs are not increasing.

    Notes
    -----
    The crossing is interpolated rather than snapped to the nearer sample.
    Snapping would quantise the phase origin to the sampling interval, and the
    whole point of the index is to align two arcs to better than the
    hundred-second drift it is correcting --- a coarser origin would put the
    quantisation error back on top of the offset it removes.

    A sample exactly at zero is treated as belonging to the interval that
    *ends* there, so a crossing is reported once and not twice.
    """
    times = np.asarray(times, dtype=float)
    rate = np.asarray(rate, dtype=float)
    if times.shape != rate.shape or times.ndim != 1:
        raise ValueError(
            f"times {times.shape} and rate {rate.shape} must be equal 1-D")
    if times.size < 2:
        raise ValueError("at least two samples are needed to find a crossing")
    if np.any(np.diff(times) <= 0.0):
        raise ValueError("times must be strictly increasing")

    before, after = rate[:-1], rate[1:]
    if falling:
        crossing = (before > 0.0) & (after <= 0.0)
    else:
        crossing = (before < 0.0) & (after >= 0.0)
    idx = np.flatnonzero(crossing)
    if idx.size == 0:
        return np.empty(0, dtype=float)
    span = after[idx] - before[idx]
    fraction = np.where(span == 0.0, 0.0, -before[idx] / span)
    return times[idx] + fraction * (times[idx + 1] - times[idx])


@dataclass(frozen=True, slots=True)
class RevolutionIndex:
    """Maps an epoch to a continuous ``revolution + phase`` coordinate.

    Attributes
    ----------
    anchors:
        Apolune epochs, ascending, shape ``(R+1,)``; ``anchors[k]`` opens
        revolution ``k``.
    spurious:
        Downward crossings discarded as too close to their predecessor.
        Reported rather than dropped silently: a large count means the arc is
        closer to the near-circular regime than the variation test suggested,
        and the plan built on it should be read accordingly.
    period_s:
        The independent period estimate the filter used.

    Notes
    -----
    The coordinate is continuous and strictly increasing across the whole arc,
    including the partial revolutions before the first anchor and after the
    last, where the phase is extrapolated with the adjacent full span.  Making
    it one monotone real rather than a pair turns the plan lookup into a plain
    one-dimensional search and removes the wrap case entirely.
    """

    anchors: Arr
    spurious: int
    period_s: float

    def __len__(self) -> int:
        """Number of complete revolutions the index spans."""
        return int(self.anchors.size - 1)

    @property
    def spans(self) -> Arr:
        """Duration of each complete revolution, shape ``(R,)``."""
        return np.diff(self.anchors)

    @classmethod
    def from_arc(cls, times: Arr, states: Arr, mu: float,
                 min_variation: float = MIN_RADIAL_VARIATION,
                 ) -> RevolutionIndex:
        """Build the index from a sampled arc.

        Parameters
        ----------
        times:
            Sample epochs, strictly increasing.
        states:
            ``[r, v]`` at those epochs, shape ``(n, 6)``.
        mu:
            Gravitational parameter, for the independent period estimate.
        min_variation:
            Relative peak-to-peak radius below which the arc is refused.

        Returns
        -------
        RevolutionIndex

        Raises
        ------
        PhaseIndexError
            If the arc is too circular for apsides to mean anything, if it is
            unbound, or if it does not contain two apolunes --- an arc shorter
            than one revolution has no phase the plan can be indexed by.
        """
        states = np.asarray(states, dtype=float)
        radius = np.linalg.norm(states[:, 0:3], axis=1)
        variation = float(radius.max() - radius.min()) / float(radius.mean())
        if variation < min_variation:
            raise PhaseIndexError(
                f"peak-to-peak radius is {variation:.3e} of the mean, below "
                f"{min_variation:.3e}: the apsides of this arc are set by the "
                "harmonics rather than by the conic. Index this arc by time; "
                "it has no perilune passage for a phase index to protect.")

        period = osculating_period(states[0], mu)
        raw = apsis_epochs(times, radial_rate(states), falling=True)
        kept: list[float] = []
        spurious = 0
        for epoch in raw:
            if kept and epoch - kept[-1] < 0.5 * period:
                spurious += 1
                continue
            kept.append(float(epoch))
        if len(kept) < 2:
            raise PhaseIndexError(
                f"found {len(kept)} apolune(s) in {times[-1] - times[0]:.6g} s "
                f"against an osculating period of {period:.6g} s; an arc "
                "shorter than one revolution cannot be phase-indexed")
        return cls(anchors=np.array(kept, dtype=float), spurious=spurious,
                   period_s=period)

    def coordinate(self, epoch: Arr) -> Arr:
        """Continuous ``revolution + phase`` at one or more epochs.

        Parameters
        ----------
        epoch:
            Scalar or array of epochs, seconds.

        Returns
        -------
        ndarray
            Same shape as ``epoch``; integer part is the revolution, fractional
            part the phase within it.

        Examples
        --------
        >>> import numpy as np
        >>> index = RevolutionIndex(np.array([0.0, 10.0, 20.0]), 0, 10.0)
        >>> index.coordinate(np.array([0.0, 5.0, 15.0, 25.0]))
        array([0. , 0.5, 1.5, 2.5])
        """
        epoch = np.asarray(epoch, dtype=float)
        anchors, spans = self.anchors, self.spans
        slot = np.clip(np.searchsorted(anchors, epoch, side="right") - 1,
                       0, spans.size - 1)
        return slot + (epoch - anchors[slot]) / spans[slot]

    def epoch(self, coordinate: Arr) -> Arr:
        """Inverse of :meth:`coordinate`.

        Extrapolates outside the anchored range with the adjacent span, so it
        is the exact inverse everywhere :meth:`coordinate` is defined.
        """
        coordinate = np.asarray(coordinate, dtype=float)
        anchors, spans = self.anchors, self.spans
        slot = np.clip(np.floor(coordinate).astype(np.int64), 0, spans.size - 1)
        return anchors[slot] + (coordinate - slot) * spans[slot]

    def covers(self, coordinate: Arr) -> NDArray[np.bool_]:
        """Whether a coordinate falls inside the anchored revolutions.

        The extrapolated tails are usable but are not evidence, and a caller
        that has drifted past the plan's last anchor is in a different
        situation from one inside it.  Reported separately so the campaign can
        count the intervals decided on an extrapolated phase.
        """
        coordinate = np.asarray(coordinate, dtype=float)
        return (coordinate >= 0.0) & (coordinate <= float(len(self)))
