"""Verification of the comparator schedules.

The properties that matter are the endpoint identities of the interpolation
family -- ``k=0`` must be the constant and ``k=1`` the radial rule, or the
family is not interpolating between the two things it claims to -- and the
budget contract, which is the tight end of the ceiling for a comparator and
must never be exceeded.

The previous campaign used this family and reported ``k=0.5``; the tests below
pin the construction rather than the numbers, since the numbers belong to the
orbits.
"""

from __future__ import annotations

import numpy as np
import pytest

from tda.policies import (
    RadialTable,
    fixed_family_envelope,
    interval_altitudes,
    nearest_constant,
    radial_family,
    snap_to_candidates,
)

CANDIDATES = tuple(range(20, 321, 10))
REFERENCE = 320
K_DEC = 24


@pytest.fixture
def setup():
    """A decision grid with a varying altitude profile and equal weights."""
    weight = np.full(K_DEC, 120.0)
    phase = np.linspace(0.0, 2.0 * np.pi, K_DEC, endpoint=False)
    altitude = 5.0e4 + 4.5e5 * (1.0 - np.cos(phase)) / 2.0    # 50..500 km
    table = RadialTable(
        bin_edges_m=np.arange(0.0, 6.0e5, 1.0e4),
        degrees=np.clip(np.round(300.0 * (5.0e4 / np.maximum(
            np.arange(0.0, 6.0e5, 1.0e4), 1.0e4)) ** 0.5), 20, 300
        ).astype(np.int64),
    )
    return weight, altitude, table


# ---------------------------------------------------------------------------
# Candidate snapping
# ---------------------------------------------------------------------------


def test_snap_picks_the_nearest_and_ties_go_low() -> None:
    got = snap_to_candidates(np.array([19.0, 24.0, 25.0, 26.0, 999.0]),
                             (20, 30, 40))
    assert list(got) == [20, 20, 20, 30, 40]


def test_snap_is_exact_on_members() -> None:
    grid = (20, 55, 300)
    assert list(snap_to_candidates(np.array(grid, dtype=float), grid)) == \
        list(grid)


@pytest.mark.parametrize(("candidates", "message"), [
    ((), "empty"),
    ((30, 20), "sorted and unique"),
    ((20, 20), "sorted and unique"),
])
def test_snap_rejects_a_bad_candidate_set(candidates, message) -> None:
    with pytest.raises(ValueError, match=message):
        snap_to_candidates(np.array([25.0]), candidates)


# ---------------------------------------------------------------------------
# F-op
# ---------------------------------------------------------------------------


def test_nearest_constant_matches_the_square_root_of_the_rate() -> None:
    """A constant degree costs ``T N^2``, so the target is ``sqrt(B/T)``."""
    duration = 600.0
    for degree in (50, 130, 300):
        assert nearest_constant(duration * degree**2, duration,
                                CANDIDATES) == degree


def test_nearest_constant_rounds_rather_than_floors() -> None:
    """It is the anchor: it defines the ceiling instead of obeying one.

    Flooring would move the ceiling whenever the candidate grid moved, which
    would make the anchor a property of the grid rather than of the budget.
    """
    duration = 100.0
    just_under = duration * 99.0**2
    assert nearest_constant(just_under, duration, CANDIDATES) == 100
    assert duration * 100.0**2 > just_under


@pytest.mark.parametrize(("budget", "duration"), [(0.0, 1.0), (1.0, 0.0)])
def test_nearest_constant_rejects_a_degenerate_budget(budget,
                                                      duration) -> None:
    with pytest.raises(ValueError, match="must be positive"):
        nearest_constant(budget, duration, CANDIDATES)


# ---------------------------------------------------------------------------
# The interpolation family
# ---------------------------------------------------------------------------


def test_k_zero_is_a_flat_profile(setup) -> None:
    """The constant endpoint: every interval gets the same degree."""
    weight, altitude, table = setup
    radial = table.at(altitude)
    result = radial_family(0.0, 150, radial, weight,
                           budget=weight.sum() * 150.0**2,
                           candidates=CANDIDATES, reference_degree=REFERENCE)
    assert len(set(result.degrees.tolist())) == 1


def test_k_one_reproduces_the_radial_shape(setup) -> None:
    """The radial endpoint: the profile is the rule's, up to one scale."""
    weight, altitude, table = setup
    radial = table.at(altitude)
    result = radial_family(1.0, 150, radial, weight,
                           budget=weight.sum() * 150.0**2,
                           candidates=CANDIDATES, reference_degree=REFERENCE)
    # Rank correlation, because the scale and the snap change the values but
    # must not reorder the intervals.
    order_rule = np.argsort(np.argsort(radial))
    order_got = np.argsort(np.argsort(result.degrees.astype(float)))
    assert np.corrcoef(order_rule, order_got)[0, 1] > 0.98


def test_span_is_multiplicative_in_k(setup) -> None:
    """The property the geometric blend was chosen for.

    ``span(k) = span(1)^k``, so ``k`` moves the aggressiveness on a scale
    whose endpoints are a flat profile and the rule's own span.  Checked
    before snapping, since rounding to a coarse candidate grid quantises the
    span and would test the grid rather than the family.
    """
    _, altitude, table = setup
    radial = table.at(altitude)
    full = radial.max() / radial.min()
    for k in (0.25, 0.5, 0.75):
        blended = 150.0 ** (1.0 - k) * radial**k
        assert blended.max() / blended.min() == pytest.approx(full**k,
                                                              rel=1e-12)


def test_the_interior_member_sits_between_the_endpoints(setup) -> None:
    """An interior member must be between the two endpoints in aggressiveness.

    Less concentrated than ``R-rad`` and less flat than the constant, or the
    name does not describe it.
    """
    weight, altitude, table = setup
    radial = table.at(altitude)
    budget = weight.sum() * 150.0**2
    spans = []
    for k in (0.0, 0.5, 1.0):
        deg = radial_family(k, 150, radial, weight, budget, CANDIDATES,
                            REFERENCE).degrees.astype(float)
        spans.append(deg.max() / deg.min())
    assert spans[0] < spans[1] < spans[2]


# ---------------------------------------------------------------------------
# The budget contract
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("k", [0.0, 0.25, 0.5, 0.75, 1.0])
def test_every_member_respects_the_ceiling(setup, k) -> None:
    """Never exceeded. A comparator over budget would win by cheating."""
    weight, altitude, table = setup
    radial = table.at(altitude)
    budget = weight.sum() * 150.0**2
    result = radial_family(k, 150, radial, weight, budget, CANDIDATES,
                           REFERENCE)
    assert result.work <= budget
    assert result.utilisation <= 1.0


@pytest.mark.parametrize("k", [0.0, 0.5, 1.0])
def test_the_ceiling_is_approached_rather_than_avoided(setup, k) -> None:
    """A comparator has no reason to leave budget unspent.

    The ceiling contract lets an *optimiser* stop short, since the objective
    is not monotone in degree.  These policies do not optimise it, so leaving
    slack would only make them artificially weak -- the calibration takes the
    largest scale that fits.
    """
    weight, altitude, table = setup
    radial = table.at(altitude)
    budget = weight.sum() * 150.0**2
    result = radial_family(k, 150, radial, weight, budget, CANDIDATES,
                           REFERENCE)
    assert result.utilisation > 0.90
    assert result.scale_limited


def test_work_uses_the_time_weights_not_a_count(setup) -> None:
    """Deliberate divergence from the archive (D60).

    On a refined grid a plain mean charges perilune for most of the arc
    merely because it is sampled most densely.  Doubling one interval's time
    weight must change the calibration; doubling a count must not enter.
    """
    _, altitude, table = setup
    radial = table.at(altitude)
    even = np.full(K_DEC, 120.0)
    skewed = even.copy()
    skewed[radial.argmax()] *= 8.0        # more time at the expensive end

    budget = even.sum() * 150.0**2
    a = radial_family(1.0, 150, radial, even, budget, CANDIDATES, REFERENCE)
    b = radial_family(1.0, 150, radial, skewed, skewed.sum() * 150.0**2,
                      CANDIDATES, REFERENCE)
    assert a.scale != b.scale


def test_an_impossible_budget_is_reported_not_absorbed(setup) -> None:
    """Below the smallest candidate there is no schedule at all.

    Saying so is a statement about the budget rather than about the policy,
    and absorbing it would return an infeasible schedule that looks fine.
    """
    weight, altitude, table = setup
    radial = table.at(altitude)
    with pytest.raises(ValueError, match="no scale fits the budget"):
        radial_family(0.5, 150, radial, weight, budget=1.0,
                      candidates=CANDIDATES, reference_degree=REFERENCE)


def test_the_reference_degree_clamps_the_schedule(setup) -> None:
    """A policy touching the reference degree has a zero defect there."""
    weight, altitude, table = setup
    radial = table.at(altitude)
    result = radial_family(1.0, 300, radial, weight,
                           budget=weight.sum() * 1.0e6,
                           candidates=CANDIDATES, reference_degree=100)
    assert result.degrees.max() <= 100


@pytest.mark.parametrize(("k", "message"), [(-0.1, r"\[0, 1\]"),
                                            (1.5, r"\[0, 1\]")])
def test_k_outside_the_family_is_refused(setup, k, message) -> None:
    weight, altitude, table = setup
    with pytest.raises(ValueError, match=message):
        radial_family(k, 150, table.at(altitude), weight, 1.0e9, CANDIDATES,
                      REFERENCE)


# ---------------------------------------------------------------------------
# Radial table and altitude reduction
# ---------------------------------------------------------------------------


def test_radial_table_clamps_rather_than_extrapolates() -> None:
    """Extrapolating a binned table would invent a rule never flown."""
    table = RadialTable(bin_edges_m=np.array([0.0, 1.0e4, 2.0e4]),
                        degrees=np.array([300, 200, 100]))
    assert list(table.at(np.array([-5.0e3, 0.0, 1.5e4, 9.9e5]))) == \
        [300.0, 300.0, 200.0, 100.0]


@pytest.mark.parametrize(("edges", "degrees", "message"), [
    (np.array([0.0, 1.0]), np.array([1, 2, 3]), "must match"),
    (np.array([]), np.array([]), "empty"),
    (np.array([1.0, 0.0]), np.array([1, 2]), "ascending"),
])
def test_radial_table_rejects_a_bad_definition(edges, degrees,
                                               message) -> None:
    with pytest.raises(ValueError, match=message):
        RadialTable(bin_edges_m=edges, degrees=degrees)


def test_interval_altitude_is_time_weighted() -> None:
    """A midpoint sample would follow wherever the refinement put a node."""
    radius = np.array([1.80e6, 1.80e6, 1.90e6])
    widths = np.array([10.0, 10.0, 1.0])
    interval_of = np.zeros(3, dtype=np.int64)
    got = interval_altitudes(radius, widths, interval_of, 1, 1.7374e6)
    expected = (10 * 1.80e6 + 10 * 1.80e6 + 1 * 1.90e6) / 21.0 - 1.7374e6
    assert got[0] == pytest.approx(expected)


def test_an_empty_interval_is_refused() -> None:
    with pytest.raises(ValueError, match="no accumulation cell"):
        interval_altitudes(np.array([1.8e6]), np.array([1.0]),
                           np.array([0]), 2, 1.7374e6)


# ---------------------------------------------------------------------------
# F-env
# ---------------------------------------------------------------------------


def test_envelope_picks_the_smallest_error() -> None:
    degree, error = fixed_family_envelope((100, 200, 300),
                                          np.array([5.0, 2.0, 7.0]))
    assert (degree, error) == (200, 2.0)


def test_envelope_cannot_be_won_by_an_unmeasured_cell() -> None:
    """A censored cell must not win by being unmeasured.

    Treating a non-finite error as small would hand the envelope to whichever
    degree failed to produce a number, which is the opposite of what a lower
    envelope means.
    """
    degree, _ = fixed_family_envelope((100, 200, 300),
                                      np.array([5.0, np.nan, 7.0]))
    assert degree == 100


def test_envelope_needs_at_least_one_measurement() -> None:
    with pytest.raises(ValueError, match="no finite error"):
        fixed_family_envelope((100, 200), np.array([np.nan, np.inf]))


def test_envelope_length_must_match() -> None:
    with pytest.raises(ValueError, match="must match"):
        fixed_family_envelope((100, 200), np.array([1.0]))


def test_a_one_parameter_family_cannot_generally_saturate(setup) -> None:
    """The scale moves every interval together, so the last unit is stranded.

    Not a defect and not fixable inside the family: closing the gap would mean
    raising one interval on its own, which is a schedule the previous campaign
    never flew.  It is reported instead, which is possible only because the
    contract is a ceiling with realized work recorded beside it rather than an
    equality (D142).
    """
    weight, altitude, table = setup
    radial = table.at(altitude)
    budget = weight.sum() * 150.0**2
    result = radial_family(0.5, 150, radial, weight, budget, CANDIDATES,
                           REFERENCE)

    assert result.scale_limited, "the scale must be at its largest feasible"
    assert not result.step_limited, "and yet one interval could still rise"
    assert 0.95 < result.utilisation < 1.0
