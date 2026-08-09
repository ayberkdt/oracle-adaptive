"""Verification of the defect table and its provenance.

Two things matter here and they are different in kind.  The arithmetic --- do
the stored defects and the recomputed transports produce the contributions the
kernel expects --- is checked against a direct construction.  The provenance
--- can a table be paired with a grid or a configuration it was not built on
--- is checked by trying to do exactly that and requiring a refusal, because a
silently mismatched table would misalign every cell and produce numbers that
look entirely reasonable.

The field is analytic, so the suite runs without the synthesis kernel.
"""

from __future__ import annotations

import json

import numpy as np
import pytest

from conftest import J2_MOON, MU_MOON, R_MOON
from tda.analytic import J2Field
from tda.config import GridConfig
from tda.grids import build_accumulation_grid, build_decision_edges
from tda.kernel import CouplingKernel
from tda.stm import position_rows, symplectic_inverse, velocity_transport
from tda.tables import (
    DefectTable,
    build_defect_table,
    grid_digest,
    load_table,
)

DURATION = 1200.0
DEGREES = (0, 1, 2)
REFERENCE = 2


@pytest.fixture
def built(tmp_path):
    """A small table built on an analytic arc, plus everything it came from."""
    cfg = GridConfig(samples_per_tau=2.0, dt_acc_min_s=20.0,
                     dt_acc_max_s=120.0, dt_dec_s=300.0)
    sample_t = np.linspace(0.0, DURATION, 601)
    radius = np.full_like(sample_t, R_MOON + 1.0e5)
    speed = np.full_like(sample_t, np.sqrt(MU_MOON / (R_MOON + 1.0e5)))
    edges = build_decision_edges(DURATION, cfg)
    grid = build_accumulation_grid(sample_t, radius, speed, 60, cfg, edges)

    epochs = grid.required_epochs()
    rng = np.random.default_rng(11)
    # A fabricated arc: positions on a circle, random invertible transports.
    # Random is stronger than physical here -- a real arc might satisfy an
    # identity the bookkeeping does not.
    angle = 2.0 * np.pi * epochs / DURATION
    r0 = R_MOON + 1.0e5
    positions = np.stack([r0 * np.cos(angle), r0 * np.sin(angle),
                          0.3 * r0 * np.sin(2.0 * angle)], axis=1)
    stm = np.stack([np.eye(6)]
                   + [rng.normal(size=(6, 6)) for _ in range(epochs.size - 1)])

    field = J2Field(mu=MU_MOON, reference_radius=R_MOON, j2=J2_MOON)
    table = build_defect_table(
        tmp_path, field, grid, epochs, positions, stm,
        DEGREES, REFERENCE, orbit="orb000", config_digest="deadbeefdeadbeef")
    return table, grid, epochs, positions, stm, field, tmp_path


# ---------------------------------------------------------------------------
# Contents
# ---------------------------------------------------------------------------


def test_shape_and_schema_agree(built) -> None:
    table, grid, *_ = built
    assert table.defect.shape == (len(grid), len(DEGREES), 3)
    assert table.schema.shape == table.defect.shape
    assert table.schema.candidate_degrees == DEGREES
    assert table.schema.n_cells == len(grid)
    assert len(table) == len(grid)


def test_defect_at_the_reference_degree_is_exactly_zero(built) -> None:
    """Not approximately: the column is written as zero and no synthesis spent.

    A tiny non-zero there would look like signal at the one degree where none
    can exist, and it would be the degree the capture fraction divides by.
    """
    table, *_ = built
    column = table.index_of(REFERENCE)
    assert np.array_equal(table.defect[:, column, :],
                          np.zeros((len(table), 3)))


def test_stored_defect_matches_a_direct_evaluation(built) -> None:
    table, grid, epochs, positions, _, field, _ = built
    nodes = grid.node_indices()
    for i in (0, len(table) // 2, len(table) - 1):
        for p, degree in enumerate(DEGREES):
            direct = field.defect(positions[nodes[i]], float(epochs[nodes[i]]),
                                  degree, REFERENCE)
            assert np.allclose(table.defect[i, p], direct, rtol=1e-13,
                               atol=0.0), (i, degree)


def test_synthesis_count_excludes_the_reference_column(built) -> None:
    """The measured cost the manifest records, not an estimate of it."""
    table, grid, *_ = built
    assert table.schema.syntheses == len(grid) * (len(DEGREES) - 1)


# ---------------------------------------------------------------------------
# Contributions
# ---------------------------------------------------------------------------


def test_contributions_match_the_definition(built) -> None:
    """``u_i = Phi(t_0, m_i) B Delta_a Delta_t``, built by hand."""
    table, grid, _, _, stm, _, _ = built
    schedule = np.full(len(table), 1, dtype=int)
    got = table.contributions(schedule)

    inv = symplectic_inverse(stm[grid.node_indices()])
    transport = velocity_transport(inv)
    column = table.index_of(1)
    expected = np.stack([
        transport[i] @ table.defect[i, column] * grid.widths[i]
        for i in range(len(table))
    ])
    assert np.allclose(got, expected, rtol=1e-12, atol=0.0)


def test_constant_schedule_shortcut_agrees_with_the_gather(built) -> None:
    table, *_ = built
    for degree in DEGREES:
        gathered = table.contributions(np.full(len(table), degree, dtype=int))
        direct = table.contributions_at(degree)
        assert np.allclose(gathered, direct, rtol=0.0, atol=0.0)


def test_reference_schedule_gives_zero_objective(built) -> None:
    """The null the whole machine must pass: no truncation, no error.

    Runs the table through the kernel, so it also checks that the transports
    the table carries are the ones the kernel expects.
    """
    table, grid, *_ = built
    kernel = CouplingKernel.from_arc(table.edge_rows, grid.outer_weights,
                                     grid.duration)
    u = table.contributions_at(REFERENCE)
    assert np.array_equal(u, np.zeros_like(u))
    assert kernel.objective(u) == 0.0


def test_the_table_reproduces_the_field_own_degree_structure(built) -> None:
    """Degrees 0 and 1 must agree, and both must differ from the reference.

    A centred field has no degree-one harmonic, so ``J2Field`` returns the
    same acceleration for both and the two columns of the table must be
    identical.  Asserting that the objective *decreases* with degree would be
    asserting a variation the fixture does not contain -- and would pass only
    because the table invented one.
    """
    table, grid, *_ = built
    kernel = CouplingKernel.from_arc(table.edge_rows, grid.outer_weights,
                                     grid.duration)
    j0 = kernel.objective(table.contributions_at(0))
    j1 = kernel.objective(table.contributions_at(1))

    assert np.array_equal(table.defect[:, table.index_of(0), :],
                          table.defect[:, table.index_of(1), :])
    assert j0 == pytest.approx(j1, rel=0.0, abs=0.0)
    assert j0 > 0.0


def test_edge_rows_are_the_position_rows_of_the_edge_transports(built) -> None:
    table, grid, _, _, stm, _, _ = built
    assert np.allclose(table.edge_rows,
                       position_rows(stm[grid.edge_indices()]),
                       rtol=0.0, atol=0.0)


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------


def test_schema_is_written_and_round_trips(built) -> None:
    table, grid, _, _, stm, _, tmp_path = built
    raw = json.loads((tmp_path / "orb000.schema.json").read_text("utf-8"))
    assert raw["orbit"] == "orb000"
    assert raw["grid_digest"] == grid_digest(grid)

    again = load_table(tmp_path, "orb000", grid, stm,
                       expect_config_digest="deadbeefdeadbeef")
    assert again.schema == table.schema
    assert np.array_equal(again.defect, table.defect)


def test_a_table_cannot_be_paired_with_another_grid(built) -> None:
    """The failure that would misalign every cell and still look plausible."""
    table, _, _, _, stm, _, tmp_path = built
    cfg = GridConfig(samples_per_tau=4.0, dt_acc_min_s=20.0,
                     dt_acc_max_s=120.0, dt_dec_s=300.0)
    sample_t = np.linspace(0.0, DURATION, 601)
    other = build_accumulation_grid(
        sample_t, np.full_like(sample_t, R_MOON + 1.0e5),
        np.full_like(sample_t, np.sqrt(MU_MOON / (R_MOON + 1.0e5))),
        60, cfg, build_decision_edges(DURATION, cfg))

    with pytest.raises(ValueError, match="different accumulation grid"):
        load_table(tmp_path, "orb000", other,
                   np.zeros((other.required_epochs().size, 6, 6)))


def test_a_wrong_config_digest_is_refused(built) -> None:
    _, grid, _, _, stm, _, tmp_path = built
    with pytest.raises(ValueError, match="config digest"):
        load_table(tmp_path, "orb000", grid, stm,
                   expect_config_digest="0000000000000000")


def test_a_stale_schema_version_is_refused(built) -> None:
    """Adapting an old layout is exactly what must not happen silently."""
    _, grid, _, _, stm, _, tmp_path = built
    path = tmp_path / "orb000.schema.json"
    raw = json.loads(path.read_text("utf-8"))
    raw["schema_version"] = 0
    path.write_text(json.dumps(raw), encoding="utf-8")
    with pytest.raises(ValueError, match="schema version"):
        load_table(tmp_path, "orb000", grid, stm)


def test_a_missing_table_is_not_silently_rebuilt(built) -> None:
    _, grid, _, _, stm, _, tmp_path = built
    (tmp_path / "orb000.schema.json").unlink()
    with pytest.raises(FileNotFoundError):
        load_table(tmp_path, "orb000", grid, stm)


def test_transports_are_recomputed_not_stored(built) -> None:
    """What lets the Q13 panel reuse an expensive table with a different Phi.

    Reloading with different state-transition matrices must give different
    transports and the *same* defects -- that separation is the reason the
    table stores the defect rather than the contribution.
    """
    table, grid, _, _, _, _, tmp_path = built
    rng = np.random.default_rng(77)
    other_stm = np.stack([np.eye(6)] + [
        rng.normal(size=(6, 6))
        for _ in range(grid.required_epochs().size - 1)])

    reloaded = load_table(tmp_path, "orb000", grid, other_stm)
    assert np.array_equal(reloaded.defect, table.defect)
    assert not np.allclose(reloaded.node_transport, table.node_transport)


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------


def test_an_arc_off_the_required_epochs_is_refused(built) -> None:
    """Interpolating a trajectory the campaign can integrate is not allowed."""
    _, grid, epochs, positions, stm, field, tmp_path = built
    shifted = epochs.copy()
    shifted[3] += 0.5
    with pytest.raises(ValueError, match="required_epochs"):
        build_defect_table(tmp_path, field, grid, shifted, positions, stm,
                           DEGREES, REFERENCE, "orb001", "d")


@pytest.mark.parametrize(("degrees", "reference", "message"), [
    ((), 2, "empty"),
    ((0, 1, 5), 2, "exceeds the reference degree"),
])
def test_a_bad_candidate_set_is_refused(built, degrees, reference,
                                        message) -> None:
    _, grid, epochs, positions, stm, field, tmp_path = built
    with pytest.raises(ValueError, match=message):
        build_defect_table(tmp_path, field, grid, epochs, positions, stm,
                           degrees, reference, "orb002", "d")


@pytest.mark.parametrize("bad", ["positions", "stm"])
def test_a_misshaped_arc_is_refused(built, bad) -> None:
    _, grid, epochs, positions, stm, field, tmp_path = built
    args = {"arc_positions": positions, "arc_stm": stm}
    args[f"arc_{bad}"] = args[f"arc_{bad}"][:-1]
    with pytest.raises(ValueError, match="must have shape"):
        build_defect_table(tmp_path, field, grid, epochs,
                           args["arc_positions"], args["arc_stm"],
                           DEGREES, REFERENCE, "orb003", "d")


def test_an_unknown_degree_is_a_key_error_not_a_neighbour(built) -> None:
    """A policy asking for an untabulated degree is a configuration error."""
    table, *_ = built
    with pytest.raises(KeyError, match="not in the candidate set"):
        table.index_of(999)


def test_a_misshaped_schedule_is_refused(built) -> None:
    table, *_ = built
    with pytest.raises(ValueError, match="must have shape"):
        table.contributions(np.zeros(3, dtype=int))


def test_the_defect_map_is_read_only(built) -> None:
    """A table is an input. Writing to one would fork the code path."""
    table, *_ = built
    assert isinstance(table, DefectTable)
    with pytest.raises(ValueError, match="read-only|assignment destination"):
        table.defect[0, 0, 0] = 1.0
