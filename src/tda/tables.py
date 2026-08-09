"""The tabulated truncation defect and its transports.

Turning the defect into the contributions the objective sums over.

This is the campaign's largest single computation and its only large artefact.
Building it costs :math:`M\\lvert\\mathcal N\\rvert` field syntheses --- some
three million per orbit at degree 300 --- which is more than the variational
integration that produces the state-transition matrices, and it dominates
stage M1 (``DECISIONS.md`` D137).

What is stored, and why not the obvious thing
---------------------------------------------
The table holds the **defect** :math:`\\Delta\\mathbf a(m_i,N)`, not the
transported contribution :math:`\\mathbf u_i(N)=\\Phi(t_0,m_i)\\mathbf B\\,
\\Delta\\mathbf a\\,\\Delta t_i`.  Three reasons, in increasing order of
importance:

* it is half the size, three components rather than six;
* it is the physical quantity --- the transport is a property of the arc, not
  of the candidate degree;
* and it separates the two things the campaign varies independently.  The
  open question of which gradient degree generates :math:`\\Phi`
  (``DECISIONS.md`` Q13) is answered by a convergence panel that changes the
  transport and leaves the defect alone.  Storing contributions would make
  that panel rebuild the three-million-synthesis table each time it moves the
  gradient degree; storing defects makes it recompute a few megabytes.

The transports are kept in the factored forms the objective actually uses ---
:math:`\\mathbf H_r\\Phi(e_j,t_0)`, which is :math:`3\\times6`, and
:math:`\\Phi(t_0,m_i)\\mathbf B`, which is :math:`6\\times3` --- rather than as
full state-transition matrices.  That halves their memory and makes the
positive semidefiniteness of :math:`\\mathbf M_j` structural rather than
incidental.

Layout
------
``(M, P, 3)`` in C order, so that all :math:`P` candidates of one cell are
contiguous.  That is the access pattern of the coordinate sweep, which
evaluates every candidate of one decision interval before moving on.  The
transposed layout would stride the whole table once per candidate.

References
----------
.. [Pines1973] S. Pines, "Uniform representation of the gravitational
   potential and its derivatives", *AIAA Journal* 11(11), 1973 -- the
   synthesis this module spends its time in.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from tda.field import GravityField
from tda.grids import AccumulationGrid
from tda.stm import position_rows, symplectic_inverse, velocity_transport

__all__ = ["DefectTable", "TableSchema", "build_defect_table", "load_table"]

Arr = NDArray[np.float64]

SCHEMA_VERSION = 1
"""Bumped when the on-disk layout changes; a mismatch refuses to load."""


@dataclass(frozen=True, slots=True)
class TableSchema:
    """Everything needed to know what a table on disk actually contains.

    Written beside the array as JSON.  A table whose schema does not match
    what the caller asked for is not silently reused: the campaign's rule is
    that every number comes from one code path, and a stale table is a second
    one.

    Attributes
    ----------
    orbit:
        Identifier of the arc, carried so a mixed-up file is caught by name.
    reference_degree:
        Adopted reference degree the defect is measured against.
    candidate_degrees:
        The set :math:`\\mathcal N`, sorted and unique.
    n_cells:
        :math:`M`.  Together with the candidate count it fixes the shape.
    dtype:
        ``float64`` by default; ``float32`` exists only so that the parity
        test of ``DECISIONS.md`` D60 can be run, and is never the default.
    config_digest, grid_digest:
        Provenance.  The first ties the table to the run configuration, the
        second to the exact edge sequence, so that a table cannot be paired
        with a grid it was not built on.
    syntheses:
        Field evaluations actually performed.  Measured, not estimated: it is
        the number the manifest records for the cost of stage M1.
    schema_version:
        On-disk layout version.
    """

    orbit: str
    reference_degree: int
    candidate_degrees: tuple[int, ...]
    n_cells: int
    dtype: str
    config_digest: str
    grid_digest: str
    syntheses: int
    schema_version: int = SCHEMA_VERSION

    @property
    def shape(self) -> tuple[int, int, int]:
        """Shape of the defect array on disk."""
        return (self.n_cells, len(self.candidate_degrees), 3)

    @property
    def nbytes(self) -> int:
        """Size of the defect array on disk."""
        return int(np.prod(self.shape)) * np.dtype(self.dtype).itemsize


def grid_digest(grid: AccumulationGrid, length: int = 16) -> str:
    """Stable short hash of an accumulation grid's edges.

    Hashing the edges rather than the settings that produced them: two runs
    can share a :class:`~tda.config.GridConfig` and still land on different
    edges if the sampled trajectory differed, and it is the edges that a table
    is tied to.
    """
    import hashlib

    edges = np.ascontiguousarray(grid.edges, dtype=np.float64)
    return hashlib.sha256(edges.tobytes()).hexdigest()[:length]


@dataclass(frozen=True, slots=True)
class DefectTable:
    """A memory-mapped defect table with the transports it needs.

    Attributes
    ----------
    defect:
        :math:`\\Delta\\mathbf a(m_i,N)`, shape ``(M, P, 3)``, memory-mapped
        read-only.
    node_transport:
        :math:`\\Phi(t_0,m_i)\\mathbf B`, shape ``(M, 6, 3)``.
    edge_rows:
        :math:`\\mathbf H_r\\Phi(e_j,t_0)`, shape ``(M+1, 3, 6)``.
    widths:
        Cell widths :math:`\\Delta t_i`, shape ``(M,)``.
    schema:
        Provenance and layout.

    Notes
    -----
    ``defect`` is a read-only memory map.  Nothing in the campaign writes to a
    table after it is built; a table is an input, and treating it as one is
    what makes a run reproducible from its schema.
    """

    defect: Arr
    node_transport: Arr
    edge_rows: Arr
    widths: Arr
    schema: TableSchema

    def __len__(self) -> int:
        """Number of cells, :math:`M`."""
        return self.schema.n_cells

    def index_of(self, degree: int) -> int:
        """Position of ``degree`` within the candidate axis.

        Raises
        ------
        KeyError
            If the degree is not a candidate.  Deliberately not a silent
            nearest-neighbour lookup: a policy asking for a degree the table
            was not built for is a configuration error, not a rounding one.
        """
        try:
            return self.schema.candidate_degrees.index(degree)
        except ValueError:
            raise KeyError(
                f"degree {degree} is not in the candidate set "
                f"{self.schema.candidate_degrees[:4]}..."
            ) from None

    def contributions(self, degree_of_cell: NDArray[np.int_]) -> Arr:
        """:math:`\\mathbf u_i` for a schedule, shape ``(M, 6)``.

        Parameters
        ----------
        degree_of_cell:
            Truncation degree at each cell, shape ``(M,)``.  A schedule that
            is piecewise constant on the decision grid is expanded to cells by
            the caller; this function does not know about decision intervals,
            which keeps it usable by the benchmark and the controller alike.

        Returns
        -------
        ndarray, shape (M, 6)

        Raises
        ------
        ValueError
            If the shape is wrong.
        KeyError
            If any degree is outside the candidate set.
        """
        degree_of_cell = np.asarray(degree_of_cell)
        if degree_of_cell.shape != (len(self),):
            raise ValueError(
                f"degree_of_cell must have shape ({len(self)},), got "
                f"{degree_of_cell.shape}")
        columns = np.array([self.index_of(int(n)) for n in degree_of_cell])
        picked = self.defect[np.arange(len(self)), columns]
        return np.einsum("iab,ib,i->ia", self.node_transport, picked,
                         self.widths)

    def contributions_at(self, degree: int) -> Arr:
        """:math:`\\mathbf u_i` for the constant-degree schedule, ``(M, 6)``.

        The comparator's case, and the one the multi-start descent begins
        from, so it avoids the fancy-index gather of :meth:`contributions`.
        """
        column = self.index_of(degree)
        return np.einsum("iab,ib,i->ia", self.node_transport,
                         self.defect[:, column, :], self.widths)


def build_defect_table(
    directory: Path | str,
    field: GravityField,
    grid: AccumulationGrid,
    arc_epochs: Arr,
    arc_positions: Arr,
    arc_stm: Arr,
    candidate_degrees: tuple[int, ...] | list[int],
    reference_degree: int,
    orbit: str,
    config_digest: str,
    dtype: str = "float64",
) -> DefectTable:
    """Tabulate the defect over the accumulation grid and write it to disk.

    Parameters
    ----------
    directory:
        Where ``<orbit>.npy`` and ``<orbit>.schema.json`` are written.
        Created if absent.
    field:
        Gravity field adapter; :meth:`~tda.field.GravityField.defect` is
        called ``M * P`` times and is the whole cost of this function.
    grid:
        The accumulation grid.  Its :meth:`~tda.grids.AccumulationGrid.required_epochs`
        must be exactly what the arc was propagated on.
    arc_epochs, arc_positions, arc_stm:
        The propagated arc, sampled on ``grid.required_epochs()``: epochs
        ``(2M+1,)``, inertial positions ``(2M+1, 3)`` and state-transition
        matrices :math:`\\Phi(t,t_0)` of shape ``(2M+1, 6, 6)``.  Passed as
        arrays rather than as a :class:`~tda.dynamics.ReferenceArc` so that
        this module does not depend on the integrator.
    candidate_degrees:
        :math:`\\mathcal N`.  Sorted and de-duplicated on the way in; every
        entry must be at most ``reference_degree``.
    reference_degree:
        The orbit's adopted reference degree.
    orbit:
        Identifier written into the schema and used for the filenames.
    config_digest:
        :meth:`tda.config.RunConfig.digest` of the run.
    dtype:
        Storage precision.  ``float32`` halves the file and exists for the
        parity test of D60; it is not the default and a run that uses it says
        so in its schema.

    Returns
    -------
    DefectTable
        Backed by the file just written, opened read-only.

    Raises
    ------
    ValueError
        If the arc was not sampled on the grid's required epochs, the shapes
        disagree, or a candidate exceeds the reference degree.

    Notes
    -----
    Written cell by cell into a memory map, so peak resident memory is one
    row rather than the whole table.  At degree 300 the file is about 74 MB
    per orbit and roughly two gigabytes across a twenty-six orbit panel; it is
    read sequentially and never needs to fit in memory.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)

    degrees = tuple(sorted({int(n) for n in candidate_degrees}))
    if not degrees:
        raise ValueError("candidate_degrees is empty")
    if degrees[-1] > reference_degree:
        raise ValueError(
            f"candidate degree {degrees[-1]} exceeds the reference degree "
            f"{reference_degree}; the defect is defined for truncations")

    expected = grid.required_epochs()
    arc_epochs = np.asarray(arc_epochs, dtype=float)
    if arc_epochs.shape != expected.shape or not np.allclose(
            arc_epochs, expected, rtol=0.0, atol=1e-9):
        raise ValueError(
            "the arc was not sampled on grid.required_epochs(); the defect "
            "must be evaluated where the quadrature reads it, not "
            "interpolated to it")
    arc_positions = np.asarray(arc_positions, dtype=float)
    arc_stm = np.asarray(arc_stm, dtype=float)
    if arc_positions.shape != (expected.size, 3):
        raise ValueError(
            f"arc_positions must have shape ({expected.size}, 3), got "
            f"{arc_positions.shape}")
    if arc_stm.shape != (expected.size, 6, 6):
        raise ValueError(
            f"arc_stm must have shape ({expected.size}, 6, 6), got "
            f"{arc_stm.shape}")

    node_idx = grid.node_indices()
    edge_idx = grid.edge_indices()
    node_transport = np.ascontiguousarray(
        velocity_transport(symplectic_inverse(arc_stm[node_idx])))
    edge_rows = np.ascontiguousarray(position_rows(arc_stm[edge_idx]))

    n_cells = len(grid)
    path = directory / f"{orbit}.npy"
    table = np.lib.format.open_memmap(
        path, mode="w+", dtype=np.dtype(dtype),
        shape=(n_cells, len(degrees), 3))

    node_times = arc_epochs[node_idx]
    node_positions = arc_positions[node_idx]
    syntheses = 0
    for i in range(n_cells):
        r, t = node_positions[i], float(node_times[i])
        for p, degree in enumerate(degrees):
            if degree == reference_degree:
                table[i, p] = 0.0          # exact, and one synthesis saved
            else:
                table[i, p] = field.defect(r, t, degree, reference_degree)
                syntheses += 1
    table.flush()
    del table

    schema = TableSchema(
        orbit=orbit,
        reference_degree=int(reference_degree),
        candidate_degrees=degrees,
        n_cells=n_cells,
        dtype=str(np.dtype(dtype)),
        config_digest=config_digest,
        grid_digest=grid_digest(grid),
        syntheses=syntheses,
    )
    (directory / f"{orbit}.schema.json").write_text(
        json.dumps(asdict(schema), indent=2), encoding="utf-8")

    return DefectTable(
        defect=np.load(path, mmap_mode="r"),
        node_transport=node_transport,
        edge_rows=edge_rows,
        widths=np.ascontiguousarray(grid.widths),
        schema=schema,
    )


def load_table(directory: Path | str, orbit: str, grid: AccumulationGrid,
               arc_stm: Arr, expect_config_digest: str | None = None
               ) -> DefectTable:
    """Reopen a table written by :func:`build_defect_table`.

    The transports are recomputed from ``arc_stm`` rather than stored, which
    is what lets the gradient-degree panel of Q13 reuse an expensive defect
    table with a differently generated :math:`\\Phi`.

    Parameters
    ----------
    directory, orbit:
        Where the table lives.
    grid:
        The grid the table must have been built on; its digest is checked.
    arc_stm:
        State-transition matrices on ``grid.required_epochs()``.
    expect_config_digest:
        If given, the schema must carry this digest.  Omit only when the
        transports are deliberately from a different configuration -- the Q13
        panel is the one case, and it records the substitution.

    Returns
    -------
    DefectTable

    Raises
    ------
    FileNotFoundError
        If either file is missing.
    ValueError
        On a schema-version, grid or configuration mismatch.
    """
    directory = Path(directory)
    schema_path = directory / f"{orbit}.schema.json"
    if not schema_path.exists():
        raise FileNotFoundError(f"no schema at {schema_path}")
    raw = json.loads(schema_path.read_text(encoding="utf-8"))
    raw["candidate_degrees"] = tuple(raw["candidate_degrees"])
    schema = TableSchema(**raw)

    if schema.schema_version != SCHEMA_VERSION:
        raise ValueError(
            f"table {orbit} has schema version {schema.schema_version}, "
            f"this code writes {SCHEMA_VERSION}; rebuild rather than adapt")
    if schema.grid_digest != grid_digest(grid):
        raise ValueError(
            f"table {orbit} was built on a different accumulation grid; "
            "pairing it with this one would misalign every cell")
    if (expect_config_digest is not None
            and schema.config_digest != expect_config_digest):
        raise ValueError(
            f"table {orbit} carries config digest {schema.config_digest}, "
            f"expected {expect_config_digest}")

    arc_stm = np.asarray(arc_stm, dtype=float)
    return DefectTable(
        defect=np.load(directory / f"{orbit}.npy", mmap_mode="r"),
        node_transport=np.ascontiguousarray(
            velocity_transport(symplectic_inverse(arc_stm[grid.node_indices()]))),
        edge_rows=np.ascontiguousarray(
            position_rows(arc_stm[grid.edge_indices()])),
        widths=np.ascontiguousarray(grid.widths),
        schema=schema,
    )
