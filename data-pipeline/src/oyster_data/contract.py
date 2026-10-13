"""Data contract validators.

The schemas below mirror docs/PLAN.md ("Data contract") exactly. The Julia
package has matching validators in src/contract.jl. Changing a schema requires
updating PLAN.md, the writer, the Julia reader and both validators together.

Column kinds and their Arrow types:

- ``string``  -> utf8 or large_utf8
- ``float64`` -> float64
- ``int8``    -> int8
- ``int``     -> int64
- ``date``    -> date32
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.ipc as ipc


class ContractError(ValueError):
    """Raised when a table does not follow the data contract."""

    def __init__(self, table_name: str, problems: list[str]) -> None:
        self.table_name = table_name
        self.problems = problems
        lines = "\n".join(f"  - {p}" for p in problems)
        super().__init__(f"{table_name} violates the data contract:\n{lines}")


@dataclass(frozen=True)
class Column:
    """One column of a contract table.

    ``lo``/``hi`` are numeric bounds; ``lo_open``/``hi_open`` make them strict.
    ``allowed`` restricts values to a fixed set (enums).
    """

    name: str
    kind: str
    nullable: bool = False
    lo: float | None = None
    hi: float | None = None
    lo_open: bool = False
    hi_open: bool = False
    allowed: tuple[str | int, ...] | None = None


@dataclass(frozen=True)
class TableSpec:
    name: str
    columns: tuple[Column, ...]
    unique_key: tuple[str, ...] = ()
    ordered_dates: tuple[str, str] | None = None  # (start, end): end >= start when both present


CLOSURE_TYPES = ("rainfall", "biotoxin", "wastewater", "pollution", "conservation", "other")

SPECS: dict[str, TableSpec] = {
    "sites": TableSpec(
        "sites",
        (
            Column("site_id", "string"),
            Column("name", "string"),
            Column("lat", "float64", lo=-90, hi=90),
            Column("lon", "float64", lo=-180, hi=180),
            Column("region", "string", allowed=("maine", "chesapeake")),
            Column("growing_area", "string", nullable=True),
            Column("notes", "string", nullable=True),
        ),
        unique_key=("site_id",),
    ),
    "forcing_daily": TableSpec(
        "forcing_daily",
        (
            Column("site_id", "string"),
            Column("date", "date"),
            Column("temp_c", "float64", lo=-2, hi=35, lo_open=True, hi_open=True),
            Column("salinity_psu", "float64", lo=0, hi=45),
            Column("chl_ug_l", "float64", nullable=True, lo=0),
            Column("pom_mg_l", "float64", nullable=True, lo=0),
            Column("source", "string"),
            Column("qc_flag", "int8", allowed=(0, 1, 2)),
        ),
        unique_key=("site_id", "date"),
    ),
    "closures": TableSpec(
        "closures",
        (
            Column("growing_area", "string"),
            Column("area_id", "string", nullable=True),
            Column("start_date", "date"),
            Column("end_date", "date", nullable=True),
            Column("closure_type", "string", allowed=CLOSURE_TYPES),
            Column("source", "string"),
        ),
        ordered_dates=("start_date", "end_date"),
    ),
    "growth_obs": TableSpec(
        "growth_obs",
        (
            Column("dataset_id", "string"),
            Column("species", "string", allowed=("C_virginica", "C_gigas")),
            Column("site_id", "string"),
            Column("cohort_id", "string"),
            Column("date", "date"),
            Column("shell_height_mm", "float64", nullable=True, lo=0),
            Column("dry_tissue_g", "float64", nullable=True, lo=0),
            Column("n", "int", lo=0),
        ),
    ),
    "disease_obs": TableSpec(
        "disease_obs",
        (
            Column("site_id", "string"),
            Column("date", "date"),
            Column("pathogen", "string", allowed=("P_marinus", "H_nelsoni")),
            Column("prevalence", "float64", lo=0, hi=1),
            Column("intensity_mean", "float64", nullable=True, lo=0),
            Column("intensity_scale", "string"),
            Column("n", "int", nullable=True, lo=0),
            Column("source", "string"),
        ),
    ),
    "economics": TableSpec(
        "economics",
        (
            Column("year", "int", lo=1800, hi=2200),
            Column("region", "string"),
            Column("value_usd", "float64", lo=0),
            Column("pounds", "float64", nullable=True, lo=0),
            Column("count", "float64", nullable=True, lo=0),
            Column("source", "string"),
        ),
    ),
    "mortality_obs": TableSpec(
        "mortality_obs",
        (
            Column("site_id", "string"),
            Column("date", "date"),
            Column("period_days", "int", nullable=True, lo=1),
            Column("mortality_fraction", "float64", lo=0, hi=1),
            Column("method", "string"),
            Column("source", "string"),
        ),
    ),
}


def _type_ok(kind: str, typ: pa.DataType) -> bool:
    if kind == "string":
        return pa.types.is_string(typ) or pa.types.is_large_string(typ)
    if kind == "float64":
        return typ == pa.float64()
    if kind == "int8":
        return typ == pa.int8()
    if kind == "int":
        return typ == pa.int64()
    if kind == "date":
        return typ == pa.date32()
    raise ValueError(f"unknown column kind {kind!r}")


def _check_column(col: Column, data: pa.ChunkedArray) -> list[str]:
    problems: list[str] = []
    if not _type_ok(col.kind, data.type):
        return [f"column {col.name!r} has type {data.type}, expected {col.kind}"]
    if not col.nullable and data.null_count > 0:
        problems.append(f"column {col.name!r} has {data.null_count} null(s) but is not nullable")
    if col.kind == "float64" and len(data) > 0:
        nan_count = pc.sum(pc.is_nan(data)).as_py() or 0
        if nan_count:
            problems.append(f"column {col.name!r} has {nan_count} NaN(s); use nulls for missing values")
    if col.lo is not None:
        below = pc.less_equal(data, col.lo) if col.lo_open else pc.less(data, col.lo)
        n = pc.sum(below).as_py() or 0
        if n:
            op = ">" if col.lo_open else ">="
            problems.append(f"column {col.name!r} has {n} value(s) not {op} {col.lo}")
    if col.hi is not None:
        above = pc.greater_equal(data, col.hi) if col.hi_open else pc.greater(data, col.hi)
        n = pc.sum(above).as_py() or 0
        if n:
            op = "<" if col.hi_open else "<="
            problems.append(f"column {col.name!r} has {n} value(s) not {op} {col.hi}")
    if col.allowed is not None:
        bad = {v for v in data.to_pylist() if v is not None and v not in col.allowed}
        if bad:
            problems.append(f"column {col.name!r} has values outside {list(col.allowed)}: {sorted(map(str, bad))}")
    return problems


def validate_table(table: pa.Table, spec: TableSpec) -> None:
    """Validate ``table`` against ``spec``; raise ContractError listing every problem."""
    problems: list[str] = []
    expected = [c.name for c in spec.columns]
    present = set(table.column_names)
    for name in expected:
        if name not in present:
            problems.append(f"missing column {name!r}")
    for name in table.column_names:
        if name not in expected:
            problems.append(f"unexpected column {name!r}")
    for col in spec.columns:
        if col.name in present:
            problems.extend(_check_column(col, table.column(col.name)))

    if not problems and spec.unique_key:
        keys = list(zip(*(table.column(k).to_pylist() for k in spec.unique_key)))
        n_dup = len(keys) - len(set(keys))
        if n_dup:
            problems.append(f"{n_dup} duplicate row(s) for key {list(spec.unique_key)}")
    if not problems and spec.ordered_dates:
        start, end = spec.ordered_dates
        n_bad = pc.sum(pc.less(table.column(end), table.column(start))).as_py() or 0
        if n_bad:
            problems.append(f"{n_bad} row(s) with {end} before {start}")

    if problems:
        raise ContractError(spec.name, problems)


def validate_sites(table: pa.Table) -> None:
    validate_table(table, SPECS["sites"])


def validate_forcing_daily(table: pa.Table) -> None:
    validate_table(table, SPECS["forcing_daily"])


def validate_closures(table: pa.Table) -> None:
    validate_table(table, SPECS["closures"])


def validate_growth_obs(table: pa.Table) -> None:
    validate_table(table, SPECS["growth_obs"])


def validate_disease_obs(table: pa.Table) -> None:
    validate_table(table, SPECS["disease_obs"])


def validate_economics(table: pa.Table) -> None:
    validate_table(table, SPECS["economics"])


def validate_mortality_obs(table: pa.Table) -> None:
    validate_table(table, SPECS["mortality_obs"])


def arrow_schema(table_name: str) -> pa.Schema:
    """The canonical Arrow schema for a contract table (strings as utf8)."""
    kinds = {"string": pa.string(), "float64": pa.float64(), "int8": pa.int8(), "int": pa.int64(), "date": pa.date32()}
    spec = SPECS[table_name]
    return pa.schema([pa.field(c.name, kinds[c.kind], nullable=c.nullable) for c in spec.columns])


def read_contract(path: str | Path, table_name: str) -> pa.Table:
    """Read an Arrow IPC file and validate it as ``table_name``."""
    with pa.OSFile(str(path), "rb") as source:
        table = ipc.open_file(source).read_all()
    validate_table(table, SPECS[table_name])
    return table


def write_contract(table: pa.Table, path: str | Path, table_name: str) -> None:
    """Validate ``table`` as ``table_name`` and write it as Arrow IPC."""
    validate_table(table, SPECS[table_name])
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with pa.OSFile(str(path), "wb") as sink, ipc.new_file(sink, table.schema) as writer:
        writer.write_table(table)
