"""Write SYNTHETIC data contract fixtures to data/fixtures/.

These files are made-up test data, not observations. They exist only so that
the Python and Julia validators can be tested against the same files. Every
file name starts with ``synthetic_``.

For each contract table this writes four files:

- ``synthetic_<table>_valid.arrow``           passes validation
- ``synthetic_<table>_missing_column.arrow``  one required column dropped
- ``synthetic_<table>_wrong_type.arrow``      one column with the wrong Arrow type
- ``synthetic_<table>_out_of_range.arrow``    one value outside its allowed range

Run from data-pipeline/:  uv run python scripts/make_fixtures.py
"""

from __future__ import annotations

import datetime as dt
from pathlib import Path

import pyarrow as pa
import pyarrow.ipc as ipc

from oyster_data.contract import arrow_schema, validate_table, SPECS

FIXTURE_DIR = Path(__file__).resolve().parents[2] / "data" / "fixtures"

D = dt.date


def _forcing_rows() -> dict[str, list]:
    """Synthetic: 60 days at one site with a smooth seasonal-looking ramp."""
    n = 60
    dates = [D(2020, 5, 1) + dt.timedelta(days=i) for i in range(n)]
    return {
        "site_id": ["synthetic_site"] * n,
        "date": dates,
        "temp_c": [8.0 + 0.2 * i for i in range(n)],
        "salinity_psu": [30.0] * n,
        "chl_ug_l": [None if i == 10 else 2.0 + 0.05 * i for i in range(n)],
        "pom_mg_l": [None] * n,
        "source": ["synthetic"] * n,
        "qc_flag": [0] * n,
    }


VALID: dict[str, dict[str, list]] = {
    "sites": {
        "site_id": ["synthetic_site", "synthetic_bay"],
        "name": ["Synthetic site A", "Synthetic site B"],
        "lat": [44.0, 38.5],
        "lon": [-69.5, -76.3],
        "region": ["maine", "chesapeake"],
        "growing_area": ["SYN1", None],
        "notes": ["synthetic fixture", None],
    },
    "forcing_daily": _forcing_rows(),
    "closures": {
        "growing_area": ["SYN1", "SYN1"],
        "area_id": [None, "SYN1-a"],
        "start_date": [D(2020, 4, 1), D(2020, 7, 10)],
        "end_date": [D(2020, 4, 5), None],
        "closure_type": ["rainfall", "biotoxin"],
        "source": ["synthetic", "synthetic"],
    },
    "growth_obs": {
        "dataset_id": ["synthetic"] * 3,
        "species": ["C_virginica"] * 3,
        "site_id": ["synthetic_site"] * 3,
        "cohort_id": ["c1"] * 3,
        "date": [D(2020, 5, 1), D(2020, 6, 1), D(2020, 7, 1)],
        "shell_height_mm": [20.0, 24.0, None],
        "dry_tissue_g": [None, 0.1, 0.15],
        "n": [30, 30, 25],
    },
    "disease_obs": {
        "site_id": ["synthetic_bay", "synthetic_bay"],
        "date": [D(2019, 10, 1), D(2019, 10, 1)],
        "pathogen": ["P_marinus", "H_nelsoni"],
        "prevalence": [0.6, 0.1],
        "intensity_mean": [2.1, None],
        "intensity_scale": ["mackin_0_5", "mackin_0_5"],
        "n": [30, None],
        "source": ["synthetic", "synthetic"],
    },
    "economics": {
        "year": [2019, 2020],
        "region": ["maine", "maine"],
        "value_usd": [1000.0, 1200.0],
        "pounds": [None, 500.0],
        "count": [4000.0, None],
        "source": ["synthetic", "synthetic"],
    },
}

# (column to drop, column to retype -> new type, (column, row, out-of-range value))
BREAKAGE: dict[str, tuple[str, tuple[str, pa.DataType], tuple[str, int, object]]] = {
    "sites": ("region", ("lat", pa.string()), ("lat", 0, 95.0)),
    "forcing_daily": ("temp_c", ("qc_flag", pa.int64()), ("temp_c", 5, 40.0)),
    "closures": ("closure_type", ("start_date", pa.string()), ("closure_type", 0, "hurricane")),
    "growth_obs": ("n", ("shell_height_mm", pa.string()), ("shell_height_mm", 0, -1.0)),
    "disease_obs": ("prevalence", ("prevalence", pa.string()), ("prevalence", 0, 1.5)),
    "economics": ("value_usd", ("year", pa.float64()), ("value_usd", 0, -5.0)),
}


def _write(table: pa.Table, path: Path) -> None:
    with pa.OSFile(str(path), "wb") as sink, ipc.new_file(sink, table.schema) as writer:
        writer.write_table(table)


def build_valid(name: str) -> pa.Table:
    return pa.Table.from_pydict(VALID[name], schema=arrow_schema(name))


def main() -> None:
    FIXTURE_DIR.mkdir(parents=True, exist_ok=True)
    for name in SPECS:
        valid = build_valid(name)
        validate_table(valid, SPECS[name])
        _write(valid, FIXTURE_DIR / f"synthetic_{name}_valid.arrow")

        drop, (retype_col, new_type), (oor_col, oor_row, oor_val) = BREAKAGE[name]
        _write(valid.drop_columns([drop]), FIXTURE_DIR / f"synthetic_{name}_missing_column.arrow")

        idx = valid.column_names.index(retype_col)
        retyped = valid.set_column(idx, pa.field(retype_col, new_type), valid.column(retype_col).cast(new_type))
        _write(retyped, FIXTURE_DIR / f"synthetic_{name}_wrong_type.arrow")

        values = valid.column(oor_col).to_pylist()
        values[oor_row] = oor_val
        idx = valid.column_names.index(oor_col)
        field = valid.schema.field(oor_col)
        out_of_range = valid.set_column(idx, field, pa.array(values, type=field.type))
        _write(out_of_range, FIXTURE_DIR / f"synthetic_{name}_out_of_range.arrow")
    print(f"wrote fixtures to {FIXTURE_DIR}")


if __name__ == "__main__":
    main()
