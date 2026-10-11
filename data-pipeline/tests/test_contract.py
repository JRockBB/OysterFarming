"""Tests for oyster_data.contract using the SYNTHETIC fixtures in data/fixtures/."""

from __future__ import annotations

import datetime as dt
from pathlib import Path

import pyarrow as pa
import pytest

from oyster_data import contract
from oyster_data.contract import SPECS, ContractError, read_contract, validate_table

FIXTURES = Path(__file__).resolve().parents[2] / "data" / "fixtures"
TABLES = sorted(SPECS)


def fixture(name: str, case: str) -> Path:
    path = FIXTURES / f"synthetic_{name}_{case}.arrow"
    assert path.exists(), f"missing fixture {path}; run scripts/make_fixtures.py"
    return path


@pytest.mark.parametrize("name", TABLES)
def test_valid_fixture_accepted(name: str) -> None:
    table = read_contract(fixture(name, "valid"), name)
    assert table.num_rows > 0


@pytest.mark.parametrize("name", TABLES)
@pytest.mark.parametrize(
    "case, message",
    [("missing_column", "missing column"), ("wrong_type", "has type"), ("out_of_range", "value")],
)
def test_broken_fixture_rejected(name: str, case: str, message: str) -> None:
    with pytest.raises(ContractError, match=message):
        read_contract(fixture(name, case), name)


@pytest.mark.parametrize("name", TABLES)
def test_named_validator_exists(name: str) -> None:
    validator = getattr(contract, f"validate_{name}")
    with pa.OSFile(str(fixture(name, "valid")), "rb") as f:
        validator(pa.ipc.open_file(f).read_all())


def _forcing() -> pa.Table:
    with pa.OSFile(str(fixture("forcing_daily", "valid")), "rb") as f:
        return pa.ipc.open_file(f).read_all()


def _replace(table: pa.Table, col: str, values: list) -> pa.Table:
    idx = table.column_names.index(col)
    field = table.schema.field(col)
    return table.set_column(idx, field, pa.array(values, type=field.type))


def test_temperature_bounds_are_strict() -> None:
    t = _forcing()
    vals = t.column("temp_c").to_pylist()
    vals[0] = -2.0
    with pytest.raises(ContractError, match="temp_c"):
        validate_table(_replace(t, "temp_c", vals), SPECS["forcing_daily"])


def test_nan_rejected() -> None:
    t = _forcing()
    vals = t.column("salinity_psu").to_pylist()
    vals[3] = float("nan")
    with pytest.raises(ContractError, match="NaN"):
        validate_table(_replace(t, "salinity_psu", vals), SPECS["forcing_daily"])


def test_null_in_required_column_rejected() -> None:
    t = _forcing()
    vals = t.column("salinity_psu").to_pylist()
    vals[3] = None
    with pytest.raises(ContractError, match="not nullable"):
        validate_table(_replace(t, "salinity_psu", vals), SPECS["forcing_daily"])


def test_duplicate_site_day_rejected() -> None:
    t = _forcing()
    dates = t.column("date").to_pylist()
    dates[1] = dates[0]
    with pytest.raises(ContractError, match="duplicate"):
        validate_table(_replace(t, "date", dates), SPECS["forcing_daily"])


def test_unexpected_column_rejected() -> None:
    t = _forcing().append_column("extra", pa.array([1] * 60))
    with pytest.raises(ContractError, match="unexpected column"):
        validate_table(t, SPECS["forcing_daily"])


def test_closure_end_before_start_rejected() -> None:
    with pa.OSFile(str(fixture("closures", "valid")), "rb") as f:
        t = pa.ipc.open_file(f).read_all()
    t = _replace(t, "end_date", [dt.date(2020, 3, 1), None])
    with pytest.raises(ContractError, match="before"):
        validate_table(t, SPECS["closures"])


def test_write_contract_round_trip(tmp_path: Path) -> None:
    t = _forcing()
    out = tmp_path / "forcing_daily.arrow"
    contract.write_contract(t, out, "forcing_daily")
    assert read_contract(out, "forcing_daily").equals(t)


def test_write_contract_refuses_invalid(tmp_path: Path) -> None:
    out = tmp_path / "bad.arrow"
    with pytest.raises(ContractError):
        contract.write_contract(_forcing().drop_columns(["temp_c"]), out, "forcing_daily")
    assert not out.exists()
