"""Tests for oyster_data.farm_ingest on a SYNTHETIC filled-in example (tests/fixtures/farm_example_synthetic)."""

import shutil
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pytest

from oyster_data import farm_ingest
from oyster_data.contract import read_contract

EXAMPLE = Path(__file__).parent / "fixtures" / "farm_example_synthetic"
TEMPLATES = Path(__file__).parents[2] / "docs" / "farm_templates"


def test_templates_match_converter_columns():
    for name, cols in farm_ingest.TEMPLATES.items():
        header = (TEMPLATES / f"{name}.csv").read_text(encoding="utf-8").splitlines()[0].split(",")
        assert header == cols


def test_example_converts_and_validates(tmp_path):
    out = farm_ingest.ingest(EXAMPLE, tmp_path)
    code = farm_ingest.site_code("Synthetic Example Farm", "North Lease")
    sites = read_contract(out["sites"], "sites").to_pylist()
    assert sites[0]["site_id"] == code and sites[0]["name"] == code
    growth = read_contract(out["growth_obs"], "growth_obs").to_pandas()
    assert list(growth["n"]) == [3, 2]                      # empty shell height not counted
    assert growth["shell_height_mm"].tolist() == [22.0, 42.0]
    assert pd.isna(growth["dry_tissue_g"].iloc[0]) and growth["dry_tissue_g"].iloc[1] == pytest.approx(0.6)
    assert set(growth["cohort_id"]) == {f"{code}_c25"}
    mort = read_contract(out["mortality_obs"], "mortality_obs").to_pandas()
    assert mort["mortality_fraction"].tolist() == [0.05, 0.1]
    assert mort["period_days"].iloc[0] == 61 and pd.isna(mort["period_days"].iloc[1])
    temp = pd.read_csv(out["farm_temperature_daily"])
    assert temp["date"].tolist() == ["2025-07-01"]            # 2025-07-02 has only 3 hours
    for f in out.values():                                     # names never written
        body = Path(f).read_bytes()
        assert b"Synthetic Example Farm" not in body and b"North Lease" not in body


def test_bad_inputs_raise(tmp_path):
    folder = tmp_path / "farm"
    shutil.copytree(EXAMPLE, folder)
    (folder / "growth_checks.csv").write_text("date,cohort,oyster_no,shell_height_mm,dry_tissue_g\n"
                                              "06/01/2025,c25,1,20,\n", encoding="utf-8")
    with pytest.raises(farm_ingest.FarmDataError, match="YEAR-MONTH-DAY"):
        farm_ingest.ingest(folder, tmp_path / "out")
    (folder / "growth_checks.csv").write_text("date,cohort,height\n2025-06-01,c25,20\n", encoding="utf-8")
    with pytest.raises(farm_ingest.FarmDataError, match="columns must be"):
        farm_ingest.ingest(folder, tmp_path / "out")
    (folder / "growth_checks.csv").unlink()
    with pytest.raises(farm_ingest.FarmDataError, match="missing required"):
        farm_ingest.ingest(folder, tmp_path / "out")
