"""Tests for the Chesapeake growth ingest, using a subset of the real Zenodo file."""

from pathlib import Path

import pandas as pd
import pytest

from oyster_data.sources import zenodo_growth as zg

FIX = Path(__file__).parent / "fixtures" / "zenodo11205266_sample.csv"


def sample() -> pd.DataFrame:
    return zg.read(FIX.read_text(encoding="utf-8"))


def test_growth_obs_means_match_raw_rows() -> None:
    df = sample()
    table = zg.growth_obs(df).to_pandas()
    assert set(table["cohort_id"]) == {c.cohort_id for c in zg.COHORTS}
    row = table[(table.cohort_id == "ks2009_patuxent_3n") & (table.date == pd.Timestamp("2005-12-13").date())].iloc[0]
    raw = df[(df.Waterbody_Name == "Patuxent River") & (df.Date_Oysters_Removed == "2005-12-13")]
    assert row.shell_height_mm == pytest.approx(raw.Total_Shell_Height_mm.mean())
    assert row.dry_tissue_g == pytest.approx(raw.Tissue_Dry_Weight_g.mean())
    assert row.n == len(raw) == 4
    assert (table.species == "C_virginica").all()


def test_parker_bricker_rows_ignored() -> None:
    table = zg.growth_obs(sample()).to_pandas()
    assert not table.dataset_id.str.contains("Parker").any()


def test_missing_cohort_raises() -> None:
    df = sample()
    with pytest.raises(ValueError, match="not found"):
        zg.growth_obs(df[df.Waterbody_Name != "York River"])


def test_fetch_saves_raw(tmp_path: Path) -> None:
    body = FIX.read_bytes()
    path = zg.fetch(raw_dir=tmp_path, get=lambda url: body)
    assert path.read_bytes() == body
    assert zg.load(raw_dir=tmp_path).shape[0] == sample().shape[0]
