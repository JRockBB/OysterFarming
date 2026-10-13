"""Tests for the Chesapeake forcing sources, using saved real responses."""

import json
from pathlib import Path

import pandas as pd
import pytest

from oyster_data.forcing import interpolate_samples
from oyster_data.sources import cbp_datahub, coops, nerrs_ncei

FIX = Path(__file__).parent / "fixtures"


def test_coops_daily_mean_and_min_obs() -> None:
    s = coops.parse((FIX / "coops_8577330_2006_head.csv").read_text())
    # 59 hourly rows: 2006-01-01 01:00 .. 2006-01-03 11:00 -> two full days, third day has 12 hours
    assert list(s.index) == [pd.Timestamp("2006-01-01"), pd.Timestamp("2006-01-02"), pd.Timestamp("2006-01-03")]
    assert 0 < s.min() < s.max() < 10


def test_coops_drops_flagged_and_rejects_bad_response(tmp_path: Path) -> None:
    text = "Date Time, Water Temperature, X, N, R\n" + "".join(
        f"2006-01-01 {h:02d}:00,{'40.0' if h == 0 else '6.0'},{1 if h == 0 else 0},0,0\n" for h in range(24))
    assert coops.parse(text).iloc[0] == pytest.approx(6.0)
    with pytest.raises(RuntimeError):
        coops.fetch_year("8577330", 2006, raw_dir=tmp_path, get=lambda url: b"Error: no data")


def test_datahub_upper_layers_only() -> None:
    """LE1.3 on 2006-07-10 reports S, AP, BP, B and M layers; only S and AP are used."""
    g = cbp_datahub.parse(json.loads((FIX / "cbp_LE1.3_sample.json").read_text()))
    row = g.loc[pd.Timestamp("2006-07-10")]
    assert row.chl_ug_l == pytest.approx((23.92 + 20.04) / 2)
    assert row.salinity_psu == pytest.approx((10.69 + 10.76) / 2)
    assert row.temp_c == pytest.approx(26.1)


def test_datahub_drops_problem_and_qualifier() -> None:
    recs = [{"SampleDate": "2005-11-07T00:00:00", "layer": "S", "parameter": "CHLA", "measureValue": None,
             "qualifier": "", "problem": "A"},
            {"SampleDate": "2005-11-08T00:00:00", "layer": "S", "parameter": "CHLA", "measureValue": 0.5,
             "qualifier": "<", "problem": ""},
            {"SampleDate": "2005-11-09T00:00:00", "layer": "S", "parameter": "CHLA", "measureValue": 4.0,
             "qualifier": "", "problem": ""}]
    g = cbp_datahub.parse(recs)
    assert list(g["chl_ug_l"].dropna().index) == [pd.Timestamp("2005-11-09")]


def test_nerrs_sonde_flags_and_daily_mean() -> None:
    text = (FIX / "nerrs_cbvgiwq2006_head.csv").read_text()
    t = nerrs_ncei.parse(text, "Temp")
    s = nerrs_ncei.parse(text, "Sal")
    assert pd.Timestamp("2006-01-01") in t.index and pd.Timestamp("2006-01-01") in s.index
    assert 0 < t.iloc[0] < 15 and 15 < s.iloc[0] < 25
    rejected = text.replace("<4> ", "<-3> ")
    assert nerrs_ncei.parse(rejected, "Temp").empty


def test_nerrs_nutrient_chl_replicates_averaged() -> None:
    s = nerrs_ncei.parse_nutrient_chl((FIX / "nerrs_cbvginut2005_novdec.csv").read_text())
    assert s[pd.Timestamp("2005-12-13")] == pytest.approx((9.4 + 8.5) / 2)
    assert s[pd.Timestamp("2005-11-11")] == pytest.approx((8.8 + 8.1) / 2)


def test_interpolate_samples_respects_max_gap() -> None:
    s = pd.Series([1.0, 3.0, 10.0], index=pd.to_datetime(["2006-01-01", "2006-01-11", "2006-04-01"]))
    d = interpolate_samples(s, 62)
    assert d[pd.Timestamp("2006-01-06")] == pytest.approx(2.0)
    assert d[pd.Timestamp("2006-02-15")] != d[pd.Timestamp("2006-02-15")]   # 80-day gap stays NaN
