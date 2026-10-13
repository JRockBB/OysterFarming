"""Tests for the Maine secondary forcing sources, using saved real responses."""

import json
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import requests

from oyster_data.sources import boothbay_dmr, dmc_dock, neracoos, viirs_chl

FIX = Path(__file__).parent / "fixtures"


def test_boothbay_parse_local_dates() -> None:
    s = boothbay_dmr.parse(json.loads((FIX / "boothbay_sst_20180101_20180110.json").read_text()))
    assert s.index[0] == pd.Timestamp("2018-01-01")      # epoch ms at local midnight -> local date
    assert len(s) == 10 and s.index.is_monotonic_increasing
    assert s.iloc[0] == pytest.approx(1.816, abs=1e-3)


def test_boothbay_fetch_pages(tmp_path: Path) -> None:
    page = json.loads((FIX / "boothbay_sst_20180101_20180110.json").read_text())
    first = dict(page, exceededTransferLimit=True)
    calls: list[str] = []

    def fake_get(url: str) -> bytes:
        calls.append(url)
        return json.dumps(first if len(calls) == 1 else page).encode()

    path = boothbay_dmr.fetch(date(2018, 1, 1), date(2018, 1, 11), raw_dir=tmp_path, get=fake_get)
    assert len(calls) == 2 and "resultOffset=2000" in calls[1]
    assert len(json.loads(path.read_text())["features"]) == 20
    with pytest.raises(FileExistsError):
        boothbay_dmr.fetch(date(2018, 1, 1), date(2018, 1, 11), raw_dir=tmp_path, get=fake_get)


def test_dmc_dock_parse_and_daily() -> None:
    samples = dmc_dock.parse((FIX / "dmc_dock_sample.csv").read_text())
    assert samples.index.is_unique                       # 2002-08-01 has two samples, averaged
    assert samples[pd.Timestamp("2002-08-01")] == pytest.approx((1.07 + 3.62) / 2)
    daily = dmc_dock.to_daily(samples)
    # 2002-07-11 -> 2002-07-15 is a 4-day gap: interpolated
    assert daily[pd.Timestamp("2002-07-13")] == pytest.approx(1.67 + (3.56 - 1.67) / 2)


def test_dmc_dock_long_gap_left_missing() -> None:
    s = pd.Series([1.0, 3.0], index=pd.to_datetime(["2020-01-01", "2020-01-20"]))
    daily = dmc_dock.to_daily(s)
    assert daily.iloc[1:-1].isna().all() and daily.iloc[0] == 1.0 and daily.iloc[-1] == 3.0


def test_dmc_dock_fetch_rejects_html(tmp_path: Path) -> None:
    with pytest.raises(RuntimeError):
        dmc_dock.fetch(raw_dir=tmp_path, get=lambda url: b"<html>")


@pytest.mark.parametrize("variable, name, lo, hi", [("salinity", "salinity_psu", 32, 33.5),
                                                      ("temperature", "temp_c", 4, 7)])
def test_e01_parse(variable: str, name: str, lo: float, hi: float) -> None:
    s = neracoos.parse((FIX / "e01_001m_20190115_16.csv").read_text(), variable)
    assert s.name == name
    assert list(s.index) == [pd.Timestamp("2019-01-15"), pd.Timestamp("2019-01-16")]
    assert ((s > lo) & (s < hi)).all()


def test_e01_qc_and_min_obs() -> None:
    text = (FIX / "e01_001m_20190115_16.csv").read_text().splitlines()
    header, units, rows = text[0], text[1], text[2:]
    bad = [r.replace(",1,", ",4,", 1) for r in rows[:48]]       # fail QC on all of day 1
    s = neracoos.parse("\n".join([header, units, *bad, *rows[48:]]), "salinity")
    assert list(s.index) == [pd.Timestamp("2019-01-16")]


def test_e01_404_means_no_data(tmp_path: Path) -> None:
    def not_found(url: str) -> bytes:
        resp = requests.Response()
        resp.status_code = 404
        raise requests.HTTPError(response=resp)

    assert neracoos.fetch_year(2014, 10, raw_dir=tmp_path, get=not_found) is None


def test_viirs_parse_spreads_composite() -> None:
    s = viirs_chl.parse((FIX / "viirs_chla8day_20190117.csv").read_text())
    # composite stamped 2019-01-21 covers 2019-01-17 .. 2019-01-24
    assert s.index.min() == pd.Timestamp("2019-01-17") and s.index.max() == pd.Timestamp("2019-01-24")
    assert s.nunique() == 1 and np.isfinite(s.iloc[0]) and s.iloc[0] > 0
