"""Tests for oyster_data.gapfill on SYNTHETIC daily series."""

import numpy as np
import pandas as pd

from oyster_data import gapfill


def synthetic_daily(dates: pd.DatetimeIndex, site: str = "s", source: str = "synthetic") -> pd.DataFrame:
    doy = dates.dayofyear.to_numpy()
    return pd.DataFrame({
        "site_id": site, "date": dates.date,
        "temp_c": 10 + 5 * np.sin(2 * np.pi * doy / 365), "salinity_psu": 30.0,
        "chl_ug_l": 2.0, "pom_mg_l": None, "source": source, "qc_flag": 0,
    })


def test_short_gap_linear() -> None:
    s = pd.Series([1.0, np.nan, np.nan, 4.0], index=pd.date_range("2020-01-01", periods=4))
    vals, meth = gapfill.fill_series(s)
    assert vals.tolist() == [1.0, 2.0, 3.0, 4.0]
    assert meth.tolist() == ["observed", "linear", "linear", "observed"]


def test_long_gap_uses_climatology_and_undefined_stays_missing() -> None:
    idx = pd.date_range("2019-01-01", "2020-12-31")
    s = pd.Series(np.where(idx.month.isin([6, 7, 8]), 20.0, np.nan), index=idx)
    july_2020 = (idx.year == 2020) & (idx.month == 7)
    s[july_2020] = np.nan                          # 31-day gap inside summer 2020
    vals, meth = gapfill.fill_series(s)
    assert (meth[july_2020] == "climatology").all()
    assert np.allclose(vals[july_2020], 20.0)
    assert (meth[idx.month == 1] == "missing").all()   # no winter data in any year


def test_gap_at_edge_not_linear() -> None:
    s = pd.Series([np.nan, 2.0, 3.0], index=pd.date_range("2020-06-01", periods=3))
    _, meth = gapfill.fill_series(s)
    assert meth.iloc[0] == "climatology"


def test_climatology_wraps_new_year() -> None:
    idx = pd.date_range("2019-12-28", "2019-12-31")
    clim = gapfill.climatology(pd.Series(5.0, index=idx))
    assert clim[1] == 5.0 and clim[365] == 5.0 and np.isnan(clim[180])


def test_fill_site_flags_and_counts() -> None:
    idx = pd.date_range("2020-05-01", "2020-09-30")
    gap = pd.date_range("2020-06-10", "2020-06-12").date
    daily = synthetic_daily(idx)
    daily = daily[~daily["date"].isin(gap)]
    out, rep = gapfill.fill_site(daily, "s")
    assert len(out) == len(idx)
    filled = out["date"].isin(gap)
    assert (out.loc[filled, "qc_flag"] == 1).all()
    assert (out.loc[filled, "source"] == "synthetic;fill=linear").all()
    assert (out.loc[~filled, "qc_flag"] == 0).all()
    assert rep.counts["temp_c"]["linear"] == 3
    assert rep.longest_climatology_run["temp_c"] == 0


def test_merge_priority() -> None:
    idx = pd.date_range("2020-06-01", periods=3)
    a = synthetic_daily(idx, source="a")
    b = synthetic_daily(idx, source="b")
    b["temp_c"] = 99.0
    merged = gapfill.merge_sources([a, b])
    assert (merged["source"] == "a").all()
    assert (merged["temp_c"] < 99).all()


def test_longest_run() -> None:
    assert gapfill.longest_run(pd.Series([True, True, False, True])) == 2


def test_report_renders() -> None:
    _, rep = gapfill.fill_site(synthetic_daily(pd.date_range("2020-06-01", periods=10)), "s")
    text = gapfill.render_report([rep])
    assert "## s" in text and "| temp_c | 100.0 |" in text
