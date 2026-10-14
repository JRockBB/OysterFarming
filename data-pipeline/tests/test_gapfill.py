"""Tests for oyster_data.gapfill on SYNTHETIC daily series."""

import numpy as np
import pandas as pd
import pytest

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
    vals, meth, _ = gapfill.fill_series(s)
    assert vals.tolist() == [1.0, 2.0, 3.0, 4.0]
    assert meth.tolist() == ["observed", "linear", "linear", "observed"]


def test_long_gap_uses_climatology_and_undefined_stays_missing() -> None:
    idx = pd.date_range("2019-01-01", "2020-12-31")
    s = pd.Series(np.where(idx.month.isin([6, 7, 8]), 20.0, np.nan), index=idx)
    july_2020 = (idx.year == 2020) & (idx.month == 7)
    s[july_2020] = np.nan                          # 31-day gap inside summer 2020
    vals, meth, _ = gapfill.fill_series(s)
    assert (meth[july_2020] == "climatology").all()
    assert np.allclose(vals[july_2020], 20.0)
    assert (meth[idx.month == 1] == "missing").all()   # no winter data in any year


def test_gap_at_edge_not_linear() -> None:
    s = pd.Series([np.nan, 2.0, 3.0], index=pd.date_range("2020-06-01", periods=3))
    _, meth, _ = gapfill.fill_series(s)
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
    assert "## s" in text and "| temp_c | observed | 10 | 100.0 |" in text


def test_regression_fill_from_secondary() -> None:
    """SYNTHETIC: site = 1 + 2 * station exactly; a 40-day gap is filled by the fit."""
    idx = pd.date_range("2020-01-01", periods=200)
    station = pd.Series(np.sin(np.arange(200) / 10.0) + 5, index=idx)
    site = 1 + 2 * station
    site.iloc[100:140] = np.nan
    vals, meth, fits = gapfill.fill_series(site, secondaries=[("stationA", station)])
    assert (meth.iloc[100:140] == "regression:stationA").all()
    # filled values are clipped to the observed range; the gap trough is ~1e-4 below the observed minimum
    assert np.allclose(vals.iloc[100:140], 1 + 2 * station.iloc[100:140], atol=1e-3)
    assert len(fits) == 1 and fits[0].r2 > 0.999 and fits[0].n_filled == 40
    assert fits[0].slope == pytest.approx(2) and fits[0].intercept == pytest.approx(1)


def test_regression_needs_overlap_and_falls_back() -> None:
    """SYNTHETIC: a station with too little overlap is skipped; the next one is used."""
    idx = pd.date_range("2020-01-01", periods=120)
    site = pd.Series(np.linspace(10, 20, 120), index=idx)
    site.iloc[50:90] = np.nan
    short = pd.Series(1.0, index=idx[:10])
    good = pd.Series(np.linspace(0, 1, 120), index=idx)
    vals, meth, fits = gapfill.fill_series(site, secondaries=[("short", short), ("good", good)])
    assert [f.station for f in fits] == ["good"]
    assert (meth.iloc[50:90] == "regression:good").all()


def test_fill_site_source_lists_methods() -> None:
    idx = pd.date_range("2020-05-01", "2020-09-30")
    daily = synthetic_daily(idx)
    drop = pd.date_range("2020-06-01", "2020-07-15").date
    daily = daily[~daily["date"].isin(drop)]
    station = pd.Series(10 + 5 * np.sin(2 * np.pi * idx.dayofyear.to_numpy() / 365), index=idx)
    out, rep = gapfill.fill_site(daily, "s", {"temp_c": [("st", station)]})
    filled = out[out["date"].isin(drop)]
    assert (filled["qc_flag"] == 1).all()
    assert filled["source"].str.contains("regression:st").all()
    assert rep.fits["temp_c"][0].n_overlap == len(idx) - len(drop)


def test_final_pass_bridges_short_gap_between_filled_values() -> None:
    """SYNTHETIC: a 5-day hole in the secondary station inside a long site gap."""
    idx = pd.date_range("2020-01-01", periods=200)
    station = pd.Series(np.linspace(0, 10, 200), index=idx)
    site = 2 * station
    site.iloc[50:150] = np.nan
    station.iloc[90:95] = np.nan
    vals, meth, _ = gapfill.fill_series(site, secondaries=[("st", station)], climatology_fill=False)
    assert (meth.iloc[90:95] == "linear").all()
    assert vals.iloc[50:150].notna().all()


def test_regression_fill_is_clipped_to_observed_range():
    # SYNTHETIC: station runs colder than anything observed at the site
    import numpy as np
    import pandas as pd
    from oyster_data import gapfill
    idx = pd.date_range("2020-01-01", periods=100, freq="D")
    station = pd.Series(np.linspace(-5.0, 15.0, 100), index=idx)
    site = station * 0.5 + 2.0
    site.iloc[:20] = np.nan                       # gap where the station is coldest
    filled, method, fits = gapfill.fill_series(site, secondaries=[("st", station)])
    assert fits and filled.iloc[:20].min() >= site.min()
