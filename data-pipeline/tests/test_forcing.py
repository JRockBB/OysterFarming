"""Tests for oyster_data.forcing. Inputs here are SYNTHETIC hourly records."""

import numpy as np
import pandas as pd

from oyster_data.forcing import hourly_to_daily, screen, to_forcing_table

DEFAULTS = {"temp_c": 15.0, "salinity_psu": 30.0, "chl_ug_l": 2.0}


def synthetic_hourly(n_hours: int, start: str = "2020-06-01") -> pd.DataFrame:
    data = {"time_utc": pd.date_range(start, periods=n_hours, freq="h", tz="UTC")}
    for k, v in DEFAULTS.items():
        data[k] = np.full(n_hours, v)
    return pd.DataFrame(data)


def test_screen_removes_implausible() -> None:
    df = synthetic_hourly(4)
    df.loc[0, "temp_c"] = 40.0
    df.loc[1, "chl_ug_l"] = -0.04
    df.loc[2, "salinity_psu"] = 50.0
    out, removed = screen(df)
    assert removed == {"temp_c": 1, "salinity_psu": 1, "chl_ug_l": 1}
    assert out.isna().sum().sum() == 3


def test_daily_mean_and_min_obs() -> None:
    df = synthetic_hourly(48)
    df.loc[24:47, "temp_c"] = 17.0          # day 2 has a different temperature
    df.loc[24:40, "chl_ug_l"] = np.nan      # day 2 keeps only 7 chl values -> null
    daily = hourly_to_daily(df, "s", "synthetic")
    assert list(daily["temp_c"]) == [15.0, 17.0]
    assert daily["chl_ug_l"].isna().tolist() == [False, True]
    assert (daily["qc_flag"] == 0).all()


def test_days_without_temperature_dropped() -> None:
    df = synthetic_hourly(48)
    df.loc[0:20, "temp_c"] = np.nan
    daily = hourly_to_daily(df, "s", "synthetic")
    assert len(daily) == 1


def test_to_forcing_table_validates() -> None:
    daily = hourly_to_daily(synthetic_hourly(72), "s", "synthetic")
    assert to_forcing_table(daily).num_rows == 3
