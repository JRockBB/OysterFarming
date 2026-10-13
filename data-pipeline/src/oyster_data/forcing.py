"""Shared steps for building ``forcing_daily`` rows from sub-daily sensor records.

These are data-processing choices, not model parameters:

- Values outside physically plausible ranges are treated as sensor errors and
  dropped before averaging (same bounds as the data contract; chlorophyll < 0).
- A daily mean is kept only if at least ``MIN_OBS_PER_DAY`` valid values exist
  for that day; otherwise that variable is null for the day.
- Days without both temperature and salinity are omitted, because those columns
  are not nullable in the contract. Gap filling (PLAN.md P1.4) restores them.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pyarrow as pa

from oyster_data.contract import SPECS, arrow_schema, validate_table

MIN_OBS_PER_DAY = 12

PLAUSIBLE = {
    "temp_c": (-2.0, 35.0, True),   # (lo, hi, open interval)
    "salinity_psu": (0.0, 45.0, False),
    "chl_ug_l": (0.0, np.inf, False),
}

VARIABLES = ("temp_c", "salinity_psu", "chl_ug_l")


def screen(df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    """Set implausible values to NaN. Returns the screened copy and counts removed per column."""
    out = df.copy()
    removed: dict[str, int] = {}
    for col, (lo, hi, open_) in PLAUSIBLE.items():
        if col not in out:
            continue
        v = out[col]
        bad = (v <= lo) | (v >= hi) if open_ else (v < lo) | (v > hi)
        bad &= v.notna()
        removed[col] = int(bad.sum())
        out.loc[bad, col] = np.nan
    return out, removed


def hourly_to_daily(df: pd.DataFrame, site_id: str, source: str,
                    min_obs: int = MIN_OBS_PER_DAY) -> pd.DataFrame:
    """Average sub-daily records to UTC calendar days.

    ``df`` needs a tz-aware ``time_utc`` column plus any of ``temp_c``,
    ``salinity_psu``, ``chl_ug_l``. Returns contract-shaped daily rows.
    """
    screened, _ = screen(df)
    screened = screened.assign(date=screened["time_utc"].dt.tz_convert("UTC").dt.date)
    rows = []
    for day, g in screened.groupby("date", sort=True):
        rec = {"site_id": site_id, "date": day}
        for col in VARIABLES:
            vals = g[col].dropna() if col in g else pd.Series(dtype=float)
            rec[col] = float(vals.mean()) if len(vals) >= min_obs else None
        rows.append(rec)
    daily = pd.DataFrame(rows, columns=["site_id", "date", *VARIABLES])
    daily = daily.dropna(subset=["temp_c", "salinity_psu"]).reset_index(drop=True)
    daily["pom_mg_l"] = None
    daily["source"] = source
    daily["qc_flag"] = 0
    return daily


def interpolate_samples(samples: pd.Series, max_gap_days: int) -> pd.Series:
    """Daily series from discrete samples: linear interpolation between samples at
    most ``max_gap_days`` apart; longer gaps stay missing. Sample days keep their value."""
    if samples.empty:
        return samples
    full = pd.date_range(samples.index.min(), samples.index.max(), freq="D")
    daily = samples.reindex(full)
    known = daily.notna()
    interp = daily.interpolate(method="time", limit_area="inside")
    prev_obs = pd.Series(full.where(known), index=full).ffill()
    next_obs = pd.Series(full.where(known), index=full).bfill()
    gap = (next_obs - prev_obs).dt.days
    ok = known | (gap <= max_gap_days)
    out = interp.where(ok)
    out.name = samples.name
    return out


def to_forcing_table(daily: pd.DataFrame) -> pa.Table:
    """Convert daily rows to a validated ``forcing_daily`` Arrow table."""
    schema = arrow_schema("forcing_daily")
    cols = {}
    for field in schema:
        values = daily[field.name].tolist()
        values = [None if (isinstance(v, float) and np.isnan(v)) else v for v in values]
        cols[field.name] = pa.array(values, type=field.type)
    table = pa.Table.from_pydict(cols, schema=schema)
    validate_table(table, SPECS["forcing_daily"])
    return table
