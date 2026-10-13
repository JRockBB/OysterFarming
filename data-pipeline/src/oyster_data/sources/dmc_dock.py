"""Darling Marine Center dock chlorophyll (secondary chlorophyll source).

Brady lab, University of Maine, "DMC Dock Sampling"
(https://umaine.edu/bradylab/dmc-dock-sampling/): surface water from the
Darling Marine Center dock on the Damariscotta River, sampled 2-5 mornings a
week since 2002; extracted chlorophyll-a and phaeopigments (ug/L). Data file:
https://www.dropbox.com/s/0oru9p1b7v5pc1r/data_dock.csv?dl=1
Columns: Date, Ave. Chl, Ave. Phaeo, Temp (_C). The temperature column is in
whole degrees and is not used.

Samples are made daily by linear interpolation across gaps of at most
MAX_INTERP_DAYS; longer gaps stay missing. Same-day duplicate samples are averaged.
"""

from __future__ import annotations

import io
from collections.abc import Callable
from pathlib import Path

import pandas as pd

from oyster_data.paths import RAW_DIR
from oyster_data.raw import http_get, save_raw

URL = "https://www.dropbox.com/s/0oru9p1b7v5pc1r/data_dock.csv?dl=1"
STATION = "dmc_dock_chl"
RAW_FILE = "dmc_dock/data_dock.csv"
MAX_INTERP_DAYS = 7


def fetch(raw_dir: Path = RAW_DIR, get: Callable[[str], bytes] = http_get, overwrite: bool = False) -> Path:
    body = get(URL)
    if not body.startswith(b"Date,"):
        raise RuntimeError("unexpected DMC dock file format (no 'Date,' header)")
    return save_raw(body, raw_dir / RAW_FILE, overwrite=overwrite)


def parse(text: str) -> pd.Series:
    """Sample-day mean chlorophyll (ug/L) indexed by date."""
    df = pd.read_csv(io.StringIO(text))
    if "Date" not in df or "Ave. Chl" not in df:
        raise ValueError("DMC dock file lacks Date / Ave. Chl columns")
    df = df[["Date", "Ave. Chl"]].rename(columns={"Date": "date", "Ave. Chl": "chl_ug_l"})
    df["date"] = pd.to_datetime(df["date"], errors="coerce")
    df["chl_ug_l"] = pd.to_numeric(df["chl_ug_l"], errors="coerce")
    df = df.dropna()
    df = df[df["chl_ug_l"] >= 0]
    s = df.groupby("date")["chl_ug_l"].mean().astype(float)
    s.name = "chl_ug_l"
    return s


def to_daily(samples: pd.Series, max_gap_days: int = MAX_INTERP_DAYS) -> pd.Series:
    """Daily series: linear interpolation between samples no more than max_gap_days apart."""
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
    out.name = "chl_ug_l"
    return out


def load(raw_dir: Path = RAW_DIR) -> pd.Series:
    return to_daily(parse((raw_dir / RAW_FILE).read_text(encoding="utf-8")))
