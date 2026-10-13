"""NASA VIIRS (Suomi-NPP) 8-day 4 km chlorophyll-a near the Damariscotta mouth.

Secondary chlorophyll source, used after the DMC dock samples. ERDDAP griddap
dataset erdVH2018chla8day (NOAA CoastWatch West Coast node, "VIIRSN, Suomi-NPP,
Level-3 SMI, NASA, Global, 4km, Chlorophyll a, OCI Algorithm, R2018",
2012-01 to 2022-07), variable chla (mg m^-3 = ug/L).

Each 8-day composite is reduced to the mean of valid pixels in BOX (needs at
least MIN_VALID_PIXELS) and assigned to the 8 days it covers. ERDDAP stamps each
composite near its middle (the 1-8 Jan 2012 composite is stamped 2012-01-05), so
the covered days are stamp - 4 days through stamp + 3 days. Satellite chlorophyll near the coast is biased (CoastWatch values
were higher than dock samples), so it is only used through the regression on
overlapping days in oyster_data.gapfill, never directly.
"""

from __future__ import annotations

import io
from collections.abc import Callable
from pathlib import Path

import pandas as pd
import requests

from oyster_data.paths import RAW_DIR
from oyster_data.raw import http_get, save_raw

URL = "https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdVH2018chla8day.csv"
STATION = "viirs_8day_chl_box"
BOX = {"lat_max": 43.95, "lat_min": 43.70, "lon_min": -69.70, "lon_max": -69.45}
COMPOSITE_DAYS = 8
MIN_VALID_PIXELS = 5
RAW_SUBDIR = "viirs_chl"


def build_url(year: int) -> str:
    b = BOX
    return (f"{URL}?chla%5B({year}-01-01):1:({year}-12-31)%5D"
            f"%5B({b['lat_max']}):1:({b['lat_min']})%5D%5B({b['lon_min']}):1:({b['lon_max']})%5D")


def raw_path(year: int, raw_dir: Path = RAW_DIR) -> Path:
    return raw_dir / RAW_SUBDIR / f"erdVH2018chla8day_{year}.csv"


def fetch_year(year: int, raw_dir: Path = RAW_DIR, get: Callable[[str], bytes] = http_get,
               overwrite: bool = False) -> Path | None:
    """Download one year; returns None when ERDDAP reports no data for the year (HTTP 404)."""
    try:
        body = get(build_url(year))
    except requests.HTTPError as err:
        if err.response is not None and err.response.status_code == 404:
            return None
        raise
    return save_raw(body, raw_path(year, raw_dir), overwrite=overwrite)


def parse(text: str) -> pd.Series:
    """Daily chlorophyll (ug/L): box mean of each composite spread over its 8 days."""
    df = pd.read_csv(io.StringIO(text), skiprows=[1])
    if not {"time", "chla"} <= set(df.columns):
        raise ValueError("ERDDAP response lacks time/chla columns")
    df["chla"] = pd.to_numeric(df["chla"], errors="coerce")
    df = df.dropna(subset=["chla"])
    df = df[df["chla"] >= 0]
    stamp = pd.to_datetime(df["time"], utc=True).dt.tz_localize(None).dt.normalize()
    df["start"] = stamp - pd.Timedelta(days=COMPOSITE_DAYS // 2)
    g = df.groupby("start")["chla"]
    comp = g.mean()[g.count() >= MIN_VALID_PIXELS]
    days = {}
    for start, value in comp.items():
        for k in range(COMPOSITE_DAYS):
            days.setdefault(start + pd.Timedelta(days=k), float(value))
    s = pd.Series(days, dtype=float).sort_index()
    s.name = "chl_ug_l"
    return s


def load(raw_dir: Path = RAW_DIR) -> pd.Series:
    files = sorted((raw_dir / RAW_SUBDIR).glob("erdVH2018chla8day_*.csv"))
    if not files:
        raise FileNotFoundError(f"no VIIRS raw files in {raw_dir / RAW_SUBDIR}")
    s = pd.concat([parse(f.read_text(encoding="utf-8")) for f in files])
    return s[~s.index.duplicated()].sort_index()
