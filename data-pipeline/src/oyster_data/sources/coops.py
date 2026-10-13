"""NOAA CO-OPS water temperature (daily means of hourly records).

API: https://api.tidesandcurrents.noaa.gov/api/prod/datagetter
     ?product=water_temperature&station=<id>&begin_date=YYYYMMDD&end_date=YYYYMMDD
     &interval=h&units=metric&time_zone=gmt&format=csv&application=OysterFarming
One request covers at most one year. CSV columns: Date Time, Water Temperature,
X, N, R (flags). Records with X = 1 (exceeded the gross range limits) are dropped.
Days need at least MIN_OBS_PER_DAY hourly values.
"""

from __future__ import annotations

import io
from collections.abc import Callable
from pathlib import Path

import pandas as pd

from oyster_data.paths import RAW_DIR
from oyster_data.raw import http_get, save_raw

BASE = "https://api.tidesandcurrents.noaa.gov/api/prod/datagetter"
RAW_SUBDIR = "coops"
MIN_OBS_PER_DAY = 12

STATIONS = {
    "8577330": "Solomons Island, MD (Chesapeake Biological Laboratory pier)",
    "8637689": "Yorktown USCG Training Center / Gloucester Point, VA",
}


def station_name(station: str) -> str:
    return f"coops_{station}_water_temperature"


def build_url(station: str, year: int) -> str:
    return (f"{BASE}?product=water_temperature&station={station}&begin_date={year}0101&end_date={year}1231"
            "&interval=h&units=metric&time_zone=gmt&format=csv&application=OysterFarming")


def raw_path(station: str, year: int, raw_dir: Path = RAW_DIR) -> Path:
    return raw_dir / RAW_SUBDIR / f"{station}_water_temperature_{year}.csv"


def fetch_year(station: str, year: int, raw_dir: Path = RAW_DIR, get: Callable[[str], bytes] = http_get,
               overwrite: bool = False) -> Path:
    body = get(build_url(station, year))
    if not body.lstrip().startswith(b"Date Time"):
        raise RuntimeError(f"CO-OPS {station} {year}: unexpected response {body[:120]!r}")
    return save_raw(body, raw_path(station, year, raw_dir), overwrite=overwrite)


def parse(text: str) -> pd.Series:
    """Daily mean water temperature (degC) indexed by UTC date."""
    df = pd.read_csv(io.StringIO(text), skipinitialspace=True)
    df.columns = [c.strip() for c in df.columns]
    if "Water Temperature" not in df or "Date Time" not in df:
        raise ValueError("CO-OPS response lacks Date Time / Water Temperature")
    df["temp_c"] = pd.to_numeric(df["Water Temperature"], errors="coerce")
    if "X" in df:
        df = df[pd.to_numeric(df["X"], errors="coerce").fillna(0) != 1]
    df = df.dropna(subset=["temp_c"])
    df = df[(df["temp_c"] > -2) & (df["temp_c"] < 35)]
    df["date"] = pd.to_datetime(df["Date Time"]).dt.normalize()
    g = df.groupby("date")["temp_c"]
    s = g.mean()[g.count() >= MIN_OBS_PER_DAY].astype(float)
    s.name = "temp_c"
    return s


def load(station: str, raw_dir: Path = RAW_DIR) -> pd.Series:
    files = sorted((raw_dir / RAW_SUBDIR).glob(f"{station}_water_temperature_*.csv"))
    if not files:
        raise FileNotFoundError(f"no CO-OPS files for {station}")
    s = pd.concat([parse(f.read_text(encoding="utf-8")) for f in files])
    return s[~s.index.duplicated()].sort_index()
