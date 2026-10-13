"""NERACOOS buoy E01 (Central Maine Shelf): secondary salinity and temperature source.

ERDDAP: https://data.neracoos.org/erddap/tabledap/E01_ocean_001m and E01_ocean_010m
Variables used: time (UTC), salinity (PSU), temperature (degC) and their
QARTOD aggregate flags (salinity_qc_agg, temperature_qc_agg). Only records with
flag 1 (pass) are kept. Half-hourly records are averaged to UTC days with at
least MIN_OBS_PER_DAY values.

Buoy position about 43.715 N, -69.355 W, roughly 30 km offshore of the
Damariscotta mouth, so values are mapped to sites by regression on overlapping
days (oyster_data.gapfill), never used directly. The 10 m sensor is a fallback
for periods when the 1 m sensor has no data (e.g. winter 2022).
"""

from __future__ import annotations

import io
from collections.abc import Callable
from pathlib import Path

import pandas as pd
import requests

from oyster_data.paths import RAW_DIR
from oyster_data.raw import http_get, save_raw

BASE = "https://data.neracoos.org/erddap/tabledap"
DATASETS = {1: "E01_ocean_001m", 10: "E01_ocean_010m"}
VARIABLES = {"salinity": ("salinity_psu", 0.0, 45.0), "temperature": ("temp_c", -2.0, 35.0)}
RAW_SUBDIR = "neracoos"
MIN_OBS_PER_DAY = 24  # of 48 half-hourly records


def station_name(variable: str, depth_m: int) -> str:
    return f"neracoos_e01_{depth_m}m_{variable}"


def build_url(year: int, depth_m: int) -> str:
    cols = "time,salinity,salinity_qc_agg,temperature,temperature_qc_agg"
    return (f"{BASE}/{DATASETS[depth_m]}.csv?{cols}"
            f"&time>={year}-01-01T00:00:00Z&time<{year + 1}-01-01T00:00:00Z")


def raw_path(year: int, depth_m: int, raw_dir: Path = RAW_DIR) -> Path:
    return raw_dir / RAW_SUBDIR / f"{DATASETS[depth_m]}_{year}.csv"


def fetch_year(year: int, depth_m: int = 1, raw_dir: Path = RAW_DIR, get: Callable[[str], bytes] = http_get,
               overwrite: bool = False) -> Path | None:
    """Download one year; returns None when ERDDAP reports no data for the year (HTTP 404)."""
    try:
        body = get(build_url(year, depth_m))
    except requests.HTTPError as err:
        if err.response is not None and err.response.status_code == 404:
            return None
        raise
    return save_raw(body, raw_path(year, depth_m, raw_dir), overwrite=overwrite)


def parse(text: str, variable: str) -> pd.Series:
    """Daily mean of a QC-passed variable ('salinity' or 'temperature') indexed by UTC date."""
    name, lo, hi = VARIABLES[variable]
    df = pd.read_csv(io.StringIO(text), skiprows=[1])   # second line holds units
    flag = f"{variable}_qc_agg"
    for col in ("time", variable, flag):
        if col not in df:
            raise ValueError(f"ERDDAP response lacks column {col!r}")
    df = df[df[flag] == 1].copy()
    df["date"] = pd.to_datetime(df["time"], utc=True).dt.tz_localize(None).dt.normalize()
    df[variable] = pd.to_numeric(df[variable], errors="coerce")
    df = df.dropna(subset=[variable])
    df = df[(df[variable] >= lo) & (df[variable] <= hi)]
    g = df.groupby("date")[variable]
    s = g.mean()[g.count() >= MIN_OBS_PER_DAY].astype(float)
    s.name = name
    return s


def load(variable: str, depth_m: int = 1, raw_dir: Path = RAW_DIR) -> pd.Series:
    files = sorted((raw_dir / RAW_SUBDIR).glob(f"{DATASETS[depth_m]}_*.csv"))
    if not files:
        raise FileNotFoundError(f"no {DATASETS[depth_m]} raw files in {raw_dir / RAW_SUBDIR}")
    s = pd.concat([parse(f.read_text(encoding="utf-8"), variable) for f in files])
    return s[~s.index.duplicated()].sort_index()
