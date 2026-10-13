"""NERRS water-quality sondes from the NOAA NCEI archive (no CDMO registration needed).

Archive: https://www.ncei.noaa.gov/data/oceans/archive/arc0023/0052765/2.2/data/0-data/NERRS_CDMO_Archive_8.3.11/
Files: <station><year>.csv, e.g. cbvgiwq2006.csv (Chesapeake Bay NERR Virginia,
Goodwin Islands), 15-minute records with columns including Temp, F_Temp, Sal, F_Sal.
QC flags: <0> passed, <4> historical data predating automated QA/QC (all
2005-2006 data), <1> suspect, <-1> not collected, <-2> missing, <-3> rejected.
Only <0> and <4> are kept. Days need at least MIN_OBS_PER_DAY records (of 96).

Nutrient files (<station>nut<year>.csv, e.g. cbvginut2005.csv) hold monthly
grab samples; CHLA_N is lab chlorophyll-a (ug/L). Replicates on the same day are
averaged (parse_nutrient_chl).
"""

from __future__ import annotations

import io
from collections.abc import Callable
from pathlib import Path

import pandas as pd

from oyster_data.paths import RAW_DIR
from oyster_data.raw import http_get, save_raw

BASE = ("https://www.ncei.noaa.gov/data/oceans/archive/arc0023/0052765/2.2/data/0-data/"
        "NERRS_CDMO_Archive_8.3.11")
RAW_SUBDIR = "nerrs_ncei"
GOOD_FLAGS = {"0", "4"}
MIN_OBS_PER_DAY = 48
VARIABLES = {"Temp": ("temp_c", -2.0, 35.0), "Sal": ("salinity_psu", 0.0, 45.0)}


def station_name(station: str, variable: str) -> str:
    return f"nerrs_{station}_{VARIABLES[variable][0]}"


def build_url(station: str, year: int) -> str:
    return f"{BASE}/{station}{year}.csv"


def raw_path(station: str, year: int, raw_dir: Path = RAW_DIR) -> Path:
    return raw_dir / RAW_SUBDIR / f"{station}{year}.csv"


def fetch_year(station: str, year: int, raw_dir: Path = RAW_DIR, get: Callable[[str], bytes] = http_get,
               overwrite: bool = False) -> Path:
    return save_raw(get(build_url(station, year)), raw_path(station, year, raw_dir), overwrite=overwrite)


def _flag(series: pd.Series) -> pd.Series:
    return series.astype(str).str.extract(r"<(-?\d+)>", expand=False)


def parse(text: str, variable: str) -> pd.Series:
    """Daily mean of a QC-accepted variable ('Temp' or 'Sal') indexed by date (local standard time)."""
    name, lo, hi = VARIABLES[variable]
    df = pd.read_csv(io.StringIO(text), low_memory=False)
    df.columns = [c.strip() for c in df.columns]
    flag_col = f"F_{variable}"
    for col in ("DateTimeStamp", variable, flag_col):
        if col not in df:
            raise ValueError(f"NERRS file lacks column {col!r}")
    df = df[_flag(df[flag_col]).isin(GOOD_FLAGS)].copy()
    df[variable] = pd.to_numeric(df[variable], errors="coerce")
    df = df.dropna(subset=[variable])
    df = df[(df[variable] >= lo) & (df[variable] <= hi)]
    df["date"] = pd.to_datetime(df["DateTimeStamp"], format="%m/%d/%Y %H:%M").dt.normalize()
    g = df.groupby("date")[variable]
    s = g.mean()[g.count() >= MIN_OBS_PER_DAY].astype(float)
    s.name = name
    return s


def parse_nutrient_chl(text: str) -> pd.Series:
    """Sample-day mean lab chlorophyll-a (ug/L) from a NERRS nutrient file."""
    df = pd.read_csv(io.StringIO(text), low_memory=False)
    df.columns = [c.strip() for c in df.columns]
    for col in ("DateTimeStamp", "CHLA_N", "F_CHLA_N"):
        if col not in df:
            raise ValueError(f"NERRS nutrient file lacks column {col!r}")
    df = df[_flag(df["F_CHLA_N"]).isin(GOOD_FLAGS)].copy()
    df["CHLA_N"] = pd.to_numeric(df["CHLA_N"], errors="coerce")
    df = df[(df["CHLA_N"] >= 0)]
    df["date"] = pd.to_datetime(df["DateTimeStamp"], format="%m/%d/%Y %H:%M").dt.normalize()
    s = df.groupby("date")["CHLA_N"].mean().astype(float)
    s.name = "chl_ug_l"
    return s


def load_nutrient_chl(station: str, raw_dir: Path = RAW_DIR) -> pd.Series:
    files = sorted((raw_dir / RAW_SUBDIR).glob(f"{station}*.csv"))
    if not files:
        raise FileNotFoundError(f"no NERRS nutrient files for {station}")
    s = pd.concat([parse_nutrient_chl(f.read_text(encoding="utf-8", errors="replace")) for f in files])
    return s[~s.index.duplicated()].sort_index()


def load(station: str, variable: str, raw_dir: Path = RAW_DIR) -> pd.Series:
    files = sorted((raw_dir / RAW_SUBDIR).glob(f"{station}*.csv"))
    if not files:
        raise FileNotFoundError(f"no NERRS files for {station}")
    s = pd.concat([parse(f.read_text(encoding="utf-8", errors="replace"), variable) for f in files])
    return s[~s.index.duplicated()].sort_index()
