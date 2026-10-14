"""Chesapeake Bay Program DataHub water-quality grab samples (temperature, salinity, chlorophyll-a).

API (v5): https://datahub.chesapeakebay.net/api/WaterQuality/{start mm-dd-yyyy}/{end mm-dd-yyyy}/
          0,1/{program ids}/{project ids}/Station/{station ids}/{substance ids}
Programs 2,4,6 (NTWQM, SWM, TWQM); projects 2,3,7,23,24 (CMON, DFLO, MAIN, TRIB, ...);
substances 123 = WTEMP (degC), 83 = SALINITY (ppt), 21 = CHLA (ug/L, extracted).
Response: JSON list of records with SampleDate, layer, parameter, measureValue,
qualifier and problem fields.

Layer rule (a processing choice, recorded in PLAN.md P2.6): at stratified
mid-channel stations that report an above-pycnocline layer (AP), only surface
(S) and above-pycnocline (AP) samples are used, as they represent water at
shallow subtidal cage sites better than water below the pycnocline. At other
stations all layers are used. Samples with a qualifier (e.g. "<" below
detection) or a problem code are dropped. One value per sampling date: the mean
over the selected samples.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from datetime import date
from pathlib import Path

import pandas as pd

from oyster_data.paths import RAW_DIR
from oyster_data.raw import http_get, save_raw

BASE = "https://datahub.chesapeakebay.net/api/WaterQuality"
RAW_SUBDIR = "cbp_datahub"
SUBSTANCES = {"WTEMP": ("temp_c", 123), "SALINITY": ("salinity_psu", 83), "CHLA": ("chl_ug_l", 21)}
# DataHub station ids from https://datahub.chesapeakebay.net/api/GeographicAttribute/Station
STATIONS = {
    "LE1.3": 1309,
    "YRK005.40": 1991,
    "WT7.1": 1730,
    "LE4.2": 1325,
    "CB5.1W": 1171,
    "LE1.1": 1307,
    "CB3.2": 1153,
    "CB3.3W": 1156,
    "ET5.2": 1257,
    "EE3.1": 1235,
    "EE3.0": 1234,
    "CB4.1W": 1162,
    "CB4.2W": 1165,
    "CB4.2E": 1164,
    "CB4.2C": 1163,
    "CB4.4": 1169,
    "LE1.2": 1308,
    "CB5.2": 1172,
    "LE2.3": 1312,
    "ET4.2": 1253,
    "XGG8251": 1888,
    "EE1.1": 1231,
    "CB4.1E": 1161,
    "CB4.1C": 1160,
    "EE2.1": 1232,
    "EE2.2": 1233,
    "CB4.3E": 1167,
    "ET6.2": 1259,
    "ET8.1": 1261,
    "ET9.1": 1262,
    "EE3.2": 1236,
    "EE3.4": 1238,
    "EE3.3": 1237,
    "LE2.2": 1311,
    "LE1.4": 1310,
    "RET2.4": 1479,
}
UPPER_LAYERS = {"S", "AP"}


def build_url(station_id: int, start: date, end: date) -> str:
    subs = ",".join(str(v[1]) for v in SUBSTANCES.values())
    return f"{BASE}/{start:%m-%d-%Y}/{end:%m-%d-%Y}/0,1/2,4,6/2,3,7,23,24/Station/{station_id}/{subs}"


def raw_path(station: str, start: date, end: date, raw_dir: Path = RAW_DIR) -> Path:
    return raw_dir / RAW_SUBDIR / f"{station}_{start:%Y%m%d}_{end:%Y%m%d}.json"


MAX_YEARS_PER_QUERY = 10  # API limit: "The maximum query year range is 10 years" (HTTP 400 beyond it)


def date_chunks(start: date, end: date, max_years: int = MAX_YEARS_PER_QUERY) -> list[tuple[date, date]]:
    """Split [start, end] into pieces spanning at most `max_years` calendar years."""
    if end < start:
        raise ValueError("end before start")
    chunks, a = [], start
    while a <= end:
        b = min(end, date(a.year + max_years - 1, 12, 31))
        chunks.append((a, b))
        a = date(b.year + 1, 1, 1)
    return chunks


def fetch(station: str, start: date, end: date, raw_dir: Path = RAW_DIR,
          get: Callable[[str], bytes] = http_get, overwrite: bool = False) -> list[Path]:
    """Download [start, end] in API-sized pieces, one raw file per piece.

    Pieces already on disk are kept unless `overwrite`; raises FileExistsError
    if every piece was already present."""
    saved = []
    for a, b in date_chunks(start, end):
        path = raw_path(station, a, b, raw_dir)
        if path.exists() and not overwrite:
            continue
        body = get(build_url(STATIONS[station], a, b))
        if not isinstance(json.loads(body), list):
            raise RuntimeError(f"DataHub {station}: unexpected response")
        saved.append(save_raw(body, path, overwrite=overwrite))
    if not saved:
        raise FileExistsError(f"DataHub {station} {start}..{end}: all pieces already saved")
    return saved


def parse(records: list[dict]) -> pd.DataFrame:
    """Per sampling date means: columns temp_c, salinity_psu, chl_ug_l (NaN where not sampled)."""
    df = pd.DataFrame(records)
    if df.empty:
        return pd.DataFrame(columns=[v[0] for v in SUBSTANCES.values()])
    df = df[df["parameter"].isin(SUBSTANCES)]
    clean = df["qualifier"].fillna("").astype(str).str.strip().eq("") & df["problem"].fillna("").astype(str).str.strip().eq("")
    df = df[clean].copy()
    if (df["layer"] == "AP").any():
        df = df[df["layer"].isin(UPPER_LAYERS)]
    df["measureValue"] = pd.to_numeric(df["measureValue"], errors="coerce")
    df["date"] = pd.to_datetime(df["SampleDate"].str[:10])
    wide = df.pivot_table(index="date", columns="parameter", values="measureValue", aggfunc="mean")
    wide = wide.rename(columns={k: v[0] for k, v in SUBSTANCES.items()})
    for _, (col, _) in SUBSTANCES.items():
        if col not in wide:
            wide[col] = float("nan")
    return wide[[v[0] for v in SUBSTANCES.values()]].sort_index()


def load(station: str, raw_dir: Path = RAW_DIR) -> pd.DataFrame:
    files = sorted((raw_dir / RAW_SUBDIR).glob(f"{station}_*.json"))
    if not files:
        raise FileNotFoundError(f"no DataHub files for {station}")
    records = [r for f in files for r in json.loads(f.read_text(encoding="utf-8"))]
    return parse(records)
