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
STATIONS = {"LE1.3": 1309, "YRK005.40": 1991, "WT7.1": 1730, "LE4.2": 1325}
UPPER_LAYERS = {"S", "AP"}


def build_url(station_id: int, start: date, end: date) -> str:
    subs = ",".join(str(v[1]) for v in SUBSTANCES.values())
    return f"{BASE}/{start:%m-%d-%Y}/{end:%m-%d-%Y}/0,1/2,4,6/2,3,7,23,24/Station/{station_id}/{subs}"


def raw_path(station: str, start: date, end: date, raw_dir: Path = RAW_DIR) -> Path:
    return raw_dir / RAW_SUBDIR / f"{station}_{start:%Y%m%d}_{end:%Y%m%d}.json"


def fetch(station: str, start: date, end: date, raw_dir: Path = RAW_DIR,
          get: Callable[[str], bytes] = http_get, overwrite: bool = False) -> Path:
    body = get(build_url(STATIONS[station], start, end))
    data = json.loads(body)
    if not isinstance(data, list):
        raise RuntimeError(f"DataHub {station}: unexpected response")
    return save_raw(body, raw_path(station, start, end, raw_dir), overwrite=overwrite)


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
