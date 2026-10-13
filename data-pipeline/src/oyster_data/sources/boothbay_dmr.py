"""Maine DMR Boothbay Harbor daily sea surface temperature (secondary temperature source).

Service: Maine DMR Open Data, ArcGIS FeatureServer
  https://services1.arcgis.com/RbMX0mRVOFNTdLzd/arcgis/rest/services/
  MaineDMR_Boothbay_Harbor_Sea_Surface_Temperatures/FeatureServer/0
Fields: COLLECTION_DATE (epoch ms at local midnight, US Eastern),
SEA_SURFACE_TEMP_MAX_C / _MIN_C / _AVG_C. Record 1905-present.
Station: West Boothbay Harbor pier, about 43.844 N, -69.641 W, sensor at
-5.5 ft MLW (Maine DMR Boothbay Harbor environmental data page). About 11 km
from the Damariscotta sites, so values are mapped to each site by regression on
overlapping days (oyster_data.gapfill), never used directly.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from datetime import date
from pathlib import Path
from urllib.parse import urlencode

import pandas as pd

from oyster_data.paths import RAW_DIR
from oyster_data.raw import http_get, save_raw

LAYER_URL = ("https://services1.arcgis.com/RbMX0mRVOFNTdLzd/arcgis/rest/services/"
             "MaineDMR_Boothbay_Harbor_Sea_Surface_Temperatures/FeatureServer/0")
STATION = "boothbay_dmr_sst"
PAGE_SIZE = 2000
RAW_SUBDIR = "boothbay_dmr"


def query_url(start: date, end: date, offset: int = 0) -> str:
    where = f"COLLECTION_DATE >= DATE '{start:%Y-%m-%d}' AND COLLECTION_DATE < DATE '{end:%Y-%m-%d}'"
    params = {"where": where, "outFields": "COLLECTION_DATE,SEA_SURFACE_TEMP_AVG_C",
              "orderByFields": "COLLECTION_DATE", "resultOffset": offset,
              "resultRecordCount": PAGE_SIZE, "f": "json"}
    return LAYER_URL + "/query?" + urlencode(params)


def raw_path(start: date, end: date, raw_dir: Path = RAW_DIR) -> Path:
    return raw_dir / RAW_SUBDIR / f"sst_{start:%Y%m%d}_{end:%Y%m%d}.json"


def fetch(start: date, end: date, raw_dir: Path = RAW_DIR, get: Callable[[str], bytes] = http_get,
          overwrite: bool = False) -> Path:
    """Download all features in [start, end) (paged) into one raw JSON file."""
    features: list[dict] = []
    offset = 0
    while True:
        page = json.loads(get(query_url(start, end, offset)))
        if "error" in page:
            raise RuntimeError(f"ArcGIS error: {page['error']}")
        features.extend(page.get("features", []))
        if not page.get("exceededTransferLimit"):
            break
        offset += PAGE_SIZE
    body = json.dumps({"layer": LAYER_URL, "start": str(start), "end": str(end), "features": features})
    return save_raw(body.encode("utf-8"), raw_path(start, end, raw_dir), overwrite=overwrite)


def parse(payload: dict) -> pd.Series:
    """Daily mean temperature (degC) indexed by local calendar date."""
    rows = [(f["attributes"]["COLLECTION_DATE"], f["attributes"]["SEA_SURFACE_TEMP_AVG_C"])
            for f in payload["features"]]
    df = pd.DataFrame(rows, columns=["ms", "temp_c"]).dropna()
    local = pd.to_datetime(df["ms"], unit="ms", utc=True).dt.tz_convert("America/New_York")
    df["date"] = local.dt.normalize().dt.tz_localize(None)
    s = df.groupby("date")["temp_c"].mean().astype(float)
    s = s[(s > -2) & (s < 35)]
    s.name = "temp_c"
    return s


def load(raw_dir: Path = RAW_DIR) -> pd.Series:
    files = sorted((raw_dir / RAW_SUBDIR).glob("sst_*.json"))
    if not files:
        raise FileNotFoundError(f"no Boothbay raw files in {raw_dir / RAW_SUBDIR}")
    parts = [parse(json.loads(f.read_text(encoding="utf-8"))) for f in files]
    s = pd.concat(parts)
    return s[~s.index.duplicated()].sort_index()
