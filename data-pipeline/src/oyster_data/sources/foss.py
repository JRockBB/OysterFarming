"""Maine eastern oyster landings from NOAA Fisheries One Stop Shop (FOSS), PLAN.md P1.6.

API: https://apps-st.fisheries.noaa.gov/ods/foss/landings/ (Oracle REST, JSON).
Query used: state_name = MAINE and species name like %OYSTER%. Only
"OYSTER, EASTERN" (Crassostrea virginica, TSN 79872) is kept; "OYSTER, EDIBLE"
is the European flat oyster.

``pounds`` are written as reported by FOSS. NOAA commercial landings for
molluscs are customarily reported as meat weight; check this before computing a
price per oyster. FOSS gives no oyster counts (``tot_count`` is null), so
``count`` is null. Records with no reported dollar value are dropped because
``value_usd`` is required by the contract.
"""

from __future__ import annotations

import json
from pathlib import Path
from urllib.parse import urlencode

import pyarrow as pa

from oyster_data.contract import arrow_schema, validate_table, SPECS
from oyster_data.paths import RAW_DIR
from oyster_data.raw import Getter, http_get, save_raw

BASE_URL = "https://apps-st.fisheries.noaa.gov/ods/foss/landings/"
SPECIES = "OYSTER, EASTERN"
RAW_FILE = "foss/maine_oyster_landings.json"


def build_url(limit: int = 500) -> str:
    q = json.dumps({"state_name": "MAINE", "ts_afs_name": {"$like": "%OYSTER%"}}, separators=(",", ":"))
    return BASE_URL + "?" + urlencode({"q": q, "limit": limit})


def fetch(raw_dir: Path = RAW_DIR, get: Getter = http_get, overwrite: bool = False) -> Path:
    body = get(build_url())
    payload = json.loads(body)
    if payload.get("hasMore"):
        raise RuntimeError("FOSS response is paginated (hasMore=true); raise the limit")
    return save_raw(body, raw_dir / RAW_FILE, overwrite=overwrite)


def to_economics(payload: dict) -> pa.Table:
    """Contract ``economics`` rows (region ``maine``) from a FOSS JSON payload.

    Some historical years (before 2000, source ACCSPHIST) come as several
    partial records; they are summed into one annual total, and ``source``
    records how many records were combined.
    """
    by_year: dict[int, list[dict]] = {}
    for item in payload["items"]:
        if item.get("ts_afs_name") != SPECIES or item.get("dollars") is None:
            continue
        by_year.setdefault(int(item["year"]), []).append(item)
    rows = []
    for year in sorted(by_year):
        items = by_year[year]
        pounds = [i.get("pounds") for i in items]
        origins = sorted({str(i.get("source", "")) for i in items})
        rows.append({
            "year": year,
            "region": "maine",
            "value_usd": float(sum(i["dollars"] for i in items)),
            "pounds": None if any(v is None for v in pounds) else float(sum(pounds)),
            "count": None,
            "source": f"noaa_foss:{'+'.join(origins)}:{SPECIES}:records={len(items)}",
        })
    table = pa.Table.from_pylist(rows, schema=arrow_schema("economics"))
    validate_table(table, SPECS["economics"])
    return table


def clean(raw_dir: Path = RAW_DIR) -> pa.Table:
    return to_economics(json.loads((raw_dir / RAW_FILE).read_text(encoding="utf-8")))
