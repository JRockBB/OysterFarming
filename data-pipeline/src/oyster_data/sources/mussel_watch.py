"""NOAA Mussel Watch per-oyster Dermo scores, Chesapeake Bay sites (PLAN.md P4.2j).

Source: NOAA National Status and Trends Mussel Watch, East Coast oyster histopathology
(``MW_EastCoast_Oyster_Histopath.tsv``, NOAA NCCOS). The Chesapeake and Potomac rows are kept in
``data/literature/dermo/mussel_watch_chesapeake_oysters.csv`` (one row per oyster, 1988-2009, 9 sites
within 15 km of a Chesapeake Bay Program station). Each oyster carries a Dermo intensity on the Mackin
0-5 scale in thirds (``dermo_value``); usually 12 oysters per site and fiscal year.

Dates: the data give only the fiscal year. NS&T bivalve collections are made in the cold season; the
nominal date 1 January of the fiscal year is used (uncertainty about +/- 2 months), recorded in the
source string of every disease row.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from datetime import date
from pathlib import Path

import pyarrow as pa

from oyster_data.contract import SPECS, arrow_schema, validate_table

LIT_FILE = Path(__file__).resolve().parents[4] / "data" / "literature" / "dermo" / "mussel_watch_chesapeake_oysters.csv"
NOMINAL_MONTH, NOMINAL_DAY = 1, 1


def read_oysters(path: Path = LIT_FILE) -> list[dict]:
    with open(path, encoding="utf-8") as f:
        rows = [r for r in csv.DictReader(l for l in f if not l.startswith("#"))]
    for r in rows:
        r["dermo_value"] = float(r["dermo_value"])
        r["fiscal_year"] = int(r["fiscal_year"])
        r["lat"] = float(r["lat"])
        r["lon"] = float(r["lon"])
    return rows


def site_id(code: str) -> str:
    return "mw_" + code.lower()


def site_rows(path: Path = LIT_FILE) -> list[dict]:
    seen = {}
    for r in read_oysters(path):
        seen.setdefault(r["nst_site"], r)
    return [{"site_id": site_id(code), "name": f"Mussel Watch {r['specific_location']} ({code})", "lat": r["lat"],
             "lon": r["lon"], "region": "chesapeake", "growing_area": None,
             "notes": "NOAA Mussel Watch oyster histopathology site; per-oyster Dermo scores in data/literature/dermo"}
            for code, r in sorted(seen.items())]


def disease_obs(path: Path = LIT_FILE) -> pa.Table:
    """One row per site and fiscal year: prevalence (score > 0), mean score (mackin_0_5), n."""
    groups = defaultdict(list)
    for r in read_oysters(path):
        groups[(r["nst_site"], r["fiscal_year"])].append(r["dermo_value"])
    rows = []
    for (code, fy), vals in sorted(groups.items()):
        rows.append({"site_id": site_id(code), "date": date(fy, NOMINAL_MONTH, NOMINAL_DAY), "pathogen": "P_marinus",
                     "prevalence": sum(v > 0 for v in vals) / len(vals), "intensity_mean": sum(vals) / len(vals),
                     "intensity_scale": "mackin_0_5", "n": len(vals),
                     "source": "noaa_mussel_watch:oyster_histopath:fiscal_year_nominal_date_01-01"})
    t = pa.Table.from_pylist(rows, schema=arrow_schema("disease_obs"))
    validate_table(t, SPECS["disease_obs"])
    return t
