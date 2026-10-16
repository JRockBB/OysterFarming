"""Great Bay, New Hampshire: Dermo (Perkinsus marinus) survey data and NERRS sonde forcing (PLAN.md P4.2d).

Disease: New Hampshire Fish and Game, Grout, D.E. and Smith, B.W. (2015) "Testing of Great Bay Oysters for Two
Protozoan Pathogens", UNH PREP Reports & Publications 356, Table 2 "Dermo Test Results - 1996 - 2014" (RFTM,
Haskin Shellfish Research Laboratory; counts of oysters per Mackin class 0.5, 1, 2, 3, 4, 5).
URL: ``REPORT_URL`` (Internet Archive copy; the repository page is behind a bot check).
Transcribed by hand from the PDF text (rows checked against the page). Kept: the five natural beds that have a
NERRS sonde within the same water body. Left out: rows whose class counts are inconsistent with n or with the
reported prevalence in the source (2004-11-18 and 2008-10-22 Nannie Island, 2008-10-27 Squamscott River), the
Piscataqua River bed (no sonde nearby), occasional beds and aquaculture lots.
Prevalence = infected / n; intensity_mean = sum(class x count) / n (zeros included), scale ``mackin_0_5``.

Bed locations: Adams Point, Nannie Island, Woodman Point from OpenStreetMap place names (Nominatim, 2026-10-14),
matching the report's Figure 1. Oyster River and Squamscott River beds: the report's Figure 1 places them in
those rivers; the NERRS sonde position in the same river is used (approximate).

Forcing: Great Bay NERR (Wells-style SWMP) sondes from the NCEI NERRS archive (``oyster_data.sources.nerrs_ncei``):
grbgbwq Great Bay (43 04 20 N, 70 52 10 W; 1995-2011), grborwq Oyster River (43.134 N, 70.911 W; 2002-2011),
grbsqwq Squamscott River (43 02 30 N, 70 55 20 W; 1997-2011), grblrwq Lamprey River (43 04 48 N, 70 56 04 W;
1998-2011); coordinates from the NERRS FGDC metadata (grbwq04-12.05m.xml). Sondes are removed in winter
(about 60 % of November-April days have no reading). Gaps are filled only from measured data: linear <= 7 days, then
regression on measured overlap days with, in order, the other Great Bay sondes, the Wells NERR inlet sonde
(``welinwq``, water temperature measured year-round 1995-2011, about 35 km south-east) and Great Bay air temperature
(``grbgbmet`` 2001-2005, ``grbglmet`` 2005-2010; the archive's grbglmet2011.csv is a mis-filed
water-quality file and is not used); no climatology. Salinity is filled from the other Great Bay sondes, then Wells inlet salinity (regression).
"""

from __future__ import annotations

from datetime import date
from pathlib import Path

import pandas as pd
import pyarrow as pa

from oyster_data import gapfill
from oyster_data.contract import SPECS, arrow_schema, validate_table
from oyster_data.paths import RAW_DIR
from oyster_data.raw import http_get, save_raw
from oyster_data.sources import nerrs_ncei

REPORT_URL = ("https://web.archive.org/web/20240226195037id_/"
              "https://scholars.unh.edu/cgi/viewcontent.cgi?article=1359&context=prep")
RAW_FILE = "nh_greatbay/nhfg_2015_testing_great_bay_oysters.pdf"
REGION = "great_bay_nh"
CLASSES = (0.5, 1, 2, 3, 4, 5)
FIRST_DATE, LAST_DATE = "1995-01-01", "2011-12-31"

SONDES = {"grbgbwq": (43.0722, -70.8694), "grborwq": (43.134, -70.911),
          "grbsqwq": (43.0417, -70.9222), "grblrwq": (43.08, -70.9344)}
SONDE_YEARS = {"grbgbwq": range(1995, 2012), "grborwq": range(2002, 2012),
               "grbsqwq": range(1997, 2012), "grblrwq": range(1998, 2012)}
# Winter temperature secondaries (station, variable, years): Wells inlet water temperature and Great Bay air temperature.
WINTER_SECONDARIES = [("welinwq", "Temp", range(1995, 2012)), ("grbgbmet", "ATemp", range(2001, 2006)),
                      ("grbglmet", "ATemp", range(2005, 2011))]   # the archive's grbglmet2011.csv is a mis-filed water-quality file

# site_id: (name, lat, lon, primary sonde, location note)
SITES = {
    "nh_adams_point": ("Adams Point (Great Bay)", 43.0908, -70.8643, "grbgbwq", "OpenStreetMap Adams Point, Durham"),
    "nh_nannie_island": ("Nannie Island (Great Bay)", 43.0691, -70.8628, "grbgbwq", "OpenStreetMap Nannie Island, Newington"),
    "nh_woodman_point": ("Woodman Point (Great Bay)", 43.0714, -70.8616, "grbgbwq", "OpenStreetMap Woodman Point, Newington"),
    "nh_oyster_river": ("Oyster River (Great Bay estuary)", 43.134, -70.911, "grborwq", "approximate: grborwq sonde position"),
    "nh_squamscott_river": ("Squamscott River (Great Bay estuary)", 43.0417, -70.9222, "grbsqwq", "approximate: grbsqwq sonde position"),
}

# (date, site_id, n, counts per Mackin class 0.5, 1, 2, 3, 4, 5) from Table 2
DERMO_ROWS = [
    ("1996-12-16", "nh_nannie_island", 25, (1, 0, 0, 0, 0, 0)),
    ("1997-08-17", "nh_adams_point", 25, (4, 0, 0, 0, 0, 0)),
    ("1997-08-14", "nh_oyster_river", 25, (1, 0, 0, 0, 0, 0)),
    ("1997-08-14", "nh_nannie_island", 25, (1, 0, 0, 0, 0, 0)),
    ("1997-09-08", "nh_squamscott_river", 25, (1, 0, 0, 0, 0, 0)),
    ("1997-11-17", "nh_adams_point", 25, (1, 0, 0, 0, 0, 0)),
    ("1997-11-17", "nh_nannie_island", 25, (0, 0, 0, 0, 0, 0)),
    ("1997-11-17", "nh_oyster_river", 25, (0, 0, 0, 0, 0, 0)),
    ("1998-12-09", "nh_adams_point", 25, (0, 0, 0, 0, 0, 0)),
    ("1998-12-09", "nh_nannie_island", 25, (0, 0, 0, 0, 0, 0)),
    ("1998-12-09", "nh_squamscott_river", 25, (0, 0, 0, 0, 0, 0)),
    ("1999-10-21", "nh_nannie_island", 20, (0, 0, 0, 0, 0, 0)),
    ("2000-11-04", "nh_adams_point", 20, (0, 0, 0, 0, 0, 0)),
    ("2000-11-04", "nh_nannie_island", 20, (0, 0, 0, 0, 0, 0)),
    ("2000-11-15", "nh_oyster_river", 20, (0, 0, 0, 0, 0, 0)),
    ("2001-10-10", "nh_nannie_island", 25, (0, 0, 0, 0, 0, 0)),
    ("2001-11-04", "nh_oyster_river", 20, (0, 0, 0, 0, 0, 0)),
    ("2001-11-04", "nh_adams_point", 20, (0, 0, 0, 0, 0, 0)),
    ("2002-10-14", "nh_adams_point", 20, (1, 2, 0, 0, 0, 0)),
    ("2002-10-14", "nh_oyster_river", 20, (0, 0, 0, 0, 0, 0)),
    ("2002-10-31", "nh_nannie_island", 24, (2, 0, 0, 0, 0, 0)),
    ("2003-10-28", "nh_nannie_island", 25, (2, 1, 0, 2, 0, 0)),
    ("2004-10-27", "nh_oyster_river", 25, (2, 0, 2, 0, 0, 0)),
    ("2004-11-19", "nh_adams_point", 20, (3, 4, 2, 4, 0, 0)),
    ("2005-11-06", "nh_oyster_river", 20, (3, 3, 5, 0, 2, 0)),
    ("2005-11-14", "nh_adams_point", 20, (6, 7, 3, 1, 1, 0)),
    ("2005-11-16", "nh_woodman_point", 20, (4, 4, 8, 2, 0, 0)),
    ("2005-11-17", "nh_squamscott_river", 20, (0, 1, 0, 0, 0, 0)),
    ("2006-11-01", "nh_oyster_river", 20, (3, 3, 4, 6, 0, 0)),
    ("2006-11-02", "nh_woodman_point", 20, (3, 8, 8, 1, 0, 0)),
    ("2006-11-07", "nh_squamscott_river", 39, (3, 1, 1, 0, 0, 0)),
    ("2006-11-22", "nh_adams_point", 20, (2, 8, 4, 5, 1, 0)),
    ("2006-12-07", "nh_nannie_island", 20, (2, 5, 4, 0, 1, 0)),
    ("2007-10-23", "nh_oyster_river", 20, (7, 1, 5, 4, 2, 1)),
    ("2007-10-24", "nh_woodman_point", 20, (3, 6, 1, 4, 3, 1)),
    ("2007-11-21", "nh_nannie_island", 20, (2, 0, 3, 0, 2, 0)),
    ("2007-12-07", "nh_adams_point", 20, (1, 1, 5, 2, 1, 1)),
    ("2008-10-08", "nh_adams_point", 20, (3, 3, 4, 4, 1, 1)),
    ("2008-10-09", "nh_woodman_point", 20, (1, 5, 0, 1, 0, 1)),
    ("2008-10-10", "nh_oyster_river", 20, (6, 7, 1, 2, 1, 0)),
    ("2009-11-04", "nh_oyster_river", 20, (3, 4, 5, 2, 3, 3)),
    ("2009-11-06", "nh_adams_point", 20, (3, 2, 6, 3, 1, 3)),
    ("2009-11-12", "nh_nannie_island", 20, (3, 9, 4, 0, 0, 0)),
    ("2009-11-13", "nh_woodman_point", 20, (0, 6, 4, 2, 1, 2)),
    ("2010-10-21", "nh_oyster_river", 20, (3, 6, 6, 2, 2, 0)),
    ("2010-10-19", "nh_adams_point", 20, (2, 7, 3, 1, 3, 2)),
    ("2010-10-20", "nh_nannie_island", 20, (1, 2, 8, 3, 1, 0)),
    ("2010-10-18", "nh_woodman_point", 20, (2, 4, 5, 3, 3, 2)),
    ("2010-11-16", "nh_squamscott_river", 20, (8, 3, 0, 0, 0, 0)),
    ("2011-10-21", "nh_adams_point", 20, (2, 4, 9, 1, 0, 1)),
    ("2011-10-26", "nh_oyster_river", 20, (3, 8, 2, 3, 2, 2)),
    ("2011-10-28", "nh_woodman_point", 20, (4, 5, 4, 6, 1, 0)),
    ("2011-11-04", "nh_nannie_island", 20, (6, 7, 4, 0, 1, 0)),
    ("2011-11-07", "nh_squamscott_river", 20, (9, 1, 3, 2, 1, 0)),
    ("2012-10-19", "nh_nannie_island", 10, (0, 1, 3, 3, 1, 0)),
    ("2012-10-25", "nh_woodman_point", 10, (0, 1, 2, 4, 1, 2)),
    ("2012-11-02", "nh_oyster_river", 10, (1, 3, 1, 2, 1, 1)),
    ("2012-11-19", "nh_adams_point", 10, (4, 1, 1, 0, 2, 0)),
    ("2012-12-04", "nh_squamscott_river", 10, (3, 0, 1, 1, 0, 0)),
    ("2013-10-29", "nh_woodman_point", 10, (3, 1, 3, 2, 0, 0)),
    ("2013-10-30", "nh_adams_point", 10, (3, 4, 1, 0, 0, 0)),
    ("2013-10-31", "nh_oyster_river", 10, (2, 5, 1, 0, 0, 0)),
    ("2013-11-04", "nh_nannie_island", 10, (3, 1, 4, 1, 0, 1)),
    ("2014-10-16", "nh_adams_point", 10, (0, 2, 2, 3, 1, 1)),
    ("2014-10-17", "nh_nannie_island", 9, (0, 0, 0, 3, 4, 1)),
    ("2014-10-28", "nh_oyster_river", 10, (0, 1, 5, 2, 1, 1)),
    ("2014-10-30", "nh_woodman_point", 10, (1, 3, 3, 1, 1, 1)),
]


def fetch_report(raw_dir: Path = RAW_DIR, get=http_get, overwrite: bool = False) -> Path:
    body = get(REPORT_URL)
    if not body.startswith(b"%PDF"):
        raise RuntimeError("NH F&G report download is not a PDF")
    return save_raw(body, raw_dir / RAW_FILE, overwrite=overwrite)


def fetch_sondes(raw_dir: Path = RAW_DIR, overwrite: bool = False) -> list[Path]:
    out = []
    for st, years in list(SONDE_YEARS.items()) + [(st, yrs) for st, _, yrs in WINTER_SECONDARIES]:
        for y in years:
            try:
                out.append(nerrs_ncei.fetch_year(st, y, raw_dir=raw_dir, overwrite=overwrite))
            except FileExistsError:
                pass
    return out


def check_rows() -> None:
    for d, s, n, counts in DERMO_ROWS:
        if s not in SITES:
            raise ValueError(f"unknown site {s}")
        if sum(counts) > n:
            raise ValueError(f"{d} {s}: {sum(counts)} infected > n = {n}")


def site_rows() -> list[dict]:
    return [{"site_id": k, "name": v[0], "lat": v[1], "lon": v[2], "region": REGION, "growing_area": None,
             "notes": f"Location: {v[4]}. Forcing: NERRS sonde {v[3]}."} for k, v in SITES.items()]


def disease_obs() -> pa.Table:
    check_rows()
    rows = []
    for d, s, n, counts in DERMO_ROWS:
        inf = sum(counts)
        rows.append({"site_id": s, "date": date.fromisoformat(d), "pathogen": "P_marinus", "prevalence": inf / n,
                     "intensity_mean": sum(c * k for c, k in zip(CLASSES, counts)) / n,
                     "intensity_scale": "mackin_0_5", "n": n,
                     "source": "nhfg:grout_smith_2015:Table2"})
    t = pa.Table.from_pylist(rows, schema=arrow_schema("disease_obs"))
    validate_table(t, SPECS["disease_obs"])
    return t


def build_forcing(raw_dir: Path = RAW_DIR) -> tuple[pd.DataFrame, list[gapfill.SiteReport]]:
    days = pd.date_range(FIRST_DATE, LAST_DATE, freq="D")
    load = {st: {v: nerrs_ncei.load(st, var, raw_dir).reindex(days)
                 for v, var in (("temp_c", "Temp"), ("salinity_psu", "Sal"))} for st in SONDES}
    parts, reports = [], []
    for sid, (_, _, _, primary, _) in SITES.items():
        daily = pd.DataFrame({"date": days, "site_id": sid, "temp_c": load[primary]["temp_c"].to_numpy(),
                              "salinity_psu": load[primary]["salinity_psu"].to_numpy(),
                              "chl_ug_l": float("nan"), "pom_mg_l": float("nan"),
                              "source": f"nerrs:{primary}", "qc_flag": 0})
        secondaries = {v: [(nerrs_ncei.station_name(st, var), load[st][v]) for st in SONDES if st != primary]
                       for v, var in (("temp_c", "Temp"), ("salinity_psu", "Sal"))}
        secondaries["temp_c"] += [(nerrs_ncei.station_name(st, var), nerrs_ncei.load(st, var, raw_dir).reindex(days))
                                  for st, var, _ in WINTER_SECONDARIES]
        # winter salinity: the contract needs salinity on every row, and the Great Bay sondes are all out in winter;
        # Wells inlet salinity is the only measured year-round series (fit statistics in forcing_report.md)
        secondaries["salinity_psu"].append((nerrs_ncei.station_name("welinwq", "Sal"),
                                            nerrs_ncei.load("welinwq", "Sal", raw_dir).reindex(days)))
        df, rep = gapfill.fill_site(daily, sid, secondaries=secondaries, climatology=False)
        parts.append(df)
        reports.append(rep)
    return pd.concat(parts, ignore_index=True), reports
