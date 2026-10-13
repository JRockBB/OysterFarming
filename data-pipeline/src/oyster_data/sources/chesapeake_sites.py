"""Daily forcing for the Chesapeake growth-calibration sites (PLAN.md P2.6).

Per site (choices made from the 2026-10-13 forcing audit; all inputs are measurements):
- temperature: daily means of the nearest NOAA CO-OPS gauge; gaps longer than
  7 days are filled by regression on a NERRS sonde in the same estuary
  (fitted on overlap days, statistics in forcing_report.md).
- salinity and chlorophyll-a: Chesapeake Bay Program grab samples at the
  nearest monitoring station (lab-extracted chlorophyll), linearly interpolated
  between sampling dates up to GRAB_MAX_GAP_DAYS apart. Sample days count as
  observed (qc_flag 0); interpolated days are filled (qc_flag 1). Where the
  primary station has a seasonal gap, a second CBP station's grab series is
  used through regression on overlap days. The limit is about two months so
  that one skipped monthly sample does not break the series.

Severn River (ches_severn_ks) is not built: no station operated in the Severn
in 2005-2007 and the nearest winter measurements are about 30 km away.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from oyster_data import gapfill
from oyster_data.forcing import interpolate_samples
from oyster_data.paths import RAW_DIR
from oyster_data.sources import cbp_datahub, coops, nerrs_ncei

GRAB_MAX_GAP_DAYS = 62


@dataclass(frozen=True)
class ChesSite:
    site_id: str
    start: str
    end: str
    coops_station: str
    cbp_station: str
    nerrs_temp_station: str | None = None
    cbp_secondary: str | None = None
    nerrs_nut_station: str | None = None
    notes: tuple[str, ...] = field(default_factory=tuple)


SITES = [
    ChesSite("ches_patuxent_ks", "2005-10-01", "2007-10-31", "8577330", "LE1.3", "cbmrrwq",
             notes=("CO-OPS Solomons 3.4 km", "CBP LE1.3 0.9 km, S+AP layers",
                    "NERRS Maryland Railroad (Jug Bay, upper Patuxent) temperature for CO-OPS gaps")),
    ChesSite("ches_york_ks", "2005-10-01", "2006-10-31", "8637689", "YRK005.40", "cbvgiwq", "LE4.2", "cbvginut",
             notes=("CO-OPS Gloucester Point 0.1 km", "CBP YRK005.40 0.3 km (seasonal)",
                    "CBP LE4.2 monthly incl. winter", "NERRS Goodwin Islands ~10 km (sonde and monthly chl)")),
]


def build_site(site: ChesSite, raw_dir=RAW_DIR) -> tuple[pd.DataFrame, gapfill.SiteReport]:
    days = pd.date_range(site.start, site.end, freq="D")
    temp = coops.load(site.coops_station, raw_dir).reindex(days)
    grabs = cbp_datahub.load(site.cbp_station, raw_dir).reindex(days)
    daily = pd.DataFrame({
        "site_id": site.site_id,
        "date": days.date,
        "temp_c": temp.to_numpy(),
        "salinity_psu": grabs["salinity_psu"].to_numpy(),
        "chl_ug_l": grabs["chl_ug_l"].to_numpy(),
        "pom_mg_l": None,
        "source": f"coops:{site.coops_station}+cbp:{site.cbp_station}",
        "qc_flag": 0,
    })
    secondaries: dict[str, list[gapfill.Secondary]] = {}
    if site.nerrs_temp_station:
        secondaries["temp_c"] = [(nerrs_ncei.station_name(site.nerrs_temp_station, "Temp"),
                                  nerrs_ncei.load(site.nerrs_temp_station, "Temp", raw_dir))]
    if site.cbp_secondary:
        sec = cbp_datahub.load(site.cbp_secondary, raw_dir)
        for var in ("salinity_psu", "chl_ug_l"):
            series = interpolate_samples(sec[var].dropna(), GRAB_MAX_GAP_DAYS)
            secondaries.setdefault(var, []).append((f"cbp_{site.cbp_secondary}_{var}", series))
    if site.nerrs_nut_station:
        nut = interpolate_samples(nerrs_ncei.load_nutrient_chl(site.nerrs_nut_station, raw_dir), GRAB_MAX_GAP_DAYS)
        secondaries.setdefault("chl_ug_l", []).append((f"nerrs_{site.nerrs_nut_station}_chl_ug_l", nut))
    max_linear = {"salinity_psu": GRAB_MAX_GAP_DAYS, "chl_ug_l": GRAB_MAX_GAP_DAYS}
    return gapfill.fill_site(daily, site.site_id, secondaries, max_linear=max_linear, climatology=False)


def build_all(raw_dir=RAW_DIR) -> tuple[pd.DataFrame, list[gapfill.SiteReport]]:
    parts, reports = [], []
    for site in SITES:
        df, rep = build_site(site, raw_dir)
        parts.append(df)
        reports.append(rep)
    return pd.concat(parts, ignore_index=True), reports
