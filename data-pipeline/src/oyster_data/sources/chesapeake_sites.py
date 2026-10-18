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

Maryland DNR disease bars (P4.2): a core set chosen for nearby long records
(2026-10-13 audit): temperature from the nearest CO-OPS gauge with other gauges
as regression fills, salinity from the nearest CBP fixed station. Chlorophyll is
carried along but not needed by the Dermo model.

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
    coops_secondaries: tuple[str, ...] = ()
    notes: tuple[str, ...] = field(default_factory=tuple)


SITES = [
    ChesSite("ches_patuxent_ks", "2005-10-01", "2007-10-31", "8577330", "LE1.3", "cbmrrwq",
             notes=("CO-OPS Solomons 3.4 km", "CBP LE1.3 0.9 km, S+AP layers",
                    "NERRS Maryland Railroad (Jug Bay, upper Patuxent) temperature for CO-OPS gaps")),
    ChesSite("ches_york_ks", "2005-10-01", "2006-10-31", "8637689", "YRK005.40", "cbvgiwq", "LE4.2", "cbvginut",
             notes=("CO-OPS Gloucester Point 0.1 km", "CBP YRK005.40 0.3 km (seasonal)",
                    "CBP LE4.2 monthly incl. winter", "NERRS Goodwin Islands ~10 km (sonde and monthly chl)")),
]

BAR_START, BAR_END = "1996-10-01", "2024-11-30"
BAR_SITES = [
    ChesSite("mddnr_hog_island", BAR_START, BAR_END, "8577330", "CB5.1W", coops_secondaries=("8573364",),
             notes=("CO-OPS Solomons 5.9 km", "CBP CB5.1W 1.3 km")),
    ChesSite("mddnr_broome_island", BAR_START, BAR_END, "8577330", "LE1.1", coops_secondaries=("8573364",),
             notes=("CO-OPS Solomons 14 km", "CBP LE1.1 3.8 km")),
    ChesSite("mddnr_swan_point", BAR_START, BAR_END, "8573364", "CB3.2", coops_secondaries=("8577330",),
             notes=("CO-OPS Tolchester 9.5 km", "CBP CB3.2 2.6 km")),
    ChesSite("mddnr_hackett_point", BAR_START, BAR_END, "8573364", "CB3.3W", coops_secondaries=("8577330",),
             notes=("CO-OPS Tolchester", "CBP CB3.3W 3.4 km")),
    ChesSite("mddnr_sandy_hill", BAR_START, BAR_END, "8571892", "ET5.2", coops_secondaries=("8577330",),
             notes=("CO-OPS Cambridge 5.4 km (from 1998)", "CBP ET5.2 5.3 km")),
    ChesSite("mddnr_oyster_shell_pt", BAR_START, BAR_END, "8571892", "ET5.2", coops_secondaries=("8577330",),
             notes=("CO-OPS Cambridge 5.6 km (from 1998)", "CBP ET5.2 5.1 km")),
    ChesSite("mddnr_sharkfin_shoal", BAR_START, BAR_END, "8571421", "EE3.1", coops_secondaries=("8571892", "8577330"),
             notes=("CO-OPS Bishops Head 4.1 km (patchy before 2012)", "CBP EE3.1 2.6 km")),
    ChesSite("mddnr_goose_creek", BAR_START, BAR_END, "8571421", "EE3.0", coops_secondaries=("8571892", "8577330"),
             notes=("CO-OPS Bishops Head 7.2 km (patchy before 2012)", "CBP EE3.0 1.3 km")),
]


# All other MD DNR disease bars (P4.2c data expansion, 2026-10-14): nearest CBP DataHub station with samples in
# both summer 2000 and summer 2015 (straight-line distance, no water-connectivity check), second nearest within
# 15 km as a regression fill; nearest CO-OPS temperature gauge, second nearest as fill. Bars whose nearest
# station is farther than 15 km are left out (mddnr_lancaster, 17.0 km).
BAR_SITES_EXTRA = [
    ChesSite("mddnr_holland_point", BAR_START, BAR_END, "8575512", "CB4.1W", cbp_secondary='CB4.2W', coops_secondaries=("8571892",),
             notes=("CO-OPS 8575512 27.7 km", "CBP CB4.1W 9.5 km", "CBP CB4.2W 10.2 km (secondary)")),
    ChesSite("mddnr_stone_rock", BAR_START, BAR_END, "8571892", "CB4.2E", cbp_secondary='CB4.2C', coops_secondaries=("8575512",),
             notes=("CO-OPS 8571892 29.4 km", "CBP CB4.2E 1.97 km", "CBP CB4.2C 3.48 km (secondary)")),
    ChesSite("mddnr_flag_pond", BAR_START, BAR_END, "8577330", "CB4.4", cbp_secondary='LE1.2', coops_secondaries=("8571892",),
             notes=("CO-OPS 8577330 13.2 km", "CBP CB4.4 8.18 km", "CBP LE1.2 9.07 km (secondary)")),
    ChesSite("mddnr_butler", BAR_START, BAR_END, "8635750", "CB5.2", cbp_secondary='LE2.3', coops_secondaries=("8577330",),
             notes=("CO-OPS 8635750 17.4 km", "CBP CB5.2 9.21 km", "CBP LE2.3 9.89 km (secondary)")),
    ChesSite("mddnr_buoy_rock", BAR_START, BAR_END, "8575512", "ET4.2", cbp_secondary='XGG8251', coops_secondaries=("8573364",),
             notes=("CO-OPS 8575512 23.2 km", "CBP ET4.2 0.34 km", "CBP XGG8251 3.99 km (secondary)")),
    ChesSite("mddnr_old_field", BAR_START, BAR_END, "8573364", "ET4.2", cbp_secondary='XGG8251', coops_secondaries=("8575512",),
             notes=("CO-OPS 8573364 16.4 km", "CBP ET4.2 10.69 km", "CBP XGG8251 14.07 km (secondary)")),
    ChesSite("mddnr_bugby", BAR_START, BAR_END, "8575512", "EE1.1", cbp_secondary='XGG8251', coops_secondaries=("8573364",),
             notes=("CO-OPS 8575512 25.0 km", "CBP EE1.1 2.54 km", "CBP XGG8251 10.15 km (secondary)")),
    ChesSite("mddnr_parsons_island", BAR_START, BAR_END, "8575512", "EE1.1", cbp_secondary='XGG8251', coops_secondaries=("8573364",),
             notes=("CO-OPS 8575512 20.4 km", "CBP EE1.1 3.15 km", "CBP XGG8251 7.49 km (secondary)")),
    ChesSite("mddnr_hollicutt_noose", BAR_START, BAR_END, "8575512", "CB4.1E", cbp_secondary='CB4.1C', coops_secondaries=("8571892",),
             notes=("CO-OPS 8575512 18.2 km", "CBP CB4.1E 4.33 km", "CBP CB4.1C 5.18 km (secondary)")),
    ChesSite("mddnr_bruffs_island", BAR_START, BAR_END, "8575512", "EE1.1", cbp_secondary='XGG8251', coops_secondaries=("8571892",),
             notes=("CO-OPS 8575512 28.5 km", "CBP EE1.1 5.62 km", "CBP XGG8251 13.42 km (secondary)")),
    ChesSite("mddnr_turtle_back", BAR_START, BAR_END, "8575512", "EE1.1", cbp_secondary='CB4.1E', coops_secondaries=("8571892",),
             notes=("CO-OPS 8575512 25.3 km", "CBP EE1.1 2.95 km", "CBP CB4.1E 12.18 km (secondary)")),
    ChesSite("mddnr_long_point", BAR_START, BAR_END, "8571892", "EE1.1", cbp_secondary='EE2.1', coops_secondaries=("8575512",),
             notes=("CO-OPS 8571892 24.1 km", "CBP EE1.1 13.86 km", "CBP EE2.1 14.96 km (secondary)")),
    ChesSite("mddnr_cook_point", BAR_START, BAR_END, "8571892", "EE2.1", cbp_secondary='CB4.2E', coops_secondaries=("8577330",),
             notes=("CO-OPS 8571892 21.8 km", "CBP EE2.1 2.27 km", "CBP CB4.2E 9.68 km (secondary)")),
    ChesSite("mddnr_royston", BAR_START, BAR_END, "8571892", "EE2.1", cbp_secondary='CB4.2E', coops_secondaries=("8575512",),
             notes=("CO-OPS 8571892 20.2 km", "CBP EE2.1 4.12 km", "CBP CB4.2E 14.64 km (secondary)")),
    ChesSite("mddnr_lighthouse", BAR_START, BAR_END, "8571892", "EE2.1", cbp_secondary='ET5.2', coops_secondaries=("8575512",),
             notes=("CO-OPS 8571892 14.6 km", "CBP EE2.1 6.51 km", "CBP ET5.2 14.21 km (secondary)")),
    ChesSite("mddnr_tilghman_wharf", BAR_START, BAR_END, "8571892", "EE2.1", cbp_secondary='CB4.2E', coops_secondaries=("8575512",),
             notes=("CO-OPS 8571892 27.4 km", "CBP EE2.1 8.12 km", "CBP CB4.2E 10.3 km (secondary)")),
    ChesSite("mddnr_deep_neck", BAR_START, BAR_END, "8571892", "EE2.1", cbp_secondary='CB4.1E', coops_secondaries=("8575512",),
             notes=("CO-OPS 8571892 24.2 km", "CBP EE2.1 9.44 km", "CBP CB4.1E 14.29 km (secondary)")),
    ChesSite("mddnr_double_mills", BAR_START, BAR_END, "8571892", "EE2.1", cbp_secondary=None, coops_secondaries=("8575512",),
             notes=("CO-OPS 8571892 18.8 km", "CBP EE2.1 13.6 km")),
    ChesSite("mddnr_cason", BAR_START, BAR_END, "8571892", "EE2.2", cbp_secondary='CB4.3E', coops_secondaries=("8577330",),
             notes=("CO-OPS 8571892 16.0 km", "CBP EE2.2 5.7 km", "CBP CB4.3E 13.47 km (secondary)")),
    ChesSite("mddnr_ragged_point_little_choptank", BAR_START, BAR_END, "8571892", "EE2.2", cbp_secondary='CB4.3E', coops_secondaries=("8577330",),
             notes=("CO-OPS 8571892 20.8 km", "CBP EE2.2 1.49 km", "CBP CB4.3E 8.41 km (secondary)")),
    ChesSite("mddnr_norman_addition", BAR_START, BAR_END, "8571421", "EE3.0", cbp_secondary=None, coops_secondaries=("8577330",),
             notes=("CO-OPS 8571421 9.5 km", "CBP EE3.0 11.46 km")),
    ChesSite("mddnr_wilson_shoals", BAR_START, BAR_END, "8571421", "ET6.2", cbp_secondary='EE3.0', coops_secondaries=("8571892",),
             notes=("CO-OPS 8571421 13.0 km", "CBP ET6.2 6.1 km", "CBP EE3.0 7.85 km (secondary)")),
    ChesSite("mddnr_georges", BAR_START, BAR_END, "8571421", "ET8.1", cbp_secondary='ET9.1', coops_secondaries=("8571892",),
             notes=("CO-OPS 8571421 19.1 km", "CBP ET8.1 4.03 km", "CBP ET9.1 9.07 km (secondary)")),
    ChesSite("mddnr_holland_straits", BAR_START, BAR_END, "8571421", "EE3.1", cbp_secondary='CB5.2', coops_secondaries=("8635750",),
             notes=("CO-OPS 8571421 12.4 km", "CBP EE3.1 12.95 km", "CBP CB5.2 13.65 km (secondary)")),
    ChesSite("mddnr_back_cove", BAR_START, BAR_END, "8571421", "EE3.2", cbp_secondary=None, coops_secondaries=("8635750",),
             notes=("CO-OPS 8571421 20.4 km", "CBP EE3.2 8.97 km")),
    ChesSite("mddnr_piney_island_east", BAR_START, BAR_END, "8571421", "EE3.2", cbp_secondary='ET9.1', coops_secondaries=("8635750",),
             notes=("CO-OPS 8571421 18.1 km", "CBP EE3.2 10.24 km", "CBP ET9.1 13.9 km (secondary)")),
    ChesSite("mddnr_old_womans_leg", BAR_START, BAR_END, "8571421", "EE3.2", cbp_secondary=None, coops_secondaries=("8635750",),
             notes=("CO-OPS 8571421 29.1 km", "CBP EE3.2 4.74 km")),
    ChesSite("mddnr_marumsco", BAR_START, BAR_END, "8571421", "EE3.4", cbp_secondary='EE3.3', coops_secondaries=("8635750",),
             notes=("CO-OPS 8571421 39.3 km", "CBP EE3.4 7.45 km", "CBP EE3.3 7.6 km (secondary)")),
    ChesSite("mddnr_chicken_cock", BAR_START, BAR_END, "8635750", "LE2.3", cbp_secondary='LE2.2', coops_secondaries=("8577330",),
             notes=("CO-OPS 8635750 14.3 km", "CBP LE2.3 13.73 km", "CBP LE2.2 14.6 km (secondary)")),
    ChesSite("mddnr_pagan", BAR_START, BAR_END, "8577330", "LE1.4", cbp_secondary='LE2.2', coops_secondaries=("8635750",),
             notes=("CO-OPS 8577330 13.9 km", "CBP LE1.4 13.51 km", "CBP LE2.2 14.06 km (secondary)")),
    ChesSite("mddnr_mills_west", BAR_START, BAR_END, "8635027", "RET2.4", cbp_secondary=None, coops_secondaries=("8577330",),
             notes=("CO-OPS 8635027 15.7 km", "CBP RET2.4 11.94 km")),
    ChesSite("mddnr_cornfield_harbor", BAR_START, BAR_END, "8635750", "LE2.3", cbp_secondary='CB5.2', coops_secondaries=("8577330",),
             notes=("CO-OPS 8635750 12.8 km", "CBP LE2.3 3.2 km", "CBP CB5.2 13.55 km (secondary)")),
    ChesSite("mddnr_ragged_point_potomac", BAR_START, BAR_END, "8635750", "LE2.2", cbp_secondary=None, coops_secondaries=("8577330",),
             notes=("CO-OPS 8635750 21.8 km", "CBP LE2.2 0.99 km")),
    ChesSite("mddnr_lower_cedar_point", BAR_START, BAR_END, "8635027", "RET2.4", cbp_secondary=None, coops_secondaries=("8577330",),
             notes=("CO-OPS 8635027 5.1 km", "CBP RET2.4 3.4 km")),
]


# NOAA Mussel Watch oyster sites in the Chesapeake (P4.2j; per-oyster Dermo scores 1988-2009): nearest long-record CBP
# station within 15 km, second nearest within 15 km as a regression fill, nearest CO-OPS gauges (straight-line distances).
MW_START, MW_END = "1987-01-01", "2010-12-31"
MW_SITES = [
    ChesSite("mw_cbbo", MW_START, MW_END, "8573364", "CB3.2", cbp_secondary="WT5.1", coops_secondaries=("8574680",),
             notes=("CO-OPS Tolchester 15 km", "CBP CB3.2 8.5 km", "CBP WT5.1 11.9 km (secondary)")),
    ChesSite("mw_cbcc", MW_START, MW_END, "8632200", "CB7.3E", cbp_secondary="CB7.2E", coops_secondaries=("8637689",),
             notes=("CO-OPS Kiptopeke 13 km", "CBP CB7.3E 7.1 km", "CBP CB7.2E 14.1 km (secondary)")),
    ChesSite("mw_cbcp", MW_START, MW_END, "8571892", "ET5.2", cbp_secondary="EE2.1", coops_secondaries=("8577330",),
             notes=("CO-OPS Cambridge 6 km", "CBP ET5.2 6.1 km", "CBP EE2.1 13.6 km (secondary)")),
    ChesSite("mw_cbdp", MW_START, MW_END, "8638610", "WE4.3", coops_secondaries=("8637689",),
             notes=("CO-OPS Sewells Point 18 km", "CBP WE4.3 11.1 km")),
    ChesSite("mw_cbhg", MW_START, MW_END, "8577330", "LE1.4", cbp_secondary="CB5.1W", coops_secondaries=("8571421",),
             notes=("CO-OPS Solomons 5 km", "CBP LE1.4 2.1 km", "CBP CB5.1W 2.4 km (secondary)")),
    ChesSite("mw_cbhp", MW_START, MW_END, "8575512", "CB3.3W", cbp_secondary="CB3.3C", coops_secondaries=("8573364",),
             notes=("CO-OPS Annapolis 6 km", "CBP CB3.3W 4.5 km", "CBP CB3.3C 5.6 km (secondary)")),
    ChesSite("mw_cbmp", MW_START, MW_END, "8575512", "CB3.3W", cbp_secondary="WT6.1", coops_secondaries=("8573364",),
             notes=("CO-OPS Annapolis 11 km", "CBP CB3.3W 7.8 km", "CBP WT6.1 8.4 km (secondary)")),
    ChesSite("mw_prrp", MW_START, MW_END, "8635750", "LE2.2", coops_secondaries=("8577330",),
             notes=("CO-OPS Lewisetta 21 km", "CBP LE2.2 0.4 km")),
    ChesSite("mw_prsp", MW_START, MW_END, "8635027", "RET2.4", coops_secondaries=("8577330",),
             notes=("CO-OPS Dahlgren 10 km", "CBP RET2.4 10.3 km")),
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
    for st in site.coops_secondaries:
        secondaries.setdefault("temp_c", []).append((coops.station_name(st), coops.load(st, raw_dir)))
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
    for site in SITES + BAR_SITES + BAR_SITES_EXTRA + MW_SITES:
        df, rep = build_site(site, raw_dir)
        parts.append(df)
        reports.append(rep)
    return pd.concat(parts, ignore_index=True), reports
