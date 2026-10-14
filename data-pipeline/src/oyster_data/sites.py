"""Site registry (PLAN.md P1.1). Writes ``data/processed/sites.arrow``.

Only sites with a verified position and a forcing source are listed. Each
site's ``notes`` names its forcing source and where its coordinates came from.
"""

from __future__ import annotations

from pathlib import Path

import pyarrow as pa

from oyster_data.contract import arrow_schema, write_contract
from oyster_data.paths import PROCESSED_DIR

SITES: list[dict[str, object]] = [
    {
        "site_id": "dam_lower",
        "name": "Lower Damariscotta River (LOBO node 65)",
        "lat": 43.93087,
        "lon": -69.58012,
        "region": "maine",
        "growing_area": "WQ",  # Maine DMR NSSP classification layers 2021-2025: area WQ, Approved
        "notes": "forcing: lobo:node0065 (Maine EPSCoR LOBO, http://maine.loboviz.com). "
                 "Position: median of the buoy's own GPS fixes 2015-2023. "
                 "Gap-fill stations (regression on overlap days, see forcing_report.md): Maine DMR Boothbay Harbor SST; NERACOOS E01 1 m and 10 m temperature/salinity; DMC dock chlorophyll; VIIRS 8-day chlorophyll.",
    },
    {
        "site_id": "dam_upper",
        "name": "Upper Damariscotta Estuary, Mook Sea Farm/Perkins Point (LOBO node 66)",
        "lat": 43.99956,
        "lon": -69.54104,
        "region": "maine",
        "growing_area": "WQ",  # Maine DMR NSSP classification layers 2021-2025: area WQ, Approved
        "notes": "forcing: lobo:node0066 (Maine EPSCoR LOBO, http://maine.loboviz.com). "
                 "Position: median of the buoy's own GPS fixes 2015-2022. "
                 "Gap-fill stations (regression on overlap days, see forcing_report.md): Maine DMR Boothbay Harbor SST; NERACOOS E01 1 m and 10 m temperature/salinity; DMC dock chlorophyll; VIIRS 8-day chlorophyll.",
    },
    # Chesapeake growth-calibration sites (PLAN.md P2.6). Coordinates are the
    # Latitude/Longitude given for these cohorts in Zenodo 10.5281/zenodo.11205266.
    {
        "site_id": "ches_patuxent_ks",
        "name": "Patuxent River, Kingsley-Smith triploid cage site",
        "lat": 38.342046,
        "lon": -76.475113,
        "region": "chesapeake",
        "growing_area": None,
        "notes": "growth: zenodo11205266 (Kingsley-Smith); forcing: see PLAN.md P2.6.",
    },
    {
        "site_id": "ches_severn_ks",
        "name": "Severn River, Kingsley-Smith triploid cage site",
        "lat": 38.981848,
        "lon": -76.459049,
        "region": "chesapeake",
        "growing_area": None,
        "notes": "growth: zenodo11205266 (Kingsley-Smith); forcing: see PLAN.md P2.6.",
    },
    {
        "site_id": "ches_york_ks",
        "name": "York River, Kingsley-Smith triploid cage site",
        "lat": 37.2471916666667,
        "lon": -76.5002,
        "region": "chesapeake",
        "growing_area": None,
        "notes": "growth: zenodo11205266 (Kingsley-Smith); forcing: see PLAN.md P2.6.",
    },
]


def sites_table() -> pa.Table:
    from oyster_data.sources.mddnr_fall_survey import site_rows as disease_bar_rows
    from oyster_data.sources.nh_greatbay import site_rows as greatbay_rows

    schema = arrow_schema("sites")
    return pa.Table.from_pylist(SITES + disease_bar_rows() + greatbay_rows(), schema=schema)


def write_sites(dest: Path = PROCESSED_DIR / "sites.arrow") -> Path:
    write_contract(sites_table(), dest, "sites")
    return dest
