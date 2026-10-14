"""Dated Maine DMR closure events for the Maine sites' growing areas, transcribed from archived notices
(PLAN.md P1.5, biotoxin and emergency part).

DMR publishes only currently active closures, so past events were reconstructed from Internet Archive
captures of DMR closure pages, legal notices and newsletters (search done 2026-10-14). Only events whose
notices name the growing area are listed. Searched and not listed:
- Biotoxin (PSP/ASP) inventories 2009-2021 (pspclosures.htm, psp.html, hotline.html; about 140 captures):
  American oysters are never listed for the Damariscotta. The 2017 PSP oyster closure named the Johns and
  Pemaquid Rivers and Casco Bay, not the Damariscotta; the 2017-18 ASP closure was Casco Bay only; the 2019
  oyster PSP event was Harpswell/Brunswick.
- 2025 PSP oyster closure ("ME/NH border through Pemaquid", late May to 2025-07-11; DMR Summer 2025
  newsletter): whether it covered the Damariscotta is unknown (no archived 100-A notice for those dates).
"""

from __future__ import annotations

from datetime import date

import pyarrow as pa

from oyster_data.contract import SPECS, arrow_schema, validate_table

EVENTS = [
    {
        "growing_area": "WQ",
        "area_id": "WQ:upper_river_oil_2021",
        "start_date": date(2021, 6, 25),
        "end_date": date(2021, 6, 28),
        "closure_type": "pollution",
        "source": ("maine_dmr:WQ notice 2021-06-28 (reclassified upper Damariscotta River P4, Little Point Rd "
                   "(Edgecomb)-Hall Pt line to Perkins Pt Rd (Edgecomb)-Tudor Rd (Bristol) line, Prohibited to "
                   "Approved; oil from a sunken vessel), "
                   "https://web.archive.org/web/20220123024224/https://www.maine.gov/dmr/shellfish-sanitation-management/closures/documents/WQ.pdf; "
                   "start date from Lincoln County News, https://lcnme.com/homepage/boat-sinks-near-oyster-farm-in-damariscotta-river"),
    },
]

NOTES = [
    "WQ:upper_river_oil_2021: dam_lower (43.931 N) is well south of the closed reach. dam_upper (43.9996 N) lies "
    "at the southern boundary line (Perkins Point Rd, about 43.9999 N in OpenStreetMap), so whether it was "
    "inside is not resolved.",
    "Biotoxin: no American-oyster biotoxin closure naming the Damariscotta found 2009-2021; 2025 PSP closure "
    "coverage of the Damariscotta unknown. Settling 2017 and 2025 needs DMR records (FOAA request to "
    "DMRPublicHealthDiv@maine.gov).",
]


def to_closures(growing_areas: list[str]) -> pa.Table:
    rows = [e for e in EVENTS if e["growing_area"] in growing_areas]
    table = pa.Table.from_pylist(rows, schema=arrow_schema("closures"))
    validate_table(table, SPECS["closures"])
    return table
