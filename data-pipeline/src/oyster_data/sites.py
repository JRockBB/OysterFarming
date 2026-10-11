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
        "growing_area": None,  # TODO: DMR growing area code, from P1.5
        "notes": "forcing: lobo:node0065 (Maine EPSCoR LOBO, http://maine.loboviz.com). "
                 "Position: median of the buoy's own GPS fixes 2015-2023.",
    },
    {
        "site_id": "dam_upper",
        "name": "Upper Damariscotta Estuary, Mook Sea Farm/Perkins Point (LOBO node 66)",
        "lat": 43.99956,
        "lon": -69.54104,
        "region": "maine",
        "growing_area": None,  # TODO: DMR growing area code, from P1.5
        "notes": "forcing: lobo:node0066 (Maine EPSCoR LOBO, http://maine.loboviz.com). "
                 "Position: median of the buoy's own GPS fixes 2015-2022.",
    },
]


def sites_table() -> pa.Table:
    schema = arrow_schema("sites")
    return pa.Table.from_pylist(SITES, schema=schema)


def write_sites(dest: Path = PROCESSED_DIR / "sites.arrow") -> Path:
    write_contract(sites_table(), dest, "sites")
    return dest
