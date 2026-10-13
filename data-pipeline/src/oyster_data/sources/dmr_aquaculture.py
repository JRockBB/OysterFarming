"""Maine DMR aquaculture landings of American/eastern oysters (PLAN.md P1.6, P5.1).

Source: Maine Department of Marine Resources, "Maine Aquaculture Landings Report"
(posted April 2026), https://www.maine.gov/dmr/sites/maine.gov.dmr/files/2026-04/MAINE%20AQUACULTURE%20LANDINGS%20REPORT.pdf
Aquaculture licence holders report trip-level data from 2021, so these are farmed
oysters only, unlike the NOAA FOSS series (oyster_data.sources.foss).

The PDF's text layer is unusable (characters extracted out of order), so the
oyster rows were transcribed from the rendered page: Table 1 (live pounds),
Table 1A (count) and Table 3 (value). `check_totals` verifies the transcription
against the report's own "Total" column.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import pyarrow as pa

from oyster_data.contract import SPECS, arrow_schema, validate_table
from oyster_data.paths import RAW_DIR
from oyster_data.raw import http_get, save_raw

URL = ("https://www.maine.gov/dmr/sites/maine.gov.dmr/files/2026-04/"
       "MAINE%20AQUACULTURE%20LANDINGS%20REPORT.pdf")
RAW_FILE = "dmr_aquaculture/maine_aquaculture_landings_report.pdf"
REGION = "maine_aquaculture"

# year: (live pounds [Table 1], count [Table 1A], value USD [Table 3])
OYSTER_ROWS = {
    2021: (6_204_768, 22_074_315, 13_148_506),
    2022: (5_154_434, 18_387_472, 13_612_203),
    2023: (6_349_604, 22_663_843, 17_489_690),
    2024: (6_560_902, 23_401_019, 18_217_846),
    2025: (7_095_903, 24_889_458, 18_265_588),
}
# The report's "Total" column for the oyster rows.
REPORT_TOTALS = {"pounds": 31_365_612, "count": 111_416_108, "value_usd": 80_733_833}


def fetch(raw_dir: Path = RAW_DIR, get: Callable[[str], bytes] = http_get, overwrite: bool = False) -> Path:
    body = get(URL)
    if not body.startswith(b"%PDF"):
        raise RuntimeError("DMR aquaculture landings download is not a PDF")
    return save_raw(body, raw_dir / RAW_FILE, overwrite=overwrite)


def check_totals() -> None:
    """Raise if the transcribed rows do not sum to the report's totals."""
    sums = {"pounds": sum(r[0] for r in OYSTER_ROWS.values()),
            "count": sum(r[1] for r in OYSTER_ROWS.values()),
            "value_usd": sum(r[2] for r in OYSTER_ROWS.values())}
    # The report shows rounded yearly values but totals the unrounded ones, so each
    # sum may differ from the printed total by up to 0.5 per year (pounds and count
    # differ by 1 for 2021-2025; value matches exactly). Any transcription slip in a
    # digit is far larger than this.
    tolerance = 0.5 * len(OYSTER_ROWS)
    bad = {k: (sums[k], REPORT_TOTALS[k]) for k in sums if abs(sums[k] - REPORT_TOTALS[k]) > tolerance}
    if bad:
        raise ValueError(f"transcribed oyster rows do not match report totals: {bad}")


def to_economics() -> pa.Table:
    check_totals()
    rows = [{"year": y, "region": REGION, "value_usd": float(v), "pounds": float(lb), "count": float(n),
             "source": "maine_dmr:aquaculture_landings_report_2026-04:Tables1,1A,3"}
            for y, (lb, n, v) in sorted(OYSTER_ROWS.items())]
    table = pa.Table.from_pylist(rows, schema=arrow_schema("economics"))
    validate_table(table, SPECS["economics"])
    return table
