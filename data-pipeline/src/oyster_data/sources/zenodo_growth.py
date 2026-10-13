"""Eastern oyster growth observations from the Chesapeake Oyster BMP compilation (PLAN.md P1.9).

Source: "Morphometrics and nutrient concentration of eastern oysters (Crassostrea
virginica) sampled from Chesapeake Bay, USA", Zenodo, doi:10.5281/zenodo.11205266,
CC-BY-4.0. File Data_CBay_OysterBMP.csv: one row per individual oyster, with
shell height (mm) and tissue dry weight (g).

Only identified cohorts (one stock deployed together, sampled repeatedly) with
coordinates are used; COHORTS lists them. Individuals are aggregated to the
mean per sampling date, with ``n`` the number measured. The contract has no
ploidy column, so ploidy is part of ``cohort_id`` (``_2n`` / ``_3n``).

Known issues in the source file, recorded rather than changed:
- Date_Oysters_Deployed is stored as Excel serial day numbers (38653 = 2005-10-27).
- Kingsley-Smith York River rows give deploy serial 42309 (2015-10-31), after
  their own 2005-2006 sampling dates, so it is a data-entry error. The Patuxent
  and Severn deployments are 2005-10-27 and 2005-10-26. Kingsley-Smith et al.
  (2009, J. Shellfish Res. 28:169) describe deployment in 2004; the file dates
  are kept as given and the first sampling date is used as the model start.
"""

from __future__ import annotations

import io
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import pyarrow as pa

from oyster_data.contract import SPECS, arrow_schema, validate_table
from oyster_data.paths import RAW_DIR
from oyster_data.raw import http_get, save_raw

URL = "https://zenodo.org/api/records/11205266/files/Data_CBay_OysterBMP.csv/content"
RAW_FILE = "zenodo_growth/Data_CBay_OysterBMP.csv"
DATASET_ID = "zenodo11205266"


@dataclass(frozen=True)
class CohortSpec:
    cohort_id: str
    site_id: str
    data_source: str
    waterbody: str
    ploidy: str
    site_name: str | None = None   # lower-case Site column value, or None if blank in the file


COHORTS = [
    CohortSpec("ks2009_patuxent_3n", "ches_patuxent_ks", "Kingsley-Smith", "Patuxent River", "Triploid"),
    CohortSpec("ks2009_severn_3n", "ches_severn_ks", "Kingsley-Smith", "Severn River", "Triploid"),
    CohortSpec("ks2009_york_3n", "ches_york_ks", "Kingsley-Smith", "York River", "Triploid"),
]

# Not used: Parker-Bricker farm series (Chester River Orchard Point, Honga River,
# Potomac/Calvert, Chesapeake Bay "penny"). Their mean shell height falls between
# dates (Chester diploid 61 mm in Nov 2016 to 43 mm in Feb 2017), so they are not
# single cohorts followed through time (likely grading or harvest), and the Honga
# and Potomac rows have no coordinates. Liddel-2008 rows are wild reef samples,
# not cohorts.

def fetch(raw_dir: Path = RAW_DIR, get: Callable[[str], bytes] = http_get, overwrite: bool = False) -> Path:
    return save_raw(get(URL), raw_dir / RAW_FILE, overwrite=overwrite)


def read(text: str) -> pd.DataFrame:
    df = pd.read_csv(io.StringIO(text), low_memory=False)
    df["Site_norm"] = df["Site"].astype("string").str.strip().str.lower()
    return df


def _select(df: pd.DataFrame, c: CohortSpec) -> pd.DataFrame:
    sel = (df["Data_Source"] == c.data_source) & (df["Waterbody_Name"] == c.waterbody) & (df["Ploidy"] == c.ploidy)
    sel &= df["Site_norm"].isna() if c.site_name is None else (df["Site_norm"] == c.site_name)
    return df[sel]


def growth_obs(df: pd.DataFrame) -> pa.Table:
    """Contract growth_obs rows: cohort-date means of shell height and tissue dry weight."""
    rows = []
    for c in COHORTS:
        sub = _select(df, c)
        if sub.empty:
            raise ValueError(f"cohort {c.cohort_id} not found in source file")
        for d, g in sub.groupby("Date_Oysters_Removed"):
            sh = pd.to_numeric(g["Total_Shell_Height_mm"], errors="coerce").dropna()
            dw = pd.to_numeric(g["Tissue_Dry_Weight_g"], errors="coerce").dropna()
            rows.append({
                "dataset_id": f"{DATASET_ID}:{c.data_source}",
                "species": "C_virginica",
                "site_id": c.site_id,
                "cohort_id": c.cohort_id,
                "date": pd.Timestamp(d).date(),
                "shell_height_mm": float(sh.mean()) if len(sh) else None,
                "dry_tissue_g": float(dw.mean()) if len(dw) else None,
                "n": int(max(len(sh), len(dw))),
            })
    rows.sort(key=lambda r: (r["cohort_id"], r["date"]))
    table = pa.Table.from_pylist(rows, schema=arrow_schema("growth_obs"))
    validate_table(table, SPECS["growth_obs"])
    return table


def load(raw_dir: Path = RAW_DIR) -> pd.DataFrame:
    return read((raw_dir / RAW_FILE).read_text(encoding="utf-8"))
