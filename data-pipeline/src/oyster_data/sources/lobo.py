"""Maine EPSCoR LOBO buoys in the Damariscotta River estuary (UMaine).

Data service: LOBOViz CGI, documented at http://maine.loboviz.com/cgidoc/
Request: ``/cgi-data/nph-data.cgi?min_date=YYYYMMDD&max_date=YYYYMMDD&y=<vars>&data_format=text&node=<n>``
Response: two title lines, a tab-separated header ``date [Eastern]\\t<name> [<unit>]...``,
then one row per hourly record. The sensors carry no QC flags.

Coverage found when this module was written (hourly records per year): node 65
2015-2016 and 2021-2023; node 66 2015-2019 and 2021-2022. Deployments are
seasonal, so winters are mostly missing.
"""

from __future__ import annotations

import io
from dataclasses import dataclass
from datetime import date
from pathlib import Path

import pandas as pd

from oyster_data.forcing import hourly_to_daily
from oyster_data.paths import RAW_DIR
from oyster_data.raw import Getter, http_get, save_raw

BASE_URL = "http://maine.loboviz.com/cgi-data/nph-data.cgi"
VARIABLES = {"temperature": "temp_c", "salinity": "salinity_psu", "chlorophyll": "chl_ug_l"}
RAW_SUBDIR = "lobo"

# Times are US Eastern local time WITH daylight saving: raw files repeat the
# 01:00 hour on each November DST change (checked 2016-2022 on both nodes).
LOCAL_TZ = "America/New_York"


# Records with salinity below this are treated as out of the water (on deck or in
# a freshwater bucket before deployment) and dropped entirely, temperature
# included. Evidence: node 65 on 2016-04-18/19 reads 0.11 psu with 15-18 degC
# before the sonde went in at ~30 psu and 7 degC. Both nodes are marine sites
# whose valid daily means never fall below ~21 psu.
OUT_OF_WATER_SALINITY_PSU = 2.0


@dataclass(frozen=True)
class LoboNode:
    node: int
    site_id: str
    name: str


NODES = {
    65: LoboNode(65, "dam_lower", "Lower Damariscotta (LOBO node 65)"),
    66: LoboNode(66, "dam_upper", "Upper Damariscotta Estuary, Mook Sea Farm/Perkins Point (LOBO node 66)"),
}


def build_url(node: int, start: date, end: date, variables: tuple[str, ...] = tuple(VARIABLES)) -> str:
    return (f"{BASE_URL}?min_date={start:%Y%m%d}&max_date={end:%Y%m%d}"
            f"&y={','.join(variables)}&data_format=text&node={node}")


def raw_path(node: int, year: int, raw_dir: Path = RAW_DIR) -> Path:
    return raw_dir / RAW_SUBDIR / f"node{node:04d}_{year}.tsv"


def fetch_year(node: int, year: int, raw_dir: Path = RAW_DIR, get: Getter = http_get,
               overwrite: bool = False) -> Path:
    """Download one calendar year for ``node`` to ``data/raw/lobo/``."""
    url = build_url(node, date(year, 1, 1), date(year, 12, 31))
    return save_raw(get(url), raw_path(node, year, raw_dir), overwrite=overwrite)


def parse(text: str) -> pd.DataFrame:
    """Parse a LOBOViz text response into ``time_utc`` plus contract-named columns."""
    lines = text.splitlines()
    header_idx = next((i for i, line in enumerate(lines) if line.startswith("date [")), None)
    if header_idx is None:
        raise ValueError("LOBO response has no 'date [...]' header line")
    header = lines[header_idx].split("\t")
    if not header[0].startswith("date [Eastern]"):
        raise ValueError(f"unexpected LOBO time column {header[0]!r}")
    names = ["local_time"]
    for col in header[1:]:
        key = col.split(" [")[0].strip()
        if key not in VARIABLES:
            raise ValueError(f"unexpected LOBO column {col!r}")
        names.append(VARIABLES[key])
    body = "\n".join(lines[header_idx + 1:])
    df = pd.read_csv(io.StringIO(body), sep="\t", header=None, names=names)
    local = pd.to_datetime(df.pop("local_time"), format="%Y-%m-%d %H:%M:%S")
    try:
        aware = local.dt.tz_localize(LOCAL_TZ, ambiguous="infer", nonexistent="NaT")
    except Exception:  # repeated hour not inferable (e.g. one copy missing): drop ambiguous records
        aware = local.dt.tz_localize(LOCAL_TZ, ambiguous="NaT", nonexistent="NaT")
    df.insert(0, "time_utc", aware.dt.tz_convert("UTC"))
    df = df[df["time_utc"].notna()].reset_index(drop=True)
    for col in VARIABLES.values():
        if col not in df:
            df[col] = float("nan")
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def drop_out_of_water(hourly: pd.DataFrame) -> pd.DataFrame:
    """Remove records taken with the sonde out of the water (see OUT_OF_WATER_SALINITY_PSU)."""
    keep = ~(hourly["salinity_psu"] < OUT_OF_WATER_SALINITY_PSU)
    return hourly[keep].reset_index(drop=True)


def clean_node(node: int, raw_dir: Path = RAW_DIR) -> pd.DataFrame:
    """Daily ``forcing_daily`` rows from every raw year file saved for ``node``."""
    info = NODES[node]
    files = sorted((raw_dir / RAW_SUBDIR).glob(f"node{node:04d}_*.tsv"))
    if not files:
        raise FileNotFoundError(f"no raw LOBO files for node {node} in {raw_dir / RAW_SUBDIR}")
    hourly = pd.concat([parse(f.read_text(encoding="utf-8")) for f in files], ignore_index=True)
    hourly = hourly.drop_duplicates(subset="time_utc").sort_values("time_utc")
    hourly = drop_out_of_water(hourly)
    return hourly_to_daily(hourly, site_id=info.site_id, source=f"lobo:node{node:04d}")
