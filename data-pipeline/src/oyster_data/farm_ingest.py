"""Convert filled farm templates (docs/farm_templates/) into contract tables (PLAN.md P6.2).

Outputs, in the chosen directory:
- contract tables: ``sites.arrow`` (one row), ``growth_obs.arrow`` (cohort-date means with ``n``),
  ``mortality_obs.arrow`` (one row per count, method ``farm_count:<cohort>``);
- non-contract tables for later tasks: ``farm_temperature_daily.csv`` (daily mean logger temperature,
  days with readings in at least ``MIN_HOURS_PER_DAY`` distinct hours), ``seeding.csv``, ``harvest.csv``.

Anonymity: the site id is ``farm_`` plus the first 8 hex digits of SHA-256 of farm name and site label.
Farm and site names are not written anywhere. Cohort ids are prefixed with the site id.
Species is C_virginica (Maine farms grow eastern oysters).
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd
import pyarrow as pa

from oyster_data.contract import SPECS, arrow_schema, validate_table

TEMPLATES = {
    "farm_info": ["farm_name", "site_label", "lat", "lon", "growing_area", "notes"],
    "growth_checks": ["date", "cohort", "oyster_no", "shell_height_mm", "dry_tissue_g"],
    "mortality_counts": ["date", "cohort", "live_count", "dead_count", "days_since_last_count"],
    "seeding": ["date", "cohort", "number_seeded", "seed_shell_height_mm", "ploidy", "gear", "seed_source"],
    "harvest": ["date", "cohort", "number_harvested", "mean_shell_height_mm", "price_per_oyster_usd"],
    "temperature_logger": ["datetime_local", "temp_c", "depth_m"],
}
REQUIRED = {"farm_info", "growth_checks"}
MIN_HOURS_PER_DAY = 12


class FarmDataError(ValueError):
    pass


def site_code(farm_name: str, site_label: str) -> str:
    key = f"{farm_name.strip().lower()}|{site_label.strip().lower()}".encode("utf-8")
    return "farm_" + hashlib.sha256(key).hexdigest()[:8]


def read_templates(folder: Path) -> dict[str, pd.DataFrame]:
    out = {}
    for name, cols in TEMPLATES.items():
        path = Path(folder) / f"{name}.csv"
        if not path.exists():
            if name in REQUIRED:
                raise FarmDataError(f"missing required file {path.name}")
            continue
        df = pd.read_csv(path, dtype=str, keep_default_na=False).apply(lambda c: c.str.strip())
        if list(df.columns) != cols:
            raise FarmDataError(f"{path.name}: columns must be {cols}, got {list(df.columns)}")
        out[name] = df.replace("", pd.NA)
    return out


def _rows(mask: pd.Series) -> list[int]:
    """Spreadsheet row numbers (header is row 1)."""
    return [int(i) + 2 for i in mask[mask].index]


def _dates(s: pd.Series, file: str) -> pd.Series:
    d = pd.to_datetime(s, format="%Y-%m-%d", errors="coerce")
    if d.isna().any():
        raise FarmDataError(f"{file}: dates must be YEAR-MONTH-DAY, bad rows {_rows(d.isna())}")
    return d.dt.date


def _num(s: pd.Series, file: str, col: str, nullable: bool = False, nonneg: bool = True) -> pd.Series:
    v = pd.to_numeric(s, errors="coerce")
    bad = (v.isna() & s.notna()) if nullable else v.isna()
    if bad.any():
        raise FarmDataError(f"{file}: {col} must be a number, bad rows {_rows(bad)}")
    if nonneg and (v < 0).any():
        raise FarmDataError(f"{file}: {col} must not be negative, bad rows {_rows(v < 0)}")
    return v


def site_table(info: pd.DataFrame) -> tuple[pa.Table, str]:
    f = "farm_info.csv"
    if len(info) != 1:
        raise FarmDataError(f"{f} must have exactly one row (one lease site per folder)")
    r = info.iloc[0]
    if pd.isna(r["farm_name"]) or pd.isna(r["site_label"]):
        raise FarmDataError(f"{f}: farm_name and site_label are required")
    code = site_code(r["farm_name"], r["site_label"])
    lat = float(_num(info["lat"], f, "lat", nonneg=False).iloc[0])
    lon = float(_num(info["lon"], f, "lon", nonneg=False).iloc[0])
    row = {"site_id": code, "name": code, "lat": lat, "lon": lon, "region": "maine",
           "growing_area": None if pd.isna(r["growing_area"]) else str(r["growing_area"]),
           "notes": "farm-supplied site (P6.2); names removed"}
    t = pa.Table.from_pylist([row], schema=arrow_schema("sites"))
    validate_table(t, SPECS["sites"])
    return t, code


def growth_table(g: pd.DataFrame, code: str) -> pa.Table:
    f = "growth_checks.csv"
    if g["cohort"].isna().any():
        raise FarmDataError(f"{f}: every row needs a cohort, bad rows {_rows(g['cohort'].isna())}")
    df = pd.DataFrame({"date": _dates(g["date"], f), "cohort": g["cohort"],
                       "sh": _num(g["shell_height_mm"], f, "shell_height_mm", nullable=True),
                       "dw": _num(g["dry_tissue_g"], f, "dry_tissue_g", nullable=True)})
    agg = (df.groupby(["cohort", "date"], sort=True)
             .agg(sh=("sh", "mean"), dw=("dw", "mean"), n=("sh", "count")).reset_index())
    rows = [{"dataset_id": f"farm:{code}", "species": "C_virginica", "site_id": code,
             "cohort_id": f"{code}_{r.cohort}", "date": r.date,
             "shell_height_mm": None if pd.isna(r.sh) else float(r.sh),
             "dry_tissue_g": None if pd.isna(r.dw) else float(r.dw), "n": int(r.n)} for r in agg.itertuples()]
    t = pa.Table.from_pylist(rows, schema=arrow_schema("growth_obs"))
    validate_table(t, SPECS["growth_obs"])
    return t


def mortality_table(m: pd.DataFrame, code: str) -> pa.Table:
    f = "mortality_counts.csv"
    d = _dates(m["date"], f)
    live, dead = _num(m["live_count"], f, "live_count"), _num(m["dead_count"], f, "dead_count")
    days = _num(m["days_since_last_count"], f, "days_since_last_count", nullable=True)
    if ((live + dead) == 0).any():
        raise FarmDataError(f"{f}: no oysters counted, bad rows {_rows((live + dead) == 0)}")
    rows = [{"site_id": code, "date": d.iloc[i],
             "period_days": None if pd.isna(days.iloc[i]) else int(days.iloc[i]),
             "mortality_fraction": float(dead.iloc[i] / (live.iloc[i] + dead.iloc[i])),
             "method": f"farm_count:{m['cohort'].iloc[i]}", "source": f"farm:{code}"} for i in range(len(m))]
    t = pa.Table.from_pylist(rows, schema=arrow_schema("mortality_obs"))
    validate_table(t, SPECS["mortality_obs"])
    return t


def logger_daily(t: pd.DataFrame) -> pd.DataFrame:
    f = "temperature_logger.csv"
    ts = pd.to_datetime(t["datetime_local"], format="%Y-%m-%d %H:%M", errors="coerce")
    if ts.isna().any():
        raise FarmDataError(f"{f}: datetime_local must be YEAR-MONTH-DAY HH:MM, bad rows {_rows(ts.isna())}")
    temp = _num(t["temp_c"], f, "temp_c", nonneg=False)
    df = pd.DataFrame({"date": ts.dt.date, "hour": ts.dt.hour, "temp_c": temp})
    g = df.groupby("date").agg(temp_c=("temp_c", "mean"), hours=("hour", "nunique")).reset_index()
    return g[g["hours"] >= MIN_HOURS_PER_DAY][["date", "temp_c"]].reset_index(drop=True)


def _write_arrow(path: Path, table: pa.Table) -> None:
    with pa.OSFile(str(path), "wb") as sink, pa.ipc.new_file(sink, table.schema) as w:
        w.write_table(table)


def ingest(folder: Path, out_dir: Path) -> dict[str, Path]:
    """Convert one farm folder; returns the files written by name."""
    tables = read_templates(Path(folder))
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sites, code = site_table(tables["farm_info"])
    written: dict[str, Path] = {}
    arrow_tables = {"sites": sites, "growth_obs": growth_table(tables["growth_checks"], code)}
    if "mortality_counts" in tables and len(tables["mortality_counts"]):
        arrow_tables["mortality_obs"] = mortality_table(tables["mortality_counts"], code)
    for name, table in arrow_tables.items():
        written[name] = out_dir / f"{name}.arrow"
        _write_arrow(written[name], table)
    if "temperature_logger" in tables and len(tables["temperature_logger"]):
        p = out_dir / "farm_temperature_daily.csv"
        logger_daily(tables["temperature_logger"]).assign(site_id=code).to_csv(p, index=False)
        written["farm_temperature_daily"] = p
    for name in ("seeding", "harvest"):
        if name in tables and len(tables[name]):
            df = tables[name].copy()
            df["date"] = _dates(df["date"], f"{name}.csv")
            df["cohort"] = code + "_" + df["cohort"]
            df.insert(0, "site_id", code)
            p = out_dir / f"{name}.csv"
            df.to_csv(p, index=False)
            written[name] = p
    return written
