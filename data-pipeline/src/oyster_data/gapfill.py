"""Merge forcing sources per site and fill gaps (PLAN.md P1.4).

Rules (processing choices, recorded here and in forcing_report.md):

1. Merge: when several sources cover the same site-day, the first source in the
   given priority order wins, variable by variable.
2. Short gaps: a run of at most ``MAX_LINEAR_GAP_DAYS`` missing days with
   observations on both sides is filled by linear interpolation.
3. Long gaps: filled from the same site's day-of-year climatology, the mean of
   observed values within +/- ``CLIM_HALF_WINDOW_DAYS`` days of that day of year
   (circular across the new year), pooled over all years.
4. Days whose day-of-year window has no observations at all stay missing. Rows
   still lacking temperature or salinity are dropped (not nullable in the contract).
5. Filled rows get ``qc_flag = 1`` and ``;fill=linear`` or ``;fill=climatology``
   appended to ``source`` (climatology wins if both methods were used in a row).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

MAX_LINEAR_GAP_DAYS = 7
CLIM_HALF_WINDOW_DAYS = 7
FILL_VARIABLES = ("temp_c", "salinity_psu", "chl_ug_l")

OBSERVED, LINEAR, CLIMATOLOGY, MISSING = "observed", "linear", "climatology", "missing"


def merge_sources(frames: list[pd.DataFrame]) -> pd.DataFrame:
    """Combine daily forcing frames for one or more sites, earlier frames taking priority."""
    if not frames:
        raise ValueError("no forcing frames to merge")
    merged: pd.DataFrame | None = None
    for df in frames:
        d = df.copy()
        d["date"] = pd.to_datetime(d["date"]).dt.date
        d = d.set_index(["site_id", "date"])
        merged = d if merged is None else merged.combine_first(d)
    assert merged is not None
    return merged.reset_index().sort_values(["site_id", "date"]).reset_index(drop=True)


def _doy(dates: pd.DatetimeIndex) -> np.ndarray:
    """Day of year on a 365-day calendar (29 Feb shares 28 Feb's index)."""
    doy = dates.dayofyear.to_numpy().copy()
    leap_after_feb = dates.is_leap_year & (dates.month > 2)
    doy[leap_after_feb] -= 1
    return doy


def climatology(series: pd.Series, half_window: int = CLIM_HALF_WINDOW_DAYS) -> np.ndarray:
    """Return a length-366 array (index = day of year 1..365) of windowed means of observed values."""
    obs = series.dropna()
    clim = np.full(366, np.nan)
    if obs.empty:
        return clim
    doy = _doy(pd.DatetimeIndex(obs.index))
    sums = np.bincount(doy, weights=obs.to_numpy(), minlength=366)[1:]
    counts = np.bincount(doy, minlength=366)[1:]
    for d in range(365):
        idx = [(d + k) % 365 for k in range(-half_window, half_window + 1)]
        n = counts[idx].sum()
        if n > 0:
            clim[d + 1] = sums[idx].sum() / n
    return clim


def fill_series(series: pd.Series, max_linear: int = MAX_LINEAR_GAP_DAYS,
                half_window: int = CLIM_HALF_WINDOW_DAYS) -> tuple[pd.Series, pd.Series]:
    """Fill one daily series indexed by a complete DatetimeIndex.

    Returns ``(filled_values, method)`` where ``method`` is one of
    observed / linear / climatology / missing for each day.
    """
    values = series.astype(float).copy()
    method = pd.Series(np.where(values.notna(), OBSERVED, MISSING), index=values.index, dtype=object)
    isna = values.isna().to_numpy()
    n = len(values)
    i = 0
    interp = values.interpolate(method="linear", limit_area="inside")
    while i < n:
        if not isna[i]:
            i += 1
            continue
        j = i
        while j < n and isna[j]:
            j += 1
        run = j - i
        bounded = i > 0 and j < n
        if bounded and run <= max_linear:
            values.iloc[i:j] = interp.iloc[i:j]
            method.iloc[i:j] = LINEAR
        i = j
    clim = climatology(series, half_window)
    still = values.isna()
    if still.any():
        doy = _doy(pd.DatetimeIndex(values.index[still]))
        cvals = clim[doy]
        ok = ~np.isnan(cvals)
        idx = values.index[still][ok]
        values.loc[idx] = cvals[ok]
        method.loc[idx] = CLIMATOLOGY
    return values, method


def longest_run(mask: pd.Series) -> int:
    """Length of the longest run of True values."""
    best = cur = 0
    for v in mask.to_numpy():
        cur = cur + 1 if v else 0
        best = max(best, cur)
    return best


@dataclass
class SiteReport:
    site_id: str
    first_date: object
    last_date: object
    n_days: int
    counts: dict[str, dict[str, int]] = field(default_factory=dict)
    longest_climatology_run: dict[str, int] = field(default_factory=dict)
    rows_written: int = 0


def fill_site(daily: pd.DataFrame, site_id: str) -> tuple[pd.DataFrame, SiteReport]:
    """Gap-fill one site's merged daily rows over its first-to-last observed date range."""
    d = daily[daily["site_id"] == site_id].copy()
    if d.empty:
        raise ValueError(f"no rows for site {site_id!r}")
    d["date"] = pd.to_datetime(d["date"])
    d = d.set_index("date").sort_index()
    full = pd.date_range(d.index.min(), d.index.max(), freq="D")
    d = d.reindex(full)
    report = SiteReport(site_id, full[0].date(), full[-1].date(), len(full))
    methods = {}
    for var in FILL_VARIABLES:
        vals, meth = fill_series(d[var])
        d[var] = vals
        methods[var] = meth
        report.counts[var] = {m: int((meth == m).sum()) for m in (OBSERVED, LINEAR, CLIMATOLOGY, MISSING)}
        report.longest_climatology_run[var] = longest_run(meth == CLIMATOLOGY)
    base_source = d["source"].ffill().bfill()
    any_clim = np.zeros(len(d), dtype=bool)
    any_lin = np.zeros(len(d), dtype=bool)
    for var in FILL_VARIABLES:
        any_clim |= (methods[var] == CLIMATOLOGY).to_numpy()
        any_lin |= (methods[var] == LINEAR).to_numpy()
    suffix = np.where(any_clim, ";fill=climatology", np.where(any_lin, ";fill=linear", ""))
    d["source"] = base_source.str.split(";").str[0] + suffix
    prior_qc = d["qc_flag"].fillna(0).astype(int)
    d["qc_flag"] = np.where(any_clim | any_lin, np.maximum(prior_qc, 1), prior_qc).astype(int)
    d.loc[prior_qc == 2, "qc_flag"] = 2
    d["site_id"] = site_id
    d = d.dropna(subset=["temp_c", "salinity_psu"])
    out = d.reset_index(names="date")
    out["date"] = out["date"].dt.date
    report.rows_written = len(out)
    cols = ["site_id", "date", "temp_c", "salinity_psu", "chl_ug_l", "pom_mg_l", "source", "qc_flag"]
    return out[cols], report


def fill_all(daily: pd.DataFrame) -> tuple[pd.DataFrame, list[SiteReport]]:
    parts, reports = [], []
    for site in sorted(daily["site_id"].unique()):
        filled, rep = fill_site(daily, site)
        parts.append(filled)
        reports.append(rep)
    return pd.concat(parts, ignore_index=True), reports


def render_report(reports: list[SiteReport]) -> str:
    """Markdown coverage report (data/processed/forcing_report.md)."""
    lines = [
        "# Forcing coverage report",
        "",
        "Generated by `oyster-data forcing` (PLAN.md P1.4). Percentages are of days between each",
        "site's first and last observation.",
        "",
        f"Rules: linear interpolation for gaps of at most {MAX_LINEAR_GAP_DAYS} days; otherwise the",
        f"site's day-of-year climatology (+/- {CLIM_HALF_WINDOW_DAYS} day window, all years pooled).",
        "Days with no climatology stay missing; rows missing temperature or salinity are not written.",
        "",
    ]
    for r in reports:
        lines += [f"## {r.site_id}", "",
                  f"Range {r.first_date} to {r.last_date}: {r.n_days} days, {r.rows_written} rows written.", "",
                  "| variable | observed % | linear % | climatology % | missing % | longest climatology run (days) |",
                  "|---|---|---|---|---|---|"]
        for var, c in r.counts.items():
            pct = {k: 100 * v / r.n_days for k, v in c.items()}
            lines.append(f"| {var} | {pct[OBSERVED]:.1f} | {pct[LINEAR]:.1f} | {pct[CLIMATOLOGY]:.1f} | "
                         f"{pct[MISSING]:.1f} | {r.longest_climatology_run[var]} |")
        lines.append("")
    return "\n".join(lines)
