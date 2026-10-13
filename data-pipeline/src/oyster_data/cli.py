"""Command-line entry point: ``uv run oyster-data <command>``.

Commands
  sites                      write data/processed/sites.arrow
  lobo-fetch --node N --years 2015-2023 [--overwrite]
                             download LOBO raw files to data/raw/lobo/
  economics-fetch [--overwrite]
                             download NOAA FOSS Maine oyster landings to data/raw/foss/
  economics                  write data/processed/economics.arrow from the raw FOSS file
  secondary-fetch            download Maine secondary stations used for gap filling
                             (Boothbay DMR SST, DMC dock and VIIRS chlorophyll, NERACOOS E01)
  disease-fetch [--overwrite] download the MD DNR 2024 Fall Survey report (Dermo/MSX tables)
  disease                    write data/processed/disease_obs.arrow
  growth-fetch [--overwrite]  download the Chesapeake Oyster BMP growth compilation (Zenodo)
  growth                     write data/processed/growth_obs.arrow
  chesapeake-fetch [--overwrite]
                             download CO-OPS, CBP DataHub and NERRS data for the Chesapeake growth sites
  forcing                    clean all raw forcing sources, merge, gap-fill and write
                             data/processed/forcing_daily.arrow and forcing_report.md
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from oyster_data import gapfill
from oyster_data.forcing import to_forcing_table
from oyster_data.contract import write_contract
from oyster_data.paths import PROCESSED_DIR, RAW_DIR
from oyster_data.sites import write_sites
from datetime import date

from oyster_data.sources import (boothbay_dmr, cbp_datahub, chesapeake_sites, coops, dmc_dock, foss, lobo,
                                 mddnr_fall_survey, nerrs_ncei, neracoos, viirs_chl, zenodo_growth)


def _years(spec: str) -> list[int]:
    if "-" in spec:
        a, b = spec.split("-", 1)
        return list(range(int(a), int(b) + 1))
    return [int(y) for y in spec.split(",")]


def cmd_lobo_fetch(args: argparse.Namespace) -> int:
    for year in _years(args.years):
        try:
            path = lobo.fetch_year(args.node, year, overwrite=args.overwrite)
            print(f"saved {path}")
        except FileExistsError as e:
            print(f"skipped: {e}")
    return 0


def clean_all_sources(raw_dir: Path = RAW_DIR) -> list[pd.DataFrame]:
    """Daily rows from every source with raw files present, in merge priority order."""
    frames = []
    for node in sorted(lobo.NODES):
        if list((raw_dir / lobo.RAW_SUBDIR).glob(f"node{node:04d}_*.tsv")):
            frames.append(lobo.clean_node(node, raw_dir))
    return frames


SECONDARY_YEARS = (2014, 2023)


def cmd_secondary_fetch(overwrite: bool) -> int:
    first, last = SECONDARY_YEARS
    jobs = [lambda: boothbay_dmr.fetch(date(first, 1, 1), date(last + 1, 1, 1), overwrite=overwrite),
            lambda: dmc_dock.fetch(overwrite=overwrite)]
    for y in range(first, last + 1):
        jobs += [(lambda y=y, d=d: neracoos.fetch_year(y, d, overwrite=overwrite)) for d in neracoos.DATASETS]
        jobs.append(lambda y=y: viirs_chl.fetch_year(y, overwrite=overwrite))
    failures = 0
    for job in jobs:
        try:
            path = job()
            print(f"saved {path}" if path is not None else "no data for one requested year (HTTP 404)")
        except FileExistsError as e:
            print(f"skipped: {e}")
        except Exception as e:  # keep going; report every failure
            failures += 1
            print(f"FAILED: {e}")
    return 1 if failures else 0


def cmd_chesapeake_fetch(overwrite: bool) -> int:
    jobs = []
    for site in chesapeake_sites.SITES:
        first, last = int(site.start[:4]), int(site.end[:4])
        jobs += [(lambda st=site.coops_station, y=y: coops.fetch_year(st, y, overwrite=overwrite))
                 for y in range(first, last + 1)]
        for st in filter(None, (site.cbp_station, site.cbp_secondary)):
            jobs.append(lambda st=st: cbp_datahub.fetch(st, date(2005, 9, 1), date(2007, 11, 30),
                                                        overwrite=overwrite))
        if site.nerrs_nut_station:
            jobs += [(lambda st=site.nerrs_nut_station, y=y: nerrs_ncei.fetch_year(st, y, overwrite=overwrite))
                     for y in range(first, last + 1)]
        if site.nerrs_temp_station:
            jobs += [(lambda st=site.nerrs_temp_station, y=y: nerrs_ncei.fetch_year(st, y, overwrite=overwrite))
                     for y in range(first, last + 1)]
    failures = 0
    for job in jobs:
        try:
            print(f"saved {job()}")
        except FileExistsError as e:
            print(f"skipped: {e}")
        except Exception as e:  # keep going; report every failure
            failures += 1
            print(f"FAILED: {e}")
    return 1 if failures else 0


def load_secondaries(raw_dir: Path = RAW_DIR) -> dict[str, list[gapfill.Secondary]]:
    """Secondary stations with raw files present, keyed by forcing variable (priority order)."""
    out: dict[str, list[gapfill.Secondary]] = {}
    loaders = [
        ("temp_c", boothbay_dmr.STATION, boothbay_dmr.load),
        ("temp_c", neracoos.station_name("temperature", 1), lambda d: neracoos.load("temperature", 1, d)),
        ("salinity_psu", neracoos.station_name("salinity", 1), lambda d: neracoos.load("salinity", 1, d)),
        ("salinity_psu", neracoos.station_name("salinity", 10), lambda d: neracoos.load("salinity", 10, d)),
        ("chl_ug_l", dmc_dock.STATION, dmc_dock.load),
        ("chl_ug_l", viirs_chl.STATION, viirs_chl.load),
    ]
    for var, name, loader in loaders:
        try:
            series = loader(raw_dir)
        except FileNotFoundError:
            continue
        out.setdefault(var, []).append((name, series))
    return out


def build_forcing(raw_dir: Path = RAW_DIR, out_dir: Path = PROCESSED_DIR) -> tuple[Path, Path]:
    frames = clean_all_sources(raw_dir)
    if not frames:
        raise SystemExit("no raw forcing files found; run a fetch command first")
    merged = gapfill.merge_sources(frames)
    filled, reports = gapfill.fill_all(merged, load_secondaries(raw_dir))
    try:
        ches, ches_reports = chesapeake_sites.build_all(raw_dir)
    except FileNotFoundError as err:
        print(f"Chesapeake sites skipped (run chesapeake-fetch): {err}")
    else:
        filled = pd.concat([filled, ches], ignore_index=True)
        reports += ches_reports
    table = to_forcing_table(filled)
    out = out_dir / "forcing_daily.arrow"
    write_contract(table, out, "forcing_daily")
    report = out_dir / "forcing_report.md"
    report.write_text(gapfill.render_report(reports), encoding="utf-8", newline="\n")
    return out, report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="oyster-data")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("sites", help="write sites.arrow")
    p = sub.add_parser("lobo-fetch", help="download LOBO buoy data")
    p.add_argument("--node", type=int, required=True, choices=sorted(lobo.NODES))
    p.add_argument("--years", required=True, help="e.g. 2015-2023 or 2016,2018")
    p.add_argument("--overwrite", action="store_true")
    sub.add_parser("forcing", help="clean, merge and gap-fill forcing")
    df_ = sub.add_parser("disease-fetch", help="download MD DNR Fall Survey report")
    df_.add_argument("--overwrite", action="store_true")
    sub.add_parser("disease", help="write disease_obs.arrow")
    gf = sub.add_parser("growth-fetch", help="download Zenodo growth data")
    gf.add_argument("--overwrite", action="store_true")
    sub.add_parser("growth", help="write growth_obs.arrow")
    cf = sub.add_parser("chesapeake-fetch", help="download Chesapeake site forcing")
    cf.add_argument("--overwrite", action="store_true")
    sf = sub.add_parser("secondary-fetch", help="download secondary stations for gap filling")
    sf.add_argument("--overwrite", action="store_true")
    e = sub.add_parser("economics-fetch", help="download NOAA FOSS landings")
    e.add_argument("--overwrite", action="store_true")
    sub.add_parser("economics", help="write economics.arrow")
    args = parser.parse_args(argv)

    if args.command == "sites":
        print(f"wrote {write_sites()}")
    elif args.command == "lobo-fetch":
        return cmd_lobo_fetch(args)
    elif args.command == "economics-fetch":
        try:
            print(f"saved {foss.fetch(overwrite=args.overwrite)}")
        except FileExistsError as err:
            print(f"skipped: {err}")
    elif args.command == "economics":
        out = PROCESSED_DIR / "economics.arrow"
        write_contract(foss.clean(), out, "economics")
        print(f"wrote {out}")
    elif args.command == "disease-fetch":
        try:
            print(f"saved {mddnr_fall_survey.fetch(overwrite=args.overwrite)}")
        except FileExistsError as err:
            print(f"skipped: {err}")
    elif args.command == "disease":
        out = PROCESSED_DIR / "disease_obs.arrow"
        write_contract(mddnr_fall_survey.disease_obs(*mddnr_fall_survey.load_pages()), out, "disease_obs")
        print(f"wrote {out}")
    elif args.command == "growth-fetch":
        try:
            print(f"saved {zenodo_growth.fetch(overwrite=args.overwrite)}")
        except FileExistsError as err:
            print(f"skipped: {err}")
    elif args.command == "growth":
        out = PROCESSED_DIR / "growth_obs.arrow"
        write_contract(zenodo_growth.growth_obs(zenodo_growth.load()), out, "growth_obs")
        print(f"wrote {out}")
    elif args.command == "chesapeake-fetch":
        return cmd_chesapeake_fetch(args.overwrite)
    elif args.command == "secondary-fetch":
        return cmd_secondary_fetch(args.overwrite)
    elif args.command == "forcing":
        out, report = build_forcing()
        print(f"wrote {out}\nwrote {report}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
