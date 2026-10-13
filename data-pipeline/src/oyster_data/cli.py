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

from oyster_data.sources import boothbay_dmr, dmc_dock, foss, lobo, neracoos, viirs_chl


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
    elif args.command == "secondary-fetch":
        return cmd_secondary_fetch(args.overwrite)
    elif args.command == "forcing":
        out, report = build_forcing()
        print(f"wrote {out}\nwrote {report}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
