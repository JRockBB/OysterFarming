"""Command-line entry point: ``uv run oyster-data <command>``.

Commands
  sites                      write data/processed/sites.arrow
  lobo-fetch --node N --years 2015-2023 [--overwrite]
                             download LOBO raw files to data/raw/lobo/
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
from oyster_data.sources import lobo


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


def build_forcing(raw_dir: Path = RAW_DIR, out_dir: Path = PROCESSED_DIR) -> tuple[Path, Path]:
    frames = clean_all_sources(raw_dir)
    if not frames:
        raise SystemExit("no raw forcing files found; run a fetch command first")
    merged = gapfill.merge_sources(frames)
    filled, reports = gapfill.fill_all(merged)
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
    args = parser.parse_args(argv)

    if args.command == "sites":
        print(f"wrote {write_sites()}")
    elif args.command == "lobo-fetch":
        return cmd_lobo_fetch(args)
    elif args.command == "forcing":
        out, report = build_forcing()
        print(f"wrote {out}\nwrote {report}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
