"""Tests for the MD DNR Fall Survey disease extraction.

Fixtures are the text pdfplumber extracts from three pages of the real 2024
report. Expected values were read by hand from the report tables.
"""

from datetime import date
from pathlib import Path

import pytest

from oyster_data.sources import mddnr_fall_survey as m

FIX = Path(__file__).parent / "fixtures"


def page(name: str) -> str:
    return (FIX / name).read_text(encoding="utf-8")


def lookup(rows, site, year, pathogen):
    hits = [r for r in rows if r["site_id"] == site and r["date"] == date(year, 11, 1) and r["pathogen"] == pathogen]
    return hits[0] if hits else None


@pytest.fixture(scope="module")
def rows():
    return m.rows_from_pages([page("mddnr2024_table3_p62.txt"), page("mddnr2024_table3_p68.txt")],
                             [page("mddnr2024_table4_p69.txt")])


def test_hand_checked_dermo_values(rows) -> None:
    # Table 3 p. 62: Swan Point 1990 = 7 % / 0.1; Royston 1991 = 100 % / 4.5
    r = lookup(rows, "mddnr_swan_point", 1990, "P_marinus")
    assert r["prevalence"] == pytest.approx(0.07) and r["intensity_mean"] == pytest.approx(0.1)
    r = lookup(rows, "mddnr_royston", 1991, "P_marinus")
    assert r["prevalence"] == pytest.approx(1.0) and r["intensity_mean"] == pytest.approx(4.5)
    # Table 3 p. 68: Back Cove 2024 = 83 % / 4
    r = lookup(rows, "mddnr_back_cove", 2024, "P_marinus")
    assert r["prevalence"] == pytest.approx(0.83) and r["intensity_mean"] == pytest.approx(4.0)


def test_two_ragged_points_kept_apart(rows) -> None:
    # p. 68, 2024: Ragged Point (Little Choptank) 70 % / 2.3; Ragged Point (Potomac) 3 % / 0.03
    assert lookup(rows, "mddnr_ragged_point_little_choptank", 2024, "P_marinus")["prevalence"] == pytest.approx(0.70)
    assert lookup(rows, "mddnr_ragged_point_potomac", 2024, "P_marinus")["intensity_mean"] == pytest.approx(0.03)


def test_hand_checked_msx_values(rows) -> None:
    # Table 4 p. 69: Hackett Point 1992 = 3 %; Cook Point 1992 = 73 %
    assert lookup(rows, "mddnr_hackett_point", 1992, "H_nelsoni")["prevalence"] == pytest.approx(0.03)
    assert lookup(rows, "mddnr_cook_point", 1992, "H_nelsoni")["prevalence"] == pytest.approx(0.73)


def test_na_and_nd_not_written(rows) -> None:
    assert lookup(rows, "mddnr_royston", 1990, "P_marinus") is None      # NA on p. 62
    assert lookup(rows, "mddnr_swan_point", 1994, "H_nelsoni") is None   # ND on p. 69
    assert lookup(rows, "mddnr_cook_point", 1995, "H_nelsoni") is None   # NA on p. 69


def test_every_bar_parsed_and_contract_valid() -> None:
    table = m.disease_obs([page("mddnr2024_table3_p62.txt")], [page("mddnr2024_table4_p69.txt")])
    sites = set(table.column("site_id").to_pylist())
    assert len(sites) == 43
    assert all(n is None for n in table.column("n").to_pylist())


def test_row_order_check_catches_mismatch() -> None:
    text = page("mddnr2024_table3_p62.txt").replace("Hackett Point", "Somewhere Else")
    with pytest.raises(ValueError, match="does not end with bar"):
        m.parse_page(text, 2)


def test_fetch_rejects_non_pdf(tmp_path: Path) -> None:
    with pytest.raises(RuntimeError):
        m.fetch(raw_dir=tmp_path, get=lambda url: b"<html>")
