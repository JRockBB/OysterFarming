"""Tests for oyster_data.sources.dmr_nssp using a saved real response (2025 layer, growing area WQ)."""

import json
from datetime import date
from pathlib import Path

from oyster_data.sources import dmr_nssp

FIX = Path(__file__).parent / "fixtures" / "dmr_nssp2025_WQ_sample.json"


def test_real_sample_converts_closed_sections_only():
    rows, notes = dmr_nssp.to_rows(json.loads(FIX.read_text(encoding="utf-8")), 2025, 2025)
    assert {r["area_id"] for r in rows} == {"WQ:P1", "WQ:P2", "WQ:P3", "WQ:P4", "WQ:P5", "WQ:P6", "WQ:R1", "WQ:CR1"}
    assert all(r["closure_type"] == "pollution" and r["end_date"] is None for r in rows)
    assert next(r for r in rows if r["area_id"] == "WQ:R1")["start_date"] == date(2015, 4, 29)
    assert any("WQ:CA1" in n and "WWTP" in n for n in notes)


def test_seasonal_rule_expands_per_year():
    # SYNTHETIC feature in the layer's format
    data = {"features": [{"attributes": {"GA": "XX", "GA_SECTION": "CA2", "NSSP": "CA", "CA_TYPE": "Seasonal",
                                         "CA_DETAIL": "Closed 6/1 to 10/31", "D_EFFECTIVE": 1477872000000}}]}
    rows, _ = dmr_nssp.to_rows(data, 2024, 2018)
    assert [(r["start_date"], r["end_date"]) for r in rows] == [
        (date(2016, 10, 31), date(2016, 10, 31)), (date(2017, 6, 1), date(2017, 10, 31)),
        (date(2018, 6, 1), date(2018, 10, 31))]
