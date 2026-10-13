"""Tests for the transcribed Maine DMR aquaculture oyster landings."""

import pytest

from oyster_data.sources import dmr_aquaculture as d


def test_transcription_matches_report_totals() -> None:
    d.check_totals()            # raises if any transcribed value is wrong


def test_corrupted_value_detected(monkeypatch) -> None:
    rows = dict(d.OYSTER_ROWS)
    lb, n, v = rows[2023]
    rows[2023] = (lb, n + 100, v)   # a one-digit slip
    monkeypatch.setattr(d, "OYSTER_ROWS", rows)
    with pytest.raises(ValueError):
        d.check_totals()


def test_economics_rows() -> None:
    t = d.to_economics().to_pylist()
    assert [r["year"] for r in t] == [2021, 2022, 2023, 2024, 2025]
    r24 = t[3]
    assert r24["region"] == "maine_aquaculture" and r24["count"] == 23_401_019
    assert r24["value_usd"] / r24["count"] == pytest.approx(0.7785, abs=1e-4)   # farm-gate $/oyster
