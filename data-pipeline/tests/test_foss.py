"""Tests for the NOAA FOSS landings source using a saved real response."""

import json
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import pytest

from oyster_data.sources import foss

FIX = Path(__file__).parent / "fixtures" / "foss_maine_oyster.json"


def payload() -> dict:
    return json.loads(FIX.read_text(encoding="utf-8"))


def test_build_url_query() -> None:
    q = parse_qs(urlparse(foss.build_url()).query)
    assert json.loads(q["q"][0]) == {"state_name": "MAINE", "ts_afs_name": {"$like": "%OYSTER%"}}


def test_economics_rows_match_raw_totals() -> None:
    p = payload()
    table = foss.to_economics(p)
    eastern = [i for i in p["items"] if i["ts_afs_name"] == "OYSTER, EASTERN" and i["dollars"] is not None]
    assert table.num_rows == len({i["year"] for i in eastern})
    assert sum(table.column("value_usd").to_pylist()) == pytest.approx(sum(i["dollars"] for i in eastern))
    for row in table.to_pylist():
        recs = [i for i in eastern if i["year"] == row["year"]]
        assert row["value_usd"] == pytest.approx(sum(i["dollars"] for i in recs))
        assert row["pounds"] == pytest.approx(sum(i["pounds"] for i in recs))
        assert row["source"].endswith(f"records={len(recs)}")
        assert row["region"] == "maine" and row["count"] is None
    years = table.column("year").to_pylist()
    assert years == sorted(years) and len(years) == len(set(years))


def test_flat_oyster_excluded() -> None:
    sources = foss.to_economics(payload()).column("source").to_pylist()
    assert all("OYSTER, EASTERN" in s for s in sources)


def test_fetch_saves_raw(tmp_path: Path) -> None:
    body = FIX.read_bytes()
    path = foss.fetch(raw_dir=tmp_path, get=lambda url: body)
    assert path.read_bytes() == body
    with pytest.raises(FileExistsError):
        foss.fetch(raw_dir=tmp_path, get=lambda url: body)
    assert foss.clean(raw_dir=tmp_path).num_rows > 0


def test_fetch_refuses_paginated(tmp_path: Path) -> None:
    with pytest.raises(RuntimeError):
        foss.fetch(raw_dir=tmp_path, get=lambda url: b'{"items": [], "hasMore": true}')
