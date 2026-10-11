"""Tests for the LOBO source using saved real responses in tests/fixtures/."""

from datetime import date
from pathlib import Path

import pandas as pd
import pytest

from oyster_data.sources import lobo

FIX = Path(__file__).parent / "fixtures"


def test_build_url() -> None:
    url = lobo.build_url(66, date(2018, 1, 1), date(2018, 12, 31))
    assert url.startswith("http://maine.loboviz.com/cgi-data/nph-data.cgi?")
    assert "min_date=20180101" in url and "max_date=20181231" in url
    assert "y=temperature,salinity,chlorophyll" in url and url.endswith("node=66")


@pytest.mark.parametrize("name", ["lobo_node65_20220701_20220703.tsv", "lobo_node66_20180701_20180703.tsv"])
def test_parse_saved_response(name: str) -> None:
    df = lobo.parse((FIX / name).read_text(encoding="utf-8"))
    assert list(df.columns) == ["time_utc", "temp_c", "salinity_psu", "chl_ug_l"]
    assert len(df) > 60
    assert str(df["time_utc"].dt.tz) == "UTC"
    assert df["temp_c"].between(5, 30).all()


def test_parse_converts_eastern_daylight_time() -> None:
    df = lobo.parse((FIX / "lobo_node66_20180701_20180703.tsv").read_text(encoding="utf-8"))
    first_local = (FIX / "lobo_node66_20180701_20180703.tsv").read_text(encoding="utf-8").splitlines()[3][:19]
    expected = pd.Timestamp(first_local, tz="America/New_York").tz_convert("UTC")
    assert df["time_utc"].iloc[0] == expected
    assert expected.hour == (pd.Timestamp(first_local).hour + 4) % 24   # EDT is UTC-4


def test_parse_repeated_dst_hour() -> None:
    text = ("Sensor 0066 - X\nX\ndate [Eastern]\ttemperature [C]\n"
            "2016-11-06 00:00:00\t10.0\n2016-11-06 01:00:00\t10.1\n"
            "2016-11-06 01:00:00\t10.2\n2016-11-06 02:00:00\t10.3\n")
    df = lobo.parse(text)
    assert df["time_utc"].is_unique and df["time_utc"].is_monotonic_increasing
    assert len(df) == 4


def test_parse_rejects_unknown_format() -> None:
    with pytest.raises(ValueError):
        lobo.parse("no header here\n")
    with pytest.raises(ValueError):
        lobo.parse("t\nt\ndate [UTC]\ttemperature [C]\n")


def test_out_of_water_records_dropped() -> None:
    df = pd.DataFrame({"time_utc": pd.date_range("2016-04-19", periods=3, freq="h", tz="UTC"),
                       "temp_c": [17.8, 7.8, 7.5], "salinity_psu": [0.11, 29.3, 29.6],
                       "chl_ug_l": [0.1, 1.5, 2.0]})
    out = lobo.drop_out_of_water(df)
    assert out["salinity_psu"].tolist() == [29.3, 29.6]


def test_fetch_and_clean_with_saved_response(tmp_path: Path) -> None:
    body = (FIX / "lobo_node66_20180701_20180703.tsv").read_bytes()
    calls: list[str] = []

    def fake_get(url: str) -> bytes:
        calls.append(url)
        return body

    path = lobo.fetch_year(66, 2018, raw_dir=tmp_path, get=fake_get)
    assert path == tmp_path / "lobo" / "node0066_2018.tsv"
    with pytest.raises(FileExistsError):
        lobo.fetch_year(66, 2018, raw_dir=tmp_path, get=fake_get)
    daily = lobo.clean_node(66, raw_dir=tmp_path)
    assert set(daily["site_id"]) == {"dam_upper"}
    assert (daily["source"] == "lobo:node0066").all()
    assert len(daily) >= 2
    assert len(calls) == 2
