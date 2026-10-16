"""Tests for oyster_data.sources.nh_greatbay (transcribed NH Fish and Game Dermo table)."""

from oyster_data.contract import SPECS, validate_table
from oyster_data.sources import nh_greatbay as nh


def test_transcribed_rows_are_consistent():
    nh.check_rows()
    assert len(nh.DERMO_ROWS) == 67
    assert {r[1] for r in nh.DERMO_ROWS} == set(nh.SITES)


def test_disease_obs_values_match_source_rows():
    t = nh.disease_obs().to_pylist()
    # Table 2, 10/17/14 Nannie Island: n = 9, classes 3:3, 4:4, 5:1 -> 8/9 infected, intensity (9+16+5)/9
    r = next(x for x in t if x["site_id"] == "nh_nannie_island" and str(x["date"]) == "2014-10-17")
    assert r["n"] == 9 and abs(r["prevalence"] - 8 / 9) < 1e-12
    assert abs(r["intensity_mean"] - 30 / 9) < 1e-12 and r["intensity_scale"] == "mackin_0_5"


def test_site_rows_validate_with_new_region():
    import pyarrow as pa
    from oyster_data.contract import arrow_schema
    t = pa.Table.from_pylist(nh.site_rows(), schema=arrow_schema("sites"))
    validate_table(t, SPECS["sites"])
    assert set(t.column("region").to_pylist()) == {"great_bay_nh"}


def test_met_file_air_temperature_parses():
    # SYNTHETIC two-day meteorological file in the NERRS archive format (DatetimeStamp spelling, ATemp/F_ATemp)
    from oyster_data.sources import nerrs_ncei
    rows = ["DatetimeStamp,Historical,ProvisionalPlus,Frequency,F_Record,ATemp,F_ATemp"]
    for d, t in (("01/05/2003", -3.0), ("01/06/2003", 1.0)):
        rows += [f"{d} {h:02d}:{m:02d},1,1,15,,{t + 0.1 * h},<0>" for h in range(24) for m in (0, 15, 30, 45)]
    s = nerrs_ncei.parse("\n".join(rows), "ATemp")
    assert s.name == "air_temp_c" and len(s) == 2
    assert abs(s.iloc[0] - (-3.0 + 0.1 * 11.5)) < 1e-9
