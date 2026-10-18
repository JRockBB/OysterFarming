"""Tests for oyster_data.sources.mussel_watch (per-oyster Dermo scores, Chesapeake sites)."""

from oyster_data.contract import SPECS, arrow_schema, validate_table
from oyster_data.sources import chesapeake_sites, mussel_watch
import pyarrow as pa


def test_sites_and_disease_rows_validate():
    rows = mussel_watch.site_rows()
    assert len(rows) == 9 and all(r["region"] == "chesapeake" for r in rows)
    validate_table(pa.Table.from_pylist(rows, schema=arrow_schema("sites")), SPECS["sites"])
    t = mussel_watch.disease_obs()
    d = t.to_pylist()
    assert all(r["intensity_scale"] == "mackin_0_5" and 0 <= r["prevalence"] <= 1 for r in d)
    # Hog Point FY1988: 60 oysters in the source; the row's n and prevalence must come from the per-oyster records
    hog = next(r for r in d if r["site_id"] == "mw_cbhg" and r["date"].year == 1988)
    oy = [o for o in mussel_watch.read_oysters() if o["nst_site"] == "CBHG" and o["fiscal_year"] == 1988]
    assert hog["n"] == len(oy) and abs(hog["prevalence"] - sum(o["dermo_value"] > 0 for o in oy) / len(oy)) < 1e-12


def test_forcing_entries_match_registered_stations():
    from oyster_data.sources import cbp_datahub, coops
    ids = {r["site_id"] for r in mussel_watch.site_rows()}
    for s in chesapeake_sites.MW_SITES:
        assert s.site_id in ids and s.cbp_station in cbp_datahub.STATIONS and s.coops_station in coops.STATIONS
        assert s.cbp_secondary is None or s.cbp_secondary in cbp_datahub.STATIONS
