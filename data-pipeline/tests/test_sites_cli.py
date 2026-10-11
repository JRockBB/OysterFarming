import shutil
from pathlib import Path

from oyster_data import cli
from oyster_data.contract import read_contract
from oyster_data.sites import sites_table, write_sites

FIX = Path(__file__).parent / "fixtures"


def test_sites_table_validates(tmp_path: Path) -> None:
    out = write_sites(tmp_path / "sites.arrow")
    table = read_contract(out, "sites")
    assert table.num_rows == sites_table().num_rows >= 2


def test_build_forcing_end_to_end(tmp_path: Path) -> None:
    raw = tmp_path / "raw" / "lobo"
    raw.mkdir(parents=True)
    shutil.copy(FIX / "lobo_node66_20180701_20180703.tsv", raw / "node0066_2018.tsv")
    shutil.copy(FIX / "lobo_node65_20220701_20220703.tsv", raw / "node0065_2022.tsv")
    out, report = cli.build_forcing(raw_dir=tmp_path / "raw", out_dir=tmp_path)
    table = read_contract(out, "forcing_daily")
    assert set(table.column("site_id").to_pylist()) == {"dam_lower", "dam_upper"}
    assert "## dam_upper" in report.read_text(encoding="utf-8")


def test_years_spec() -> None:
    assert cli._years("2015-2017") == [2015, 2016, 2017]
    assert cli._years("2016,2018") == [2016, 2018]
