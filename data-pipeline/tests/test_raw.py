from pathlib import Path

import pytest

from oyster_data.raw import save_raw


def test_save_raw_refuses_overwrite(tmp_path: Path) -> None:
    dest = tmp_path / "a" / "b.txt"
    save_raw(b"one", dest)
    with pytest.raises(FileExistsError):
        save_raw(b"two", dest)
    assert dest.read_bytes() == b"one"
    save_raw(b"two", dest, overwrite=True)
    assert dest.read_bytes() == b"two"
    assert not dest.with_suffix(".txt.part").exists()
