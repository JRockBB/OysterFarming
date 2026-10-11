"""Saving raw downloads. Raw files are never silently overwritten (docs/AGENTS.md)."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import requests

Getter = Callable[[str], bytes]


def http_get(url: str, timeout_s: float = 120.0) -> bytes:
    """GET ``url`` and return the body; raise for HTTP errors."""
    response = requests.get(url, timeout=timeout_s)
    response.raise_for_status()
    return response.content


def save_raw(content: bytes, dest: Path, overwrite: bool = False) -> Path:
    """Write ``content`` to ``dest``. Refuse to replace an existing file unless ``overwrite``."""
    if dest.exists() and not overwrite:
        raise FileExistsError(f"{dest} exists; pass overwrite=True (CLI: --overwrite) to replace it")
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    tmp.write_bytes(content)
    tmp.replace(dest)
    return dest
