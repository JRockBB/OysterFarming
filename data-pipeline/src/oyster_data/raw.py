"""Saving raw downloads. Raw files are never silently overwritten (docs/AGENTS.md)."""

from __future__ import annotations

import time
from collections.abc import Callable
from pathlib import Path

import requests

Getter = Callable[[str], bytes]


def http_get(url: str, timeout_s: float = 300.0, retries: int = 4, backoff_s: float = 15.0) -> bytes:
    """GET ``url`` and return the body; raise for HTTP errors.

    Server errors (HTTP 5xx) and connection errors are retried ``retries`` times
    with linearly increasing waits; client errors (4xx) are raised at once.
    """
    for attempt in range(retries + 1):
        try:
            response = requests.get(url, timeout=timeout_s)
            if response.status_code < 500:
                response.raise_for_status()
                return response.content
            if attempt == retries:
                response.raise_for_status()
        except (requests.ConnectionError, requests.Timeout):
            if attempt == retries:
                raise
        time.sleep(backoff_s * (attempt + 1))
    raise RuntimeError('unreachable')


def save_raw(content: bytes, dest: Path, overwrite: bool = False) -> Path:
    """Write ``content`` to ``dest``. Refuse to replace an existing file unless ``overwrite``."""
    if dest.exists() and not overwrite:
        raise FileExistsError(f"{dest} exists; pass overwrite=True (CLI: --overwrite) to replace it")
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    tmp.write_bytes(content)
    tmp.replace(dest)
    return dest
