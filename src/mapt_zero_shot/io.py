"""Small TSV helpers used by the CLI."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path
from typing import Iterable


def open_text(path: str | Path, mode: str = "rt"):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, mode, newline="", encoding="utf-8", errors="replace")
    return path.open(mode, newline="", encoding="utf-8", errors="replace")


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [dict(row) for row in reader]


def write_tsv(path: str | Path, rows: Iterable[dict[str, object]], fieldnames: list[str]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})


def merge_fieldnames(*row_sets: Iterable[dict[str, object]]) -> list[str]:
    names: list[str] = []
    seen: set[str] = set()
    for rows in row_sets:
        for row in rows:
            for name in row:
                if name not in seen:
                    seen.add(name)
                    names.append(name)
    return names
