"""External predictor import utilities."""

from __future__ import annotations

from collections import defaultdict

from .io import open_text
from .variants import parse_protein_change


def _normalize_header(name: str) -> str:
    return name.lstrip("#").strip()


def load_alphamissense(
    path: str,
    transcript_id: str | None = None,
    uniprot_id: str | None = None,
    position_offset: int = 0,
) -> dict[str, dict[str, str]]:
    """Load AlphaMissense rows keyed by MAPT variant_id.

    AlphaMissense distributes both canonical and isoform TSV files. Use the
    isoform-level file when matching 441-aa Tau-F coordinates. If the source
    coordinates differ from this atlas, pass position_offset explicitly.
    """
    by_variant: dict[str, dict[str, str]] = {}
    classes: dict[str, set[str]] = defaultdict(set)

    with open_text(path, "rt") as handle:
        header_line = handle.readline().rstrip("\n")
        header = [_normalize_header(name) for name in header_line.split("\t")]
        for line in handle:
            if not line.strip():
                continue
            row = dict(zip(header, line.rstrip("\n").split("\t"), strict=False))
            if transcript_id and row.get("transcript_id") != transcript_id:
                continue
            if uniprot_id and row.get("uniprot_id") != uniprot_id:
                continue
            parsed = parse_protein_change(row.get("protein_variant", ""))
            if not parsed:
                continue
            position = parsed.position + position_offset
            if not 1 <= position <= 441:
                continue
            variant_id = f"{parsed.wt_aa}{position}{parsed.mut_aa}"
            try:
                score = float(row.get("am_pathogenicity", ""))
            except ValueError:
                continue
            classes[variant_id].add(row.get("am_class", ""))
            previous = by_variant.get(variant_id)
            if previous is None or score > float(previous["alphamissense_score"]):
                by_variant[variant_id] = {
                    "alphamissense_score": str(score),
                    "alphamissense_class": row.get("am_class", ""),
                    "alphamissense_transcript_id": row.get("transcript_id", ""),
                    "alphamissense_uniprot_id": row.get("uniprot_id", ""),
                }

    for variant_id, class_values in classes.items():
        if variant_id in by_variant:
            by_variant[variant_id]["alphamissense_class_all"] = ";".join(sorted(class_values))
    return by_variant
