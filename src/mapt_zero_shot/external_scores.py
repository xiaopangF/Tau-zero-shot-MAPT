"""External predictor import utilities."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

from .constants import MAPT_441_SEQUENCE
from .io import open_text
from .variants import MissenseVariant, parse_protein_change


@dataclass(frozen=True)
class AlphaMissenseImportResult:
    accepted_by_variant: dict[str, dict[str, str]]
    rejected_rows: list[dict[str, str]]


def _normalize_header(name: str) -> str:
    return name.lstrip("#").strip()


def _read_header_and_rows(handle):
    for line in handle:
        if not line.strip():
            continue
        raw_header = line.rstrip("\n").split("\t")
        header = [_normalize_header(name) for name in raw_header]
        if "protein_variant" in header and "am_pathogenicity" in header:
            yield header, handle
            return
    raise ValueError(
        "Could not find AlphaMissense TSV header with protein_variant/am_pathogenicity"
    )


def _reject_row(
    row: dict[str, str], reason: str, parsed: MissenseVariant | None = None
) -> dict[str, str]:
    rejected = {
        "reject_reason": reason,
        "uniprot_id": row.get("uniprot_id", ""),
        "transcript_id": row.get("transcript_id", ""),
        "protein_variant": row.get("protein_variant", ""),
        "am_pathogenicity": row.get("am_pathogenicity", ""),
        "am_class": row.get("am_class", ""),
    }
    if parsed:
        rejected.update(
            {
                "parsed_variant_id": parsed.variant_id,
                "parsed_position": str(parsed.position),
                "parsed_wt_aa": parsed.wt_aa,
                "parsed_mut_aa": parsed.mut_aa,
            }
        )
    return rejected


def load_alphamissense_with_qc(
    path: str,
    transcript_id: str | None = None,
    uniprot_id: str | None = None,
    position_offset: int = 0,
    reference_sequence: str = MAPT_441_SEQUENCE,
    require_reference_match: bool = True,
) -> AlphaMissenseImportResult:
    """Load AlphaMissense rows keyed by 441-aa Tau-F variant_id.

    AlphaMissense files begin with comment/license lines. This reader skips those
    until it finds the real TSV header. Rows are accepted only when optional
    transcript/uniprot filters pass and the parsed protein variant is compatible
    with the 441-aa Tau-F coordinate system.
    """
    by_variant: dict[str, dict[str, str]] = {}
    classes: dict[str, set[str]] = defaultdict(set)
    rejected_rows: list[dict[str, str]] = []

    with open_text(path, "rt") as handle:
        header, rows_iter = next(_read_header_and_rows(handle))
        for line in rows_iter:
            if not line.strip():
                continue
            row = dict(zip(header, line.rstrip("\n").split("\t"), strict=False))
            if transcript_id and row.get("transcript_id") != transcript_id:
                continue
            if uniprot_id and row.get("uniprot_id") != uniprot_id:
                continue

            parsed = parse_protein_change(row.get("protein_variant", ""))
            if not parsed:
                rejected_rows.append(_reject_row(row, "no_parseable_missense"))
                continue

            position = parsed.position + position_offset
            if not 1 <= position <= len(reference_sequence):
                rejected_rows.append(_reject_row(row, "outside_reference_range", parsed))
                continue

            reference_wt = reference_sequence[position - 1]
            if require_reference_match and reference_wt != parsed.wt_aa:
                rejected = _reject_row(row, "reference_wt_mismatch", parsed)
                rejected["reference_wt_aa"] = reference_wt
                rejected_rows.append(rejected)
                continue

            variant_id = f"{parsed.wt_aa}{position}{parsed.mut_aa}"
            try:
                score = float(row.get("am_pathogenicity", ""))
            except ValueError:
                rejected_rows.append(_reject_row(row, "invalid_score", parsed))
                continue

            classes[variant_id].add(row.get("am_class", ""))
            previous = by_variant.get(variant_id)
            if previous is None or score > float(previous["alphamissense_score"]):
                by_variant[variant_id] = {
                    "alphamissense_score": str(score),
                    "alphamissense_class": row.get("am_class", ""),
                    "alphamissense_transcript_id": row.get("transcript_id", ""),
                    "alphamissense_uniprot_id": row.get("uniprot_id", ""),
                    "alphamissense_coordinate_qc": "reference_wt_match",
                }

    for variant_id, class_values in classes.items():
        if variant_id in by_variant:
            by_variant[variant_id]["alphamissense_class_all"] = ";".join(sorted(class_values))
    return AlphaMissenseImportResult(accepted_by_variant=by_variant, rejected_rows=rejected_rows)


def load_alphamissense(
    path: str,
    transcript_id: str | None = None,
    uniprot_id: str | None = None,
    position_offset: int = 0,
) -> dict[str, dict[str, str]]:
    """Load accepted AlphaMissense rows keyed by 441-aa Tau-F variant_id."""
    return load_alphamissense_with_qc(
        path,
        transcript_id=transcript_id,
        uniprot_id=uniprot_id,
        position_offset=position_offset,
    ).accepted_by_variant
