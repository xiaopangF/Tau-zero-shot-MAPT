"""ClinVar import utilities for MAPT missense benchmarking."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

from .constants import MAPT_441_SEQUENCE
from .io import open_text
from .variants import MissenseVariant, parse_protein_change

PATHOGENIC_TERMS = ("pathogenic", "likely pathogenic")
BENIGN_TERMS = ("benign", "likely benign")


@dataclass(frozen=True)
class ClinVarImportResult:
    accepted_by_variant: dict[str, dict[str, str]]
    rejected_rows: list[dict[str, str]]


def normalize_clinvar_label(clinical_significance: str) -> str:
    text = clinical_significance.lower()
    if "conflicting" in text:
        return "conflicting"
    if "uncertain" in text or "vus" in text:
        return "VUS"
    has_pathogenic = any(term in text for term in PATHOGENIC_TERMS)
    has_benign = any(term in text for term in BENIGN_TERMS)
    if has_pathogenic and not has_benign:
        return "P_LP"
    if has_benign and not has_pathogenic:
        return "B_LB"
    return "other"


def _candidate_text_fields(row: dict[str, str]) -> list[str]:
    fields = []
    for name in ("Name", "ProteinChange", "HGVS", "OtherIDs"):
        value = row.get(name, "")
        if value:
            fields.append(value)
    return fields


def _find_protein_change(row: dict[str, str]) -> MissenseVariant | None:
    for value in _candidate_text_fields(row):
        parsed = parse_protein_change(value)
        if parsed:
            return parsed
    return None


def _reject_row(row: dict[str, str], reason: str, parsed: MissenseVariant | None = None) -> dict[str, str]:
    rejected = {
        "reject_reason": reason,
        "clinvar_variation_id": row.get("VariationID", ""),
        "clinvar_name": row.get("Name", ""),
        "clinvar_significance": row.get("ClinicalSignificance", ""),
        "clinvar_review_status": row.get("ReviewStatus", ""),
        "clinvar_label": normalize_clinvar_label(row.get("ClinicalSignificance", "")),
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


def load_mapt_clinvar_with_qc(
    path: str,
    reference_sequence: str = MAPT_441_SEQUENCE,
    require_reference_match: bool = True,
) -> ClinVarImportResult:
    """Load MAPT ClinVar missense rows and retain rejected-row QC details."""
    accepted: dict[str, dict[str, str]] = {}
    rejected_rows: list[dict[str, str]] = []
    evidence_by_variant: dict[str, list[str]] = defaultdict(list)

    with open_text(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            values = line.rstrip("\n").split("\t")
            row = dict(zip(header, values, strict=False))
            if row.get("GeneSymbol") != "MAPT":
                continue
            parsed = _find_protein_change(row)
            if not parsed:
                rejected_rows.append(_reject_row(row, "no_parseable_missense"))
                continue
            if not (1 <= parsed.position <= len(reference_sequence)):
                rejected_rows.append(_reject_row(row, "outside_reference_range", parsed))
                continue
            reference_wt = reference_sequence[parsed.position - 1]
            if require_reference_match and reference_wt != parsed.wt_aa:
                rejected = _reject_row(row, "reference_wt_mismatch", parsed)
                rejected["reference_wt_aa"] = reference_wt
                rejected_rows.append(rejected)
                continue

            variant_id = parsed.variant_id
            label = normalize_clinvar_label(row.get("ClinicalSignificance", ""))
            evidence_by_variant[variant_id].append(label)
            accepted[variant_id] = {
                "clinvar_variation_id": row.get("VariationID", ""),
                "clinvar_name": row.get("Name", ""),
                "clinvar_significance": row.get("ClinicalSignificance", ""),
                "clinvar_review_status": row.get("ReviewStatus", ""),
                "clinvar_label": label,
                "clinvar_coordinate_qc": "reference_wt_match",
            }

    for variant_id, labels in evidence_by_variant.items():
        unique = sorted(set(labels))
        if len(unique) > 1:
            accepted[variant_id]["clinvar_label"] = "conflicting"
            accepted[variant_id]["clinvar_label_components"] = ";".join(unique)
        else:
            accepted[variant_id]["clinvar_label_components"] = unique[0]

    return ClinVarImportResult(accepted_by_variant=accepted, rejected_rows=rejected_rows)


def load_mapt_clinvar(path: str) -> dict[str, dict[str, str]]:
    """Load accepted MAPT ClinVar rows keyed by 441-aa Tau-F variant_id."""
    return load_mapt_clinvar_with_qc(path).accepted_by_variant
