"""ClinVar import utilities for MAPT missense benchmarking."""

from __future__ import annotations

from collections import defaultdict

from .io import open_text
from .variants import MissenseVariant, parse_protein_change

PATHOGENIC_TERMS = ("pathogenic", "likely pathogenic")
BENIGN_TERMS = ("benign", "likely benign")


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


def load_mapt_clinvar(path: str) -> dict[str, dict[str, str]]:
    """Load MAPT rows from NCBI variant_summary.txt(.gz), keyed by variant_id."""
    merged: dict[str, dict[str, str]] = {}
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
                continue
            if not (1 <= parsed.position <= 441):
                continue
            variant_id = parsed.variant_id
            label = normalize_clinvar_label(row.get("ClinicalSignificance", ""))
            evidence_by_variant[variant_id].append(label)
            merged[variant_id] = {
                "clinvar_variation_id": row.get("VariationID", ""),
                "clinvar_name": row.get("Name", ""),
                "clinvar_significance": row.get("ClinicalSignificance", ""),
                "clinvar_review_status": row.get("ReviewStatus", ""),
                "clinvar_label": label,
            }

    for variant_id, labels in evidence_by_variant.items():
        unique = sorted(set(labels))
        if len(unique) > 1:
            merged[variant_id]["clinvar_label"] = "conflicting"
            merged[variant_id]["clinvar_label_components"] = ";".join(unique)
        else:
            merged[variant_id]["clinvar_label_components"] = unique[0]

    return merged

