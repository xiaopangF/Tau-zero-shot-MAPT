"""VUS and unlabeled variant prioritization utilities."""

from __future__ import annotations


def _percentile_rank_desc(scores: list[float], score: float) -> float:
    if not scores:
        return float("nan")
    better_or_equal = sum(value >= score for value in scores)
    return better_or_equal / len(scores)


def risk_tier(percentile_rank: float) -> str:
    if percentile_rank <= 0.01:
        return "top_1pct"
    if percentile_rank <= 0.05:
        return "top_5pct"
    if percentile_rank <= 0.10:
        return "top_10pct"
    return "lower_priority"


def mechanism_summary(row: dict[str, str]) -> str:
    tags = []
    for field in (
        "tau_region",
        "tau_motif",
        "charge_change",
        "special_residue_change",
    ):
        value = row.get(field, "")
        if value and value != "none":
            tags.append(value)
    if row.get("near_ptm_site_3aa") == "True":
        tags.append("PTM_proximal_3aa")
    if row.get("near_known_pathogenic_hotspot_3aa") == "True":
        tags.append("pathogenic_hotspot_proximal_3aa")
    return ";".join(tags)


def prioritize_rows(
    score_rows: list[dict[str, str]],
    annotation_rows: list[dict[str, str]],
    score_column: str,
    target_labels: set[str],
    include_unlabeled: bool = False,
) -> list[dict[str, str]]:
    annotations = {row["variant_id"]: row for row in annotation_rows}
    valid_scores = [float(row[score_column]) for row in score_rows if row.get(score_column, "") != ""]
    rows = []
    for score_row in score_rows:
        if score_row.get(score_column, "") == "":
            continue
        variant_id = score_row["variant_id"]
        annotated = dict(annotations.get(variant_id, {}))
        label = annotated.get("clinvar_label", "")
        if label not in target_labels and not (include_unlabeled and label == ""):
            continue
        score = float(score_row[score_column])
        percentile = _percentile_rank_desc(valid_scores, score)
        merged = dict(score_row)
        merged.update(annotated)
        merged["priority_percentile_rank"] = f"{percentile:.6g}"
        merged["priority_tier"] = risk_tier(percentile)
        merged["mechanism_summary"] = mechanism_summary(merged)
        rows.append(merged)
    return sorted(rows, key=lambda row: float(row[score_column]), reverse=True)
