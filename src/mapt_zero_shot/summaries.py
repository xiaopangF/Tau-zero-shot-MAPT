"""Summary tables for MAPT/Tau atlas analyses."""

from __future__ import annotations

from collections import Counter, defaultdict
from statistics import mean, median


def count_rows_by_field(rows: list[dict[str, str]], field: str) -> list[dict[str, object]]:
    counts = Counter(row.get(field, "") or "unlabeled" for row in rows)
    return [{"category": key, "count": counts[key]} for key in sorted(counts)]


def clinvar_summary_rows(
    benchmark_rows: list[dict[str, str]],
    rejected_rows: list[dict[str, str]] | None = None,
) -> list[dict[str, object]]:
    annotated = [row for row in benchmark_rows if row.get("clinvar_label")]
    rows: list[dict[str, object]] = [
        {"section": "accepted", "category": "total_annotated_variants", "count": len(annotated)}
    ]
    for item in count_rows_by_field(annotated, "clinvar_label"):
        rows.append({"section": "accepted_label", **item})

    if rejected_rows is not None:
        rows.append(
            {"section": "rejected", "category": "total_rejected_rows", "count": len(rejected_rows)}
        )
        for item in count_rows_by_field(rejected_rows, "reject_reason"):
            rows.append({"section": "rejected_reason", **item})
        for item in count_rows_by_field(rejected_rows, "clinvar_label"):
            rows.append({"section": "rejected_label", **item})
    return rows


def alphamissense_summary_rows(
    alpha_rows: list[dict[str, str]],
    rejected_rows: list[dict[str, str]] | None = None,
    score_column: str = "alphamissense_score",
) -> list[dict[str, object]]:
    scored = [row for row in alpha_rows if row.get(score_column, "") != ""]
    rows: list[dict[str, object]] = [
        {"section": "coverage", "category": "total_variants", "count": len(alpha_rows)},
        {"section": "coverage", "category": "scored_variants", "count": len(scored)},
        {
            "section": "coverage",
            "category": "scored_positions",
            "count": len({row.get("position") for row in scored}),
        },
    ]
    for item in count_rows_by_field(scored, "tau_region"):
        rows.append({"section": "scored_region", **item})
    for item in count_rows_by_field(scored, "alphamissense_class"):
        rows.append({"section": "scored_class", **item})

    if rejected_rows is not None:
        rows.append(
            {"section": "rejected", "category": "total_rejected_rows", "count": len(rejected_rows)}
        )
        for item in count_rows_by_field(rejected_rows, "reject_reason"):
            rows.append({"section": "rejected_reason", **item})
    return rows

def _float_values(rows: list[dict[str, str]], score_column: str) -> list[float]:
    values = []
    for row in rows:
        value = row.get(score_column, "")
        if value != "":
            values.append(float(value))
    return values


def _quantile(sorted_values: list[float], q: float) -> float:
    if not sorted_values:
        return float("nan")
    index = (len(sorted_values) - 1) * q
    lower = int(index)
    upper = min(lower + 1, len(sorted_values) - 1)
    fraction = index - lower
    return sorted_values[lower] * (1 - fraction) + sorted_values[upper] * fraction


def domain_summary_rows(
    score_rows: list[dict[str, str]],
    annotation_rows: list[dict[str, str]],
    score_column: str,
) -> list[dict[str, object]]:
    annotations = {row["variant_id"]: row for row in annotation_rows}
    merged_by_region: dict[str, list[dict[str, str]]] = defaultdict(list)
    global_scores = sorted(_float_values(score_rows, score_column), reverse=True)
    if not global_scores:
        return []
    top_5_cutoff = global_scores[max(0, round(len(global_scores) * 0.05) - 1)]
    top_10_cutoff = global_scores[max(0, round(len(global_scores) * 0.10) - 1)]

    for score_row in score_rows:
        annotation = annotations.get(score_row["variant_id"], {})
        region = annotation.get("tau_region", "unannotated")
        merged = dict(score_row)
        merged.update(annotation)
        merged_by_region[region].append(merged)

    rows: list[dict[str, object]] = []
    for region in sorted(merged_by_region):
        region_rows = merged_by_region[region]
        scores = _float_values(region_rows, score_column)
        if not scores:
            continue
        sorted_scores = sorted(scores)
        top_5 = sum(score >= top_5_cutoff for score in scores)
        top_10 = sum(score >= top_10_cutoff for score in scores)
        rows.append(
            {
                "tau_region": region,
                "n_variants": len(scores),
                "mean_score": mean(scores),
                "median_score": median(scores),
                "q25_score": _quantile(sorted_scores, 0.25),
                "q75_score": _quantile(sorted_scores, 0.75),
                "max_score": max(scores),
                "top_5pct_global_count": top_5,
                "top_5pct_global_fraction": top_5 / len(scores),
                "top_10pct_global_count": top_10,
                "top_10pct_global_fraction": top_10 / len(scores),
            }
        )
    return rows
