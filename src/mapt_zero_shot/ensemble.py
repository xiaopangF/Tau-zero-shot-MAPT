"""Combine scores from multiple zero-shot models or checkpoints."""

from __future__ import annotations

import math


def _mean(values: list[float]) -> float:
    return sum(values) / len(values)


def _sample_std(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    mean = _mean(values)
    return math.sqrt(sum((value - mean) ** 2 for value in values) / (len(values) - 1))


def ensemble_score_rows(
    score_tables: list[list[dict[str, str]]],
    score_column: str = "pathogenic_score",
) -> list[dict[str, object]]:
    """Average a score column across model/checkpoint TSV tables."""
    by_variant: dict[str, dict[str, object]] = {}
    values_by_variant: dict[str, list[float]] = {}
    models_by_variant: dict[str, list[str]] = {}

    for table in score_tables:
        for row in table:
            variant_id = row["variant_id"]
            value = row.get(score_column, "")
            if value == "":
                continue
            values_by_variant.setdefault(variant_id, []).append(float(value))
            models_by_variant.setdefault(variant_id, []).append(row.get("model", "unknown"))
            if variant_id not in by_variant:
                by_variant[variant_id] = {
                    "variant_id": variant_id,
                    "protein_change": row.get("protein_change", ""),
                    "position": row.get("position", ""),
                    "wt_aa": row.get("wt_aa", ""),
                    "mut_aa": row.get("mut_aa", ""),
                }

    rows: list[dict[str, object]] = []
    for variant_id in sorted(by_variant, key=lambda item: int(by_variant[item]["position"])):
        values = values_by_variant[variant_id]
        row = dict(by_variant[variant_id])
        row.update(
            {
                "ensemble_score_column": score_column,
                "n_models": len(values),
                "models": ";".join(models_by_variant[variant_id]),
                f"{score_column}_mean": _mean(values),
                f"{score_column}_std": _sample_std(values),
            }
        )
        rows.append(row)
    return rows
