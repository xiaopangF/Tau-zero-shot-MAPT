"""Label-free fusion of ESM and physicochemical MAPT scores."""

from __future__ import annotations

from math import ceil
from statistics import mean, pstdev


DEFAULT_ESM_WEIGHT = 0.75
DEFAULT_PHYSCHEM_WEIGHT = 0.25


def _zscore_map(rows: list[dict[str, str]], column: str) -> dict[str, float]:
    values = [float(row[column]) for row in rows if row.get(column, "") != ""]
    if not values:
        raise ValueError(f"No numeric values found in score column: {column}")
    center = mean(values)
    scale = pstdev(values) or 1.0
    return {
        str(row["variant_id"]): (float(row[column]) - center) / scale
        for row in rows
        if row.get(column, "") != ""
    }


def composite_zero_shot_rows(
    esm_rows: list[dict[str, str]],
    physchem_rows: list[dict[str, str]],
    esm_column: str = "pathogenic_score_mean",
    physchem_column: str = "physchem_zero_shot_score",
    esm_weight: float = DEFAULT_ESM_WEIGHT,
    physchem_weight: float = DEFAULT_PHYSCHEM_WEIGHT,
) -> list[dict[str, object]]:
    """Combine two zero-shot scores after unlabeled distribution normalization."""

    if esm_weight < 0 or physchem_weight < 0 or esm_weight + physchem_weight <= 0:
        raise ValueError("Fusion weights must be non-negative and have a positive sum")
    total_weight = esm_weight + physchem_weight
    esm_weight /= total_weight
    physchem_weight /= total_weight

    esm_z = _zscore_map(esm_rows, esm_column)
    physchem_z = _zscore_map(physchem_rows, physchem_column)
    common_ids = sorted(set(esm_z) & set(physchem_z))
    if not common_ids:
        raise ValueError("ESM and physicochemical tables have no shared variant IDs")

    esm_by_id = {str(row["variant_id"]): row for row in esm_rows}
    physchem_by_id = {str(row["variant_id"]): row for row in physchem_rows}
    rows: list[dict[str, object]] = []
    for variant_id in common_ids:
        esm_row = esm_by_id[variant_id]
        physchem_row = physchem_by_id[variant_id]
        esm_score = esm_z[variant_id]
        physchem_score = physchem_z[variant_id]
        row = {
            "variant_id": variant_id,
            "protein_change": esm_row.get("protein_change", physchem_row.get("protein_change", "")),
            "position": esm_row.get("position", physchem_row.get("position", "")),
            "wt_aa": esm_row.get("wt_aa", physchem_row.get("wt_aa", "")),
            "mut_aa": esm_row.get("mut_aa", physchem_row.get("mut_aa", "")),
            "esm_raw_score": float(esm_row[esm_column]),
            "physchem_raw_score": float(physchem_row[physchem_column]),
            "esm_standardized_score": esm_score,
            "physchem_standardized_score": physchem_score,
            "esm_weight": esm_weight,
            "physchem_weight": physchem_weight,
            "zero_shot_composite_score": esm_weight * esm_score + physchem_weight * physchem_score,
        }
        for name in ("tau_region", "tau_motif", "near_ptm_site_3aa"):
            if name in physchem_row:
                row[name] = physchem_row[name]
        rows.append(row)

    ranked = sorted(
        rows,
        key=lambda row: (-float(row["zero_shot_composite_score"]), str(row["variant_id"])),
    )
    scores = [float(row["zero_shot_composite_score"]) for row in ranked]
    for row in ranked:
        score = float(row["zero_shot_composite_score"])
        rank = 1 + sum(other > score for other in scores)
        row["zero_shot_composite_rank"] = rank
        row["zero_shot_composite_top_fraction"] = rank / len(ranked)
        row["zero_shot_composite_top_percent"] = 100 * rank / len(ranked)
    return ranked


def composite_summary_rows(
    ranked_rows: list[dict[str, object]],
    controls: tuple[str, ...],
    top_fraction: float = 0.01,
) -> list[dict[str, object]]:
    if not ranked_rows:
        raise ValueError("At least one ranked row is required")
    if not 0 < top_fraction <= 1:
        raise ValueError("top_fraction must be in (0, 1]")
    by_id = {str(row["variant_id"]): row for row in ranked_rows}
    missing = sorted(set(controls) - set(by_id))
    if missing:
        raise ValueError(f"Controls not found in fused table: {missing}")
    ranks = {variant_id: int(by_id[variant_id]["zero_shot_composite_rank"]) for variant_id in controls}
    cutoff = ceil(len(ranked_rows) * top_fraction)
    in_top = sum(rank <= cutoff for rank in ranks.values())
    return [
        {"metric": "model_type", "value": "zero_shot_esm_physchem_fusion"},
        {"metric": "uses_gold_standard_labels_for_scoring", "value": "False"},
        {"metric": "n_variants", "value": len(ranked_rows)},
        {"metric": "top_fraction", "value": top_fraction},
        {"metric": "top_rank_cutoff", "value": cutoff},
        {"metric": "gold_standard_count", "value": len(ranks)},
        {"metric": "gold_standard_in_top_fraction", "value": in_top},
        {"metric": "gold_standard_mean_rank", "value": mean(ranks.values())},
        {"metric": "gold_standard_max_rank", "value": max(ranks.values())},
        *[
            {"metric": f"rank_{variant_id}", "value": ranks[variant_id]}
            for variant_id in sorted(ranks)
        ],
    ]