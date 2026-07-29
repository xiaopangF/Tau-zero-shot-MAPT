"""Concordance and discordance analysis across MAPT variant scoring methods."""

from __future__ import annotations


def _score_map(rows: list[dict[str, str]], score_column: str) -> dict[str, float]:
    values: dict[str, float] = {}
    for row in rows:
        value = row.get(score_column, "")
        if value == "":
            continue
        values[row["variant_id"]] = float(value)
    return values


def _top_set(scores: dict[str, float], top_fraction: float) -> set[str]:
    if not scores:
        return set()
    count = max(1, round(len(scores) * top_fraction))
    ordered = sorted(scores, key=lambda variant: scores[variant], reverse=True)
    return set(ordered[:count])


def concordance_rows(
    annotations: list[dict[str, str]],
    esm_rows: list[dict[str, str]],
    heuristic_rows: list[dict[str, str]],
    alphamissense_rows: list[dict[str, str]],
    esm_column: str = "pathogenic_score_mean",
    heuristic_column: str = "heuristic_score",
    alphamissense_column: str = "alphamissense_score",
    top_fraction: float = 0.10,
) -> list[dict[str, object]]:
    """Create one row per variant with model top-fraction agreement labels."""
    esm = _score_map(esm_rows, esm_column)
    heuristic = _score_map(heuristic_rows, heuristic_column)
    alphamissense = _score_map(alphamissense_rows, alphamissense_column)

    esm_top = _top_set(esm, top_fraction)
    heuristic_top = _top_set(heuristic, top_fraction)
    alphamissense_top = _top_set(alphamissense, top_fraction)

    rows: list[dict[str, object]] = []
    for annotation in annotations:
        variant_id = annotation["variant_id"]
        esm_high = variant_id in esm_top
        heuristic_high = variant_id in heuristic_top
        alpha_has_score = variant_id in alphamissense
        alpha_high = variant_id in alphamissense_top
        high_count = sum([esm_high, heuristic_high, alpha_high])

        if esm_high and heuristic_high and alpha_high:
            category = "all_three_high"
        elif esm_high and heuristic_high and not alpha_high:
            category = "esm_heuristic_high"
        elif esm_high and alpha_high and not heuristic_high:
            category = "esm_alphamissense_high"
        elif heuristic_high and alpha_high and not esm_high:
            category = "heuristic_alphamissense_high"
        elif esm_high:
            category = "esm_only_high"
        elif heuristic_high:
            category = "heuristic_only_high"
        elif alpha_high:
            category = "alphamissense_only_high"
        else:
            category = "no_model_high"

        rows.append(
            {
                "variant_id": variant_id,
                "protein_change": annotation.get("protein_change", ""),
                "position": annotation.get("position", ""),
                "wt_aa": annotation.get("wt_aa", ""),
                "mut_aa": annotation.get("mut_aa", ""),
                "tau_region": annotation.get("tau_region", ""),
                "tau_motif": annotation.get("tau_motif", ""),
                "charge_change": annotation.get("charge_change", ""),
                "special_residue_change": annotation.get("special_residue_change", ""),
                "esm_score": esm.get(variant_id, ""),
                "heuristic_score": heuristic.get(variant_id, ""),
                "alphamissense_score": alphamissense.get(variant_id, ""),
                "esm_top_decile": esm_high,
                "heuristic_top_decile": heuristic_high,
                "alphamissense_top_decile": alpha_high,
                "alphamissense_has_score": alpha_has_score,
                "n_models_high": high_count,
                "concordance_category": category,
            }
        )
    return rows


def top_concordance_rows(
    rows: list[dict[str, object]],
    category: str | None = None,
    limit: int = 50,
) -> list[dict[str, object]]:
    selected = [row for row in rows if category is None or row["concordance_category"] == category]
    selected.sort(
        key=lambda row: (
            int(row["n_models_high"]),
            float(row["esm_score"]) if row["esm_score"] != "" else float("-inf"),
        ),
        reverse=True,
    )
    return selected[:limit]


def concordance_summary_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    counts: dict[str, int] = {}
    alpha_coverage = 0
    for row in rows:
        category = str(row["concordance_category"])
        counts[category] = counts.get(category, 0) + 1
        if row["alphamissense_has_score"]:
            alpha_coverage += 1

    summary = [
        {"section": "coverage", "category": "total_variants", "count": len(rows)},
        {"section": "coverage", "category": "alphamissense_scored", "count": alpha_coverage},
    ]
    for category in sorted(counts):
        summary.append(
            {"section": "concordance_category", "category": category, "count": counts[category]}
        )
    return summary