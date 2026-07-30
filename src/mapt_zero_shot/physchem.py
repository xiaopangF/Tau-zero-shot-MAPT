"""Zero-shot and legacy physicochemical scoring for MAPT variants.

The zero-shot model in this module never uses pathogenicity labels to select
weights. It measures how unusual a substitution is in three physicochemical
dimensions and applies a fixed Tau-region prior derived from the protein's
known domain organization.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal
from math import ceil
from statistics import mean, pstdev

from .annotations import region_for_position


KYTE_DOOLITTLE = {
    "A": 1.8,
    "R": -4.5,
    "N": -3.5,
    "D": -3.5,
    "C": 2.5,
    "Q": -3.5,
    "E": -3.5,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "L": 3.8,
    "K": -3.9,
    "M": 1.9,
    "F": 2.8,
    "P": -1.6,
    "S": -0.8,
    "T": -0.7,
    "W": -0.9,
    "Y": -1.3,
    "V": 4.2,
}

CHOU_FASMAN_BETA = {
    "A": 0.83,
    "R": 0.93,
    "N": 0.89,
    "D": 0.54,
    "C": 1.19,
    "Q": 1.10,
    "E": 0.37,
    "G": 0.75,
    "H": 0.87,
    "I": 1.60,
    "L": 1.30,
    "K": 0.74,
    "M": 1.05,
    "F": 1.38,
    "P": 0.55,
    "S": 0.75,
    "T": 1.19,
    "W": 1.37,
    "Y": 1.47,
    "V": 1.70,
}

# Simple side-chain charge approximation near physiological pH.
NET_CHARGE = {
    "D": -1.0,
    "E": -1.0,
    "K": 1.0,
    "R": 1.0,
    "H": 0.1,
    "A": 0.0,
    "C": 0.0,
    "F": 0.0,
    "G": 0.0,
    "I": 0.0,
    "L": 0.0,
    "M": 0.0,
    "N": 0.0,
    "P": 0.0,
    "Q": 0.0,
    "S": 0.0,
    "T": 0.0,
    "V": 0.0,
    "W": 0.0,
    "Y": 0.0,
}

GOLD_STANDARD_CONTROLS = ("G272V", "P301L", "V337M", "R406W", "N279K")

# Fixed biological priors. These are not parameters learned from the five
# controls. The repeat region is the main microtubule-binding and aggregation-
# related part of adult Tau.
REGION_MULTIPLIERS = {
    "N_terminal_projection": 0.85,
    "proline_rich_region": 1.00,
    "microtubule_repeat_R1": 1.35,
    "microtubule_repeat_R2_exon10": 1.45,
    "microtubule_repeat_R3": 1.35,
    "microtubule_repeat_R4": 1.35,
    "C_terminal_tail": 1.00,
    "unannotated": 1.00,
}


@dataclass(frozen=True)
class PhyschemWeights:
    hydrophobicity_delta: float
    beta_sheet_delta: float
    net_charge_delta: float


@dataclass(frozen=True)
class ZeroShotPhyschemConfig:
    """Fixed, label-free configuration for the physicochemical model."""

    hydrophobicity_weight: float = 1.0
    beta_sheet_weight: float = 1.0
    net_charge_weight: float = 1.0


@dataclass(frozen=True)
class GridSearchResult:
    weights: PhyschemWeights
    controls_in_top_fraction: int
    mean_control_rank: float
    max_control_rank: int
    control_ranks: dict[str, int]
    top_rank_cutoff: int
    n_variants: int

    @property
    def all_controls_in_top_fraction(self) -> bool:
        return self.controls_in_top_fraction == len(self.control_ranks)


def _scale_delta(scale: dict[str, float], wt_aa: str, mut_aa: str) -> float:
    try:
        return scale[mut_aa] - scale[wt_aa]
    except KeyError as exc:
        raise ValueError(f"Unknown amino-acid symbol in {wt_aa}->{mut_aa}") from exc


def feature_row(row: dict[str, str]) -> dict[str, object]:
    wt_aa = row["wt_aa"]
    mut_aa = row["mut_aa"]
    return {
        **row,
        "hydrophobicity_delta": _scale_delta(KYTE_DOOLITTLE, wt_aa, mut_aa),
        "beta_sheet_delta": _scale_delta(CHOU_FASMAN_BETA, wt_aa, mut_aa),
        "net_charge_delta": _scale_delta(NET_CHARGE, wt_aa, mut_aa),
    }


def feature_rows(variant_rows: list[dict[str, str]]) -> list[dict[str, object]]:
    return [feature_row(row) for row in variant_rows]


def _feature_statistics(
    rows: list[dict[str, object]],
) -> dict[str, tuple[float, float]]:
    if not rows:
        raise ValueError("At least one variant row is required")
    columns = ("hydrophobicity_delta", "beta_sheet_delta", "net_charge_delta")
    statistics: dict[str, tuple[float, float]] = {}
    for column in columns:
        values = [float(row[column]) for row in rows]
        standard_deviation = pstdev(values)
        statistics[column] = (mean(values), standard_deviation or 1.0)
    return statistics


def _region_multiplier(row: dict[str, object]) -> tuple[str, float]:
    position = int(row["position"])
    region = str(row.get("tau_region") or region_for_position(position))
    return region, REGION_MULTIPLIERS.get(region, REGION_MULTIPLIERS["unannotated"])


def zero_shot_physchem_rows(
    rows: list[dict[str, object]],
    config: ZeroShotPhyschemConfig = ZeroShotPhyschemConfig(),
) -> list[dict[str, object]]:
    """Score physicochemical disruption without using disease labels.

    Each feature is standardized across the complete atlas. The absolute
    standardized change measures disruption regardless of direction, which
    avoids cancelling mutations with different mechanisms.
    """

    statistics = _feature_statistics(rows)
    scored_rows: list[dict[str, object]] = []
    feature_weights = {
        "hydrophobicity_delta": config.hydrophobicity_weight,
        "beta_sheet_delta": config.beta_sheet_weight,
        "net_charge_delta": config.net_charge_weight,
    }
    for row in rows:
        scored = dict(row)
        perturbation_score = 0.0
        for column, weight in feature_weights.items():
            center, scale = statistics[column]
            value = abs((float(row[column]) - center) / scale)
            scored[f"standardized_{column}"] = value
            perturbation_score += weight * value
        region, multiplier = _region_multiplier(row)
        scored["tau_region"] = region
        scored["region_multiplier"] = multiplier
        scored["physchem_perturbation_score"] = perturbation_score
        scored["physchem_zero_shot_score"] = perturbation_score * multiplier
        scored_rows.append(scored)
    return scored_rows


def ranked_zero_shot_physchem_rows(
    rows: list[dict[str, object]],
    controls: tuple[str, ...] = GOLD_STANDARD_CONTROLS,
) -> list[dict[str, object]]:
    """Rank rows using the label-free physicochemical score."""

    if not rows:
        raise ValueError("At least one scored row is required")
    control_set = set(controls)
    ranked = sorted(
        (dict(row) for row in rows),
        key=lambda row: (
            -float(row["physchem_zero_shot_score"]),
            str(row["variant_id"]),
        ),
    )
    scores = [float(row["physchem_zero_shot_score"]) for row in ranked]
    for row in ranked:
        score = float(row["physchem_zero_shot_score"])
        rank = 1 + sum(other > score for other in scores)
        row["physchem_rank"] = rank
        row["physchem_top_fraction"] = rank / len(ranked)
        row["physchem_top_percent"] = 100 * rank / len(ranked)
        row["gold_standard_control"] = str(row["variant_id"] in control_set)
    return ranked


def zero_shot_summary_rows(
    ranked_rows: list[dict[str, object]],
    controls: tuple[str, ...] = GOLD_STANDARD_CONTROLS,
    top_fraction: float = 0.01,
) -> list[dict[str, object]]:
    """Create validation metrics for a fixed zero-shot score."""

    if not ranked_rows:
        raise ValueError("At least one ranked row is required")
    if not 0 < top_fraction <= 1:
        raise ValueError("top_fraction must be in (0, 1]")
    control_set = set(controls)
    control_rows = [row for row in ranked_rows if str(row["variant_id"]) in control_set]
    found = {str(row["variant_id"]) for row in control_rows}
    missing = sorted(control_set - found)
    if missing:
        raise ValueError(f"Gold-standard controls not found in variant list: {missing}")
    ranks = {str(row["variant_id"]): int(row["physchem_rank"]) for row in control_rows}
    top_rank_cutoff = ceil(len(ranked_rows) * top_fraction)
    controls_in_top = sum(rank <= top_rank_cutoff for rank in ranks.values())
    mean_rank = mean(ranks.values())
    return [
        {"metric": "model_type", "value": "zero_shot_physchem_perturbation"},
        {"metric": "uses_gold_standard_labels_for_scoring", "value": "False"},
        {"metric": "n_variants", "value": len(ranked_rows)},
        {"metric": "top_fraction", "value": top_fraction},
        {"metric": "top_rank_cutoff", "value": top_rank_cutoff},
        {"metric": "feature_weight_hydrophobicity", "value": 1.0},
        {"metric": "feature_weight_beta_sheet", "value": 1.0},
        {"metric": "feature_weight_net_charge", "value": 1.0},
        {"metric": "gold_standard_count", "value": len(ranks)},
        {"metric": "gold_standard_in_top_fraction", "value": controls_in_top},
        {"metric": "gold_standard_mean_rank", "value": mean_rank},
        {"metric": "gold_standard_max_rank", "value": max(ranks.values())},
        {
            "metric": "all_gold_standard_in_top_fraction",
            "value": str(controls_in_top == len(ranks)),
        },
        *[
            {"metric": f"rank_{variant_id}", "value": ranks[variant_id]}
            for variant_id in sorted(ranks)
        ],
    ]


def _score(row: dict[str, object], weights: PhyschemWeights) -> float:
    return (
        weights.hydrophobicity_delta * float(row["hydrophobicity_delta"])
        + weights.beta_sheet_delta * float(row["beta_sheet_delta"])
        + weights.net_charge_delta * float(row["net_charge_delta"])
    )


def _grid_values(start: float, stop: float, step: float) -> list[float]:
    if step <= 0:
        raise ValueError("grid step must be positive")
    if start > stop:
        raise ValueError("grid min must be <= grid max")
    values = []
    current = Decimal(str(start))
    end = Decimal(str(stop))
    increment = Decimal(str(step))
    while current <= end:
        values.append(float(current))
        current += increment
    return values


def _unique_feature_counts(
    rows: list[dict[str, object]],
) -> dict[tuple[float, float, float], int]:
    counts: dict[tuple[float, float, float], int] = {}
    for row in rows:
        key = (
            float(row["hydrophobicity_delta"]),
            float(row["beta_sheet_delta"]),
            float(row["net_charge_delta"]),
        )
        counts[key] = counts.get(key, 0) + 1
    return counts


def _control_feature_map(
    rows: list[dict[str, object]],
    controls: tuple[str, ...],
) -> dict[str, tuple[float, float, float]]:
    wanted = set(controls)
    found: dict[str, tuple[float, float, float]] = {}
    for row in rows:
        variant_id = str(row["variant_id"])
        if variant_id in wanted:
            found[variant_id] = (
                float(row["hydrophobicity_delta"]),
                float(row["beta_sheet_delta"]),
                float(row["net_charge_delta"]),
            )
    missing = sorted(wanted - set(found))
    if missing:
        raise ValueError(f"Gold-standard controls not found in variant list: {missing}")
    return found


def _feature_score(features: tuple[float, float, float], weights: PhyschemWeights) -> float:
    return (
        weights.hydrophobicity_delta * features[0]
        + weights.beta_sheet_delta * features[1]
        + weights.net_charge_delta * features[2]
    )


def _control_ranks_from_unique_features(
    feature_counts: dict[tuple[float, float, float], int],
    controls: dict[str, tuple[float, float, float]],
    weights: PhyschemWeights,
) -> dict[str, int]:
    scored_features = [
        (_feature_score(features, weights), count) for features, count in feature_counts.items()
    ]
    ranks = {}
    for variant_id, features in controls.items():
        control_score = _feature_score(features, weights)
        ranks[variant_id] = 1 + sum(
            count for score, count in scored_features if score > control_score
        )
    return ranks


def grid_search_weights(
    rows: list[dict[str, object]],
    controls: tuple[str, ...] = GOLD_STANDARD_CONTROLS,
    grid_min: float = -5.0,
    grid_max: float = 5.0,
    grid_step: float = 0.5,
    top_fraction: float = 0.01,
) -> GridSearchResult:
    """Legacy label-informed search retained for reproducibility only."""

    if not 0 < top_fraction <= 1:
        raise ValueError("top_fraction must be in (0, 1]")
    feature_counts = _unique_feature_counts(rows)
    control_features = _control_feature_map(rows, controls)
    values = _grid_values(grid_min, grid_max, grid_step)
    top_rank_cutoff = ceil(len(rows) * top_fraction)
    best: GridSearchResult | None = None

    for hydrophobicity_weight in values:
        for beta_weight in values:
            for charge_weight in values:
                if hydrophobicity_weight == beta_weight == charge_weight == 0:
                    continue
                weights = PhyschemWeights(
                    hydrophobicity_delta=hydrophobicity_weight,
                    beta_sheet_delta=beta_weight,
                    net_charge_delta=charge_weight,
                )
                ranks = _control_ranks_from_unique_features(
                    feature_counts, control_features, weights
                )
                controls_in_top = sum(rank <= top_rank_cutoff for rank in ranks.values())
                mean_rank = mean(ranks.values())
                max_rank = max(ranks.values())
                candidate = GridSearchResult(
                    weights=weights,
                    controls_in_top_fraction=controls_in_top,
                    mean_control_rank=mean_rank,
                    max_control_rank=max_rank,
                    control_ranks=ranks,
                    top_rank_cutoff=top_rank_cutoff,
                    n_variants=len(rows),
                )
                candidate_key = (
                    candidate.controls_in_top_fraction,
                    -candidate.max_control_rank,
                    -candidate.mean_control_rank,
                )
                if best is None:
                    best = candidate
                    continue
                best_key = (
                    best.controls_in_top_fraction,
                    -best.max_control_rank,
                    -best.mean_control_rank,
                )
                if candidate_key > best_key:
                    best = candidate

    if best is None:
        raise ValueError("Grid search produced no non-zero weight combinations")
    return best


def ranked_physchem_rows(
    rows: list[dict[str, object]],
    weights: PhyschemWeights,
    controls: tuple[str, ...] = GOLD_STANDARD_CONTROLS,
) -> list[dict[str, object]]:
    control_set = set(controls)
    scored_rows = []
    for row in rows:
        scored = dict(row)
        scored["physchem_score"] = _score(row, weights)
        scored["gold_standard_control"] = str(row["variant_id"] in control_set)
        scored_rows.append(scored)
    ranked = sorted(
        scored_rows,
        key=lambda row: (-float(row["physchem_score"]), str(row["variant_id"])),
    )
    scores = [float(row["physchem_score"]) for row in ranked]
    for index, row in enumerate(ranked, start=1):
        score = float(row["physchem_score"])
        rank = 1 + sum(other > score for other in scores)
        row["physchem_rank"] = rank
        row["physchem_top_fraction"] = rank / len(ranked)
        row["physchem_top_percent"] = 100 * rank / len(ranked)
    return ranked


def summary_rows(
    result: GridSearchResult,
    grid_min: float,
    grid_max: float,
    grid_step: float,
    top_fraction: float,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = [
        {"metric": "model_type", "value": "legacy_label_informed_linear_grid"},
        {"metric": "n_variants", "value": result.n_variants},
        {"metric": "top_fraction", "value": top_fraction},
        {"metric": "top_rank_cutoff", "value": result.top_rank_cutoff},
        {"metric": "grid_min", "value": grid_min},
        {"metric": "grid_max", "value": grid_max},
        {"metric": "grid_step", "value": grid_step},
        {"metric": "weight_hydrophobicity_delta", "value": result.weights.hydrophobicity_delta},
        {"metric": "weight_beta_sheet_delta", "value": result.weights.beta_sheet_delta},
        {"metric": "weight_net_charge_delta", "value": result.weights.net_charge_delta},
        {"metric": "gold_standard_count", "value": len(result.control_ranks)},
        {"metric": "gold_standard_in_top_fraction", "value": result.controls_in_top_fraction},
        {"metric": "gold_standard_mean_rank", "value": result.mean_control_rank},
        {"metric": "gold_standard_max_rank", "value": result.max_control_rank},
        {
            "metric": "all_gold_standard_in_top_fraction",
            "value": str(result.all_controls_in_top_fraction),
        },
    ]
    for variant_id in sorted(result.control_ranks):
        rows.append({"metric": f"rank_{variant_id}", "value": result.control_ranks[variant_id]})
    return rows