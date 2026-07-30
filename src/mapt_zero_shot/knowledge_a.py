"""Knowledge-driven Scheme A scoring for MAPT missense variants."""

from __future__ import annotations

from math import ceil
from statistics import mean, pstdev

from .physchem import GOLD_STANDARD_CONTROLS, KYTE_DOOLITTLE, NET_CHARGE

R_REGION_INTERVALS = ((244, 274), (275, 305), (306, 336), (337, 368))
SCHEME_A_PTM_SITES = (202, 205, 208, 210, 214, 217, 235, 262, 396, 400, 403, 404)
MICROTUBULE_INTERFACE_INTERVALS = (
    (255, 260),
    (270, 274),
    (280, 285),
    (295, 300),
    (310, 315),
    (325, 330),
    (340, 345),
    (355, 360),
)
SPLICING_CENTER = 285


def _in_r_region(position: int) -> bool:
    return any(start <= position <= end for start, end in R_REGION_INTERVALS)


def _nearest_distance(position: int, sites: tuple[int, ...]) -> int:
    return min(abs(position - site) for site in sites)


def _charge_sign(aa: str) -> int:
    value = NET_CHARGE[aa]
    if value > 0:
        return 1
    if value < 0:
        return -1
    return 0


def _charge_sign_changed(wt_aa: str, mut_aa: str) -> bool:
    return _charge_sign(wt_aa) != _charge_sign(mut_aa)


def _hydrophobicity_delta(row: dict[str, object]) -> float:
    if row.get("hydrophobicity_delta", "") != "":
        return float(row["hydrophobicity_delta"])
    return KYTE_DOOLITTLE[str(row["mut_aa"])] - KYTE_DOOLITTLE[str(row["wt_aa"])]


def _beta_delta(row: dict[str, object]) -> float:
    if row.get("beta_sheet_delta", "") == "":
        raise ValueError("Scheme A requires beta_sheet_delta in the physicochemical table")
    return float(row["beta_sheet_delta"])


def _net_charge_delta(row: dict[str, object]) -> float:
    if row.get("net_charge_delta", "") != "":
        return float(row["net_charge_delta"])
    return NET_CHARGE[str(row["mut_aa"])] - NET_CHARGE[str(row["wt_aa"])]


def directed_beta_disruption_score(row: dict[str, object]) -> float:
    """Feature 1: R-region beta-gain score with Proline-loss bonus."""

    position = int(row["position"])
    beta_delta = _beta_delta(row)
    if not _in_r_region(position) or beta_delta <= 0.2:
        return 0.0
    score = beta_delta * 2.0
    if str(row["wt_aa"]) == "P":
        score += 1.0
    return min(score, 5.0)


def ptm_microenvironment_score(row: dict[str, object]) -> float:
    """Feature 2: PTM-neighborhood charge and hydrophobicity disruption."""

    position = int(row["position"])
    if _nearest_distance(position, SCHEME_A_PTM_SITES) > 3:
        return 0.0
    score = 0.0
    if _charge_sign_changed(str(row["wt_aa"]), str(row["mut_aa"])):
        score += 2.0
    if abs(_hydrophobicity_delta(row)) > 2.0:
        score += 1.5
    return min(score, 4.0)


def splicing_proximity_score(row: dict[str, object]) -> float:
    """Feature 3: distance to the exon-10 splicing-sensitive center residue."""

    distance = abs(int(row["position"]) - SPLICING_CENTER)
    if distance <= 5:
        return 3.0
    if distance <= 10:
        return 2.0
    if distance <= 25:
        return 1.0
    return 0.0


def microtubule_interface_score(row: dict[str, object]) -> float:
    """Feature 5: fixed prior for Tau residues at the microtubule interface."""

    position = int(row["position"])
    if any(start <= position <= end for start, end in MICROTUBULE_INTERFACE_INTERVALS):
        return 1.5
    return 0.0


def _baseline_physchem_z_scores(rows: list[dict[str, object]]) -> dict[str, float]:
    features = {
        "abs_hydrophobicity_delta": [abs(_hydrophobicity_delta(row)) for row in rows],
        "abs_beta_sheet_delta": [abs(_beta_delta(row)) for row in rows],
        "abs_net_charge_delta": [abs(_net_charge_delta(row)) for row in rows],
    }
    usable: dict[str, tuple[float, float, list[float]]] = {}
    for name, values in features.items():
        standard_deviation = pstdev(values)
        if standard_deviation > 0:
            usable[name] = (mean(values), standard_deviation, values)
    if not usable:
        return {str(row["variant_id"]): 0.0 for row in rows}

    scores: dict[str, float] = {}
    for index, row in enumerate(rows):
        z_values = []
        for center, scale, values in usable.values():
            z_values.append((values[index] - center) / scale)
        scores[str(row["variant_id"])] = mean(z_values)
    return scores


def score_scheme_a_rows(
    rows: list[dict[str, object]],
    controls: tuple[str, ...] = GOLD_STANDARD_CONTROLS,
) -> list[dict[str, object]]:
    """Score and rank all rows with the fixed Scheme A formula."""

    if not rows:
        raise ValueError("At least one physicochemical feature row is required")
    baseline_scores = _baseline_physchem_z_scores(rows)
    control_set = set(controls)
    scored_rows: list[dict[str, object]] = []
    for row in rows:
        scored = dict(row)
        variant_id = str(row["variant_id"])
        beta_score = directed_beta_disruption_score(row)
        ptm_score = ptm_microenvironment_score(row)
        splicing_score = splicing_proximity_score(row)
        interface_score = microtubule_interface_score(row)
        baseline_z = baseline_scores[variant_id]
        final_score = (
            beta_score + ptm_score + splicing_score + interface_score + 0.3 * baseline_z
        )
        scored.update(
            {
                "baseline_physchem_z": baseline_z,
                "directed_beta_score": beta_score,
                "ptm_microenvironment_score": ptm_score,
                "splicing_proximity_score": splicing_score,
                "microtubule_interface_score": interface_score,
                "scheme_a_final_score": final_score,
                "gold_standard_control": str(variant_id in control_set),
            }
        )
        scored_rows.append(scored)

    ranked = sorted(
        scored_rows,
        key=lambda row: (-float(row["scheme_a_final_score"]), str(row["variant_id"])),
    )
    scores = [float(row["scheme_a_final_score"]) for row in ranked]
    for row in ranked:
        score = float(row["scheme_a_final_score"])
        rank = 1 + sum(other > score for other in scores)
        row["scheme_a_rank"] = rank
        row["scheme_a_top_fraction"] = rank / len(ranked)
        row["scheme_a_top_percent"] = 100 * rank / len(ranked)
    return ranked


def scheme_a_gold_rows(
    ranked_rows: list[dict[str, object]],
    controls: tuple[str, ...] = GOLD_STANDARD_CONTROLS,
) -> list[dict[str, object]]:
    by_id = {str(row["variant_id"]): row for row in ranked_rows}
    missing = sorted(set(controls) - set(by_id))
    if missing:
        raise ValueError(f"Gold-standard controls not found in Scheme A table: {missing}")
    return [by_id[variant_id] for variant_id in controls]


def scheme_a_summary_rows(
    ranked_rows: list[dict[str, object]],
    controls: tuple[str, ...] = GOLD_STANDARD_CONTROLS,
    top_fraction: float = 0.01,
) -> list[dict[str, object]]:
    if not ranked_rows:
        raise ValueError("At least one ranked row is required")
    if not 0 < top_fraction <= 1:
        raise ValueError("top_fraction must be in (0, 1]")
    gold_rows = scheme_a_gold_rows(ranked_rows, controls)
    ranks = {str(row["variant_id"]): int(row["scheme_a_rank"]) for row in gold_rows}
    cutoff = ceil(len(ranked_rows) * top_fraction)
    in_top = sum(rank <= cutoff for rank in ranks.values())
    rows: list[dict[str, object]] = [
        {"metric": "model_type", "value": "scheme_a_knowledge_driven_with_interface"},
        {"metric": "uses_gold_standard_labels_for_scoring", "value": "False"},
        {"metric": "n_variants", "value": len(ranked_rows)},
        {"metric": "top_fraction", "value": top_fraction},
        {"metric": "top_rank_cutoff", "value": cutoff},
        {"metric": "gold_standard_count", "value": len(gold_rows)},
        {"metric": "gold_standard_in_top_fraction", "value": in_top},
        {"metric": "gold_standard_mean_rank", "value": mean(ranks.values())},
        {"metric": "gold_standard_max_rank", "value": max(ranks.values())},
    ]
    for variant_id in sorted(ranks):
        rows.append({"metric": f"rank_{variant_id}", "value": ranks[variant_id]})
    return rows