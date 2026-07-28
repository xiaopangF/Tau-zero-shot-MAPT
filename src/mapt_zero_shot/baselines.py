"""Simple heuristic baseline scores for MAPT missense variants.

These scores are intentionally transparent non-ML baselines. They are useful for
checking whether zero-shot protein language model scores add value beyond obvious
Tau-domain and amino-acid-change rules.
"""

from __future__ import annotations

REGION_WEIGHTS = {
    "microtubule_repeat_R2_exon10": 3.0,
    "microtubule_repeat_R3": 2.5,
    "microtubule_repeat_R4": 2.5,
    "microtubule_repeat_R1": 2.0,
    "proline_rich_region": 1.0,
    "C_terminal_tail": 0.5,
    "N_terminal_projection": 0.25,
}

CHARGE_WEIGHTS = {
    "positive_to_negative": 1.25,
    "negative_to_positive": 1.25,
    "charge_gain_or_loss": 0.75,
    "none": 0.0,
}

SPECIAL_CHANGE_WEIGHTS = {
    "P_loss": 1.25,
    "P_gain": 0.75,
    "G_loss": 0.75,
    "G_gain": 0.5,
    "C_loss": 0.75,
    "C_gain": 0.75,
}


def heuristic_score_row(row: dict[str, str]) -> dict[str, object]:
    score = 0.0
    components: list[str] = []

    region = row.get("tau_region", "")
    region_score = REGION_WEIGHTS.get(region, 0.0)
    score += region_score
    components.append(f"region={region_score:g}")

    if row.get("tau_motif"):
        score += 2.0
        components.append("aggregation_motif=2")

    if row.get("near_known_pathogenic_hotspot_3aa") == "True":
        score += 2.0
        components.append("hotspot_3aa=2")

    if row.get("near_ptm_site_3aa") == "True":
        score += 1.0
        components.append("ptm_3aa=1")

    charge_change = row.get("charge_change", "none")
    charge_score = CHARGE_WEIGHTS.get(charge_change, 0.0)
    score += charge_score
    if charge_score:
        components.append(f"charge={charge_score:g}")

    special_score = 0.0
    for change in row.get("special_residue_change", "").split(";"):
        if not change:
            continue
        special_score += SPECIAL_CHANGE_WEIGHTS.get(change, 0.0)
    score += special_score
    if special_score:
        components.append(f"special_residue={special_score:g}")

    scored = dict(row)
    scored.update(
        {
            "model": "tau_heuristic_v1",
            "heuristic_score": score,
            "heuristic_components": ";".join(components),
        }
    )
    return scored


def heuristic_score_rows(annotation_rows: list[dict[str, str]]) -> list[dict[str, object]]:
    return [heuristic_score_row(row) for row in annotation_rows]
