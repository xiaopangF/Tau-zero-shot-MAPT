"""Tau-specific mechanistic annotations."""

from __future__ import annotations

from .constants import KNOWN_PATHOGENIC_HOTSPOTS, MOTIFS, PTM_SITES, REGIONS

CHARGED = set("DEKRH")
POSITIVE = set("KRH")
NEGATIVE = set("DE")
SPECIAL = set("PGC")


def interval_label(position: int, intervals: tuple[tuple[str, int, int], ...]) -> str:
    for label, start, end in intervals:
        if start <= position <= end:
            return label
    return "unannotated"


def region_for_position(position: int) -> str:
    return interval_label(position, REGIONS)


def motif_for_position(position: int) -> str:
    labels = [label for label, start, end in MOTIFS if start <= position <= end]
    return ";".join(labels)


def nearest_distance(position: int, sites: tuple[int, ...]) -> int:
    return min(abs(position - site) for site in sites)


def annotate_variant_row(row: dict[str, str | int]) -> dict[str, str | int]:
    position = int(row["position"])
    wt_aa = str(row["wt_aa"])
    mut_aa = str(row["mut_aa"])

    region = region_for_position(position)
    motif = motif_for_position(position)
    ptm_distance = nearest_distance(position, PTM_SITES)
    hotspot_distance = nearest_distance(position, KNOWN_PATHOGENIC_HOTSPOTS)

    charge_change = "none"
    if wt_aa in POSITIVE and mut_aa in NEGATIVE:
        charge_change = "positive_to_negative"
    elif wt_aa in NEGATIVE and mut_aa in POSITIVE:
        charge_change = "negative_to_positive"
    elif (wt_aa in CHARGED) != (mut_aa in CHARGED):
        charge_change = "charge_gain_or_loss"

    special_change = []
    for aa in SPECIAL:
        if wt_aa == aa and mut_aa != aa:
            special_change.append(f"{aa}_loss")
        elif wt_aa != aa and mut_aa == aa:
            special_change.append(f"{aa}_gain")

    annotated = dict(row)
    annotated.update(
        {
            "tau_region": region,
            "tau_motif": motif,
            "near_ptm_site_3aa": str(ptm_distance <= 3),
            "nearest_ptm_distance": ptm_distance,
            "near_known_pathogenic_hotspot_3aa": str(hotspot_distance <= 3),
            "nearest_hotspot_distance": hotspot_distance,
            "charge_change": charge_change,
            "special_residue_change": ";".join(special_change),
        }
    )
    return annotated

