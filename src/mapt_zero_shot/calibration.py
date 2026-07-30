"""Calibration utilities for established MAPT missense controls."""

from __future__ import annotations

from pathlib import Path
from statistics import median


DEFAULT_POSITIVE_CONTROLS = (
    {
        "variant_id": "G272V",
        "evidence_class": "established_pathogenic",
        "evidence_note": "MAPT FTDP-17 missense mutation",
    },
    {
        "variant_id": "P301L",
        "evidence_class": "established_pathogenic",
        "evidence_note": "MAPT FTDP-17 missense mutation",
    },
    {
        "variant_id": "V337M",
        "evidence_class": "established_pathogenic",
        "evidence_note": "MAPT FTDP-17 missense mutation",
    },
    {
        "variant_id": "R406W",
        "evidence_class": "established_pathogenic",
        "evidence_note": "MAPT FTDP-17 missense mutation",
    },
    {
        "variant_id": "N279K",
        "evidence_class": "established_pathogenic",
        "evidence_note": "MAPT FTDP-17 missense mutation",
    },
)


def _rank_and_fraction(scores: list[float], score: float) -> tuple[int, float]:
    rank = 1 + sum(value > score for value in scores)
    return rank, rank / len(scores)


def positive_control_rows(
    score_rows: list[dict[str, str]],
    annotation_rows: list[dict[str, str]] | None = None,
    controls: tuple[dict[str, str], ...] = DEFAULT_POSITIVE_CONTROLS,
    score_column: str = "pathogenic_score_mean",
) -> list[dict[str, object]]:
    """Return ranked ESM calibration rows for established MAPT controls."""
    score_by_id = {
        row["variant_id"]: row
        for row in score_rows
        if row.get("variant_id") and row.get(score_column, "") != ""
    }
    annotations = {row["variant_id"]: row for row in annotation_rows or []}
    valid_scores = sorted(
        (float(row[score_column]) for row in score_by_id.values()),
        reverse=True,
    )
    rows: list[dict[str, object]] = []

    for control in controls:
        variant_id = control["variant_id"]
        score_row = score_by_id.get(variant_id)
        annotation = annotations.get(variant_id, {})
        if score_row is None:
            rows.append(
                {
                    **control,
                    "found_in_atlas": "False",
                    "variant_id": variant_id,
                    "score": "",
                    "score_std": "",
                    "atlas_size": len(valid_scores),
                    "rank": "",
                    "top_fraction": "",
                    "top_percent": "",
                    "top_1pct": "False",
                    "top_5pct": "False",
                    "top_10pct": "False",
                    "calibration_interpretation": "not_found",
                    "position": annotation.get("position", ""),
                    "tau_region": annotation.get("tau_region", ""),
                }
            )
            continue

        score = float(score_row[score_column])
        rank, top_fraction = _rank_and_fraction(valid_scores, score)
        if top_fraction <= 0.01:
            interpretation = "top_1pct"
        elif top_fraction <= 0.05:
            interpretation = "top_5pct"
        elif top_fraction <= 0.10:
            interpretation = "top_10pct"
        else:
            interpretation = "outside_top_10pct"

        rows.append(
            {
                **control,
                "found_in_atlas": "True",
                "variant_id": variant_id,
                "protein_change": score_row.get("protein_change", ""),
                "position": score_row.get("position", annotation.get("position", "")),
                "tau_region": annotation.get("tau_region", ""),
                "score": score,
                "score_std": score_row.get("pathogenic_score_std", ""),
                "atlas_size": len(valid_scores),
                "rank": rank,
                "top_fraction": top_fraction,
                "top_percent": top_fraction * 100,
                "top_1pct": str(top_fraction <= 0.01),
                "top_5pct": str(top_fraction <= 0.05),
                "top_10pct": str(top_fraction <= 0.10),
                "calibration_interpretation": interpretation,
            }
        )
    return rows


def positive_control_summary_rows(
    rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    """Summarize calibration without treating it as clinical validation."""
    found = [row for row in rows if row.get("found_in_atlas") == "True"]
    top_fractions = [float(row["top_fraction"]) for row in found if row.get("top_fraction") != ""]
    return [
        {"metric": "controls_requested", "value": len(rows)},
        {"metric": "controls_found", "value": len(found)},
        {"metric": "atlas_size", "value": found[0]["atlas_size"] if found else 0},
        {"metric": "top_1pct_count", "value": sum(row["top_1pct"] == "True" for row in found)},
        {"metric": "top_5pct_count", "value": sum(row["top_5pct"] == "True" for row in found)},
        {"metric": "top_10pct_count", "value": sum(row["top_10pct"] == "True" for row in found)},
        {
            "metric": "median_top_fraction",
            "value": median(top_fractions) if top_fractions else "",
        },
        {
            "metric": "interpretation",
            "value": (
                "Known controls are not concentrated in the ESM-1v top decile; "
                "use scores for prioritization, not standalone diagnosis."
            ),
        },
    ]


def create_positive_control_plot(rows: list[dict[str, object]], out: Path) -> Path:
    """Create a percentile plot with explicit top-1/5/10% calibration cutoffs."""
    import matplotlib.pyplot as plt  # type: ignore

    plotted = [row for row in rows if row.get("top_percent") not in ("", None)]
    plotted = sorted(plotted, key=lambda row: float(row["top_percent"]), reverse=True)
    labels = [str(row["variant_id"]) for row in plotted]
    values = [float(row["top_percent"]) for row in plotted]

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.barh(labels, values, color="#80b1d3", edgecolor="#264653", linewidth=0.6)
    for cutoff, label in ((1, "top 1%"), (5, "top 5%"), (10, "top 10%")):
        ax.axvline(cutoff, color="#d95f02", linestyle="--", linewidth=1.0)
        ax.text(cutoff + 0.6, len(labels) - 0.25, label, fontsize=8, color="#8c2d04")
    for index, row in enumerate(plotted):
        ax.text(
            float(row["top_percent"]) + 0.8,
            index,
            f"rank {row['rank']}/{row['atlas_size']}",
            va="center",
            fontsize=8,
        )
    ax.set_xlabel("Fraction of atlas scoring at least as high (%)")
    ax.set_title("ESM-1v calibration against established MAPT controls")
    ax.invert_yaxis()
    ax.set_xlim(0, max(values, default=10) * 1.25)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="#dddddd", linewidth=0.6, alpha=0.8)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300)
    plt.close(fig)
    return out
