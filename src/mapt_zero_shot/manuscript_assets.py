"""Generate manuscript-oriented summary figures and markdown from result tables."""

from __future__ import annotations

from pathlib import Path
from statistics import mean

from .constants import REGIONS


REGION_COLORS = {
    "N_terminal_projection": "#8dd3c7",
    "proline_rich_region": "#ffffb3",
    "microtubule_repeat_R1": "#bebada",
    "microtubule_repeat_R2_exon10": "#fb8072",
    "microtubule_repeat_R3": "#80b1d3",
    "microtubule_repeat_R4": "#fdb462",
    "C_terminal_tail": "#b3de69",
}


def _metric_value(rows: list[dict[str, str]], metric: str, model: str | None = None) -> str:
    for row in rows:
        if row.get("metric") != metric:
            continue
        if model is not None and row.get("model") != model:
            continue
        return row.get("value", "")
    return ""


def _as_float(value: str) -> float:
    return float(value) if value not in ("", "nan") else float("nan")


def _join_by_variant(rows: list[dict[str, str]], value_column: str) -> dict[str, float]:
    values: dict[str, float] = {}
    for row in rows:
        value = row.get(value_column, "")
        if value != "":
            values[row["variant_id"]] = float(value)
    return values


def _summary_count(rows: list[dict[str, str]], section: str, category: str) -> str:
    for row in rows:
        if row.get("section") == section and row.get("category") == category:
            return row.get("count", "")
    return ""


def pearson(x_values: list[float], y_values: list[float]) -> float:
    if len(x_values) != len(y_values) or len(x_values) < 2:
        return float("nan")
    x_mean = mean(x_values)
    y_mean = mean(y_values)
    numerator = sum((x - x_mean) * (y - y_mean) for x, y in zip(x_values, y_values, strict=True))
    x_denom = sum((x - x_mean) ** 2 for x in x_values) ** 0.5
    y_denom = sum((y - y_mean) ** 2 for y in y_values) ** 0.5
    if x_denom == 0 or y_denom == 0:
        return float("nan")
    return numerator / (x_denom * y_denom)


def create_domain_summary_plot(domain_rows: list[dict[str, str]], out: Path) -> Path:
    import matplotlib.pyplot as plt  # type: ignore

    ordered = sorted(domain_rows, key=lambda row: float(row["mean_score"]), reverse=True)
    labels = [row["tau_region"].replace("microtubule_repeat_", "R") for row in ordered]
    means = [float(row["mean_score"]) for row in ordered]
    top10 = [float(row["top_10pct_global_fraction"]) for row in ordered]

    fig, ax1 = plt.subplots(figsize=(10, 4.8))
    bars = ax1.bar(labels, means, color="#386cb0")
    ax1.set_ylabel("Mean ESM-1v pathogenic score")
    ax1.set_xlabel("Tau region")
    ax1.tick_params(axis="x", rotation=35)
    ax1.bar_label(bars, fmt="%.1f", padding=2, fontsize=8)

    ax2 = ax1.twinx()
    ax2.plot(labels, top10, color="#fdb462", marker="o", linewidth=2)
    ax2.set_ylabel("Fraction in global top 10%")
    ax2.set_ylim(0, max(top10) * 1.25 if top10 else 1)

    ax1.set_title("ESM-1v ensemble sensitivity by Tau region")
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300)
    plt.close(fig)
    return out


def create_model_comparison_plot(comparison_rows: list[dict[str, str]], out: Path) -> Path:
    import matplotlib.pyplot as plt  # type: ignore

    models = []
    aurocs = []
    auprcs = []
    for row in comparison_rows:
        model = row.get("model", "")
        if model and model not in models:
            models.append(model)
    for model in models:
        aurocs.append(_as_float(_metric_value(comparison_rows, "AUROC", model)))
        auprcs.append(_as_float(_metric_value(comparison_rows, "AUPRC", model)))

    x = list(range(len(models)))
    width = 0.36
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.bar([i - width / 2 for i in x], aurocs, width, label="AUROC", color="#7fc97f")
    ax.bar([i + width / 2 for i in x], auprcs, width, label="AUPRC", color="#beaed4")
    ax.set_xticks(x)
    ax.set_xticklabels(models, rotation=20, ha="right")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Metric value")
    ax.set_title("Strict ClinVar benchmark smoke comparison")
    ax.legend(frameon=False)
    ax.text(
        0.01,
        0.02,
        "Strict P/LP vs B/LB benchmark has n=3; not a performance claim.",
        transform=ax.transAxes,
        fontsize=8,
    )
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300)
    plt.close(fig)
    return out


def create_score_scatter_plot(
    esm_rows: list[dict[str, str]],
    heuristic_rows: list[dict[str, str]],
    out: Path,
) -> tuple[Path, float, int]:
    import matplotlib.pyplot as plt  # type: ignore

    esm = _join_by_variant(esm_rows, "pathogenic_score_mean")
    heuristic = _join_by_variant(heuristic_rows, "heuristic_score")
    variants = sorted(set(esm) & set(heuristic))
    x = [heuristic[variant] for variant in variants]
    y = [esm[variant] for variant in variants]
    r = pearson(x, y)

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(x, y, s=5, alpha=0.28, linewidths=0, color="#4daf4a")
    ax.set_xlabel("Tau heuristic score")
    ax.set_ylabel("ESM-1v ensemble pathogenic score")
    ax.set_title("ESM-1v vs transparent Tau heuristic")
    ax.text(0.03, 0.96, f"n={len(variants)}\nPearson r={r:.2f}", transform=ax.transAxes, va="top")
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300)
    plt.close(fig)
    return out, r, len(variants)


def create_workflow_schematic(
    out: Path,
    atlas_variant_count: int,
    clinvar_summary_rows: list[dict[str, str]] | None = None,
    alphamissense_summary_rows: list[dict[str, str]] | None = None,
) -> Path:
    import matplotlib.pyplot as plt  # type: ignore
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # type: ignore

    clinvar_summary_rows = clinvar_summary_rows or []
    alphamissense_summary_rows = alphamissense_summary_rows or []
    clinvar_accepted = _summary_count(
        clinvar_summary_rows, "accepted", "total_annotated_variants"
    )
    alpha_scored = _summary_count(alphamissense_summary_rows, "coverage", "scored_variants")

    fig, ax = plt.subplots(figsize=(11.5, 5.4))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    boxes = [
        ((0.05, 0.58), "Tau-F reference", "441-aa 2N4R / hTau40\nP10636-8 coordinates", "#e5f5f9"),
        (
            (0.28, 0.58),
            "Missense atlas",
            f"{atlas_variant_count:,} variants\n19 substitutions/site",
            "#edf8e9",
        ),
        (
            (0.51, 0.58),
            "Zero-shot scoring",
            "ESM-1v five-checkpoint\nmasked marginal scores",
            "#fef0d9",
        ),
        ((0.74, 0.58), "Main outputs", "Domain signal\nVUS prioritization", "#f2f0f7"),
        (
            (0.28, 0.18),
            "ClinVar QC",
            f"strict WT match\naccepted: {clinvar_accepted or 'NA'}",
            "#fee8c8",
        ),
        (
            (0.51, 0.18),
            "AlphaMissense QC",
            f"strict WT match\nscored: {alpha_scored or 'NA'}",
            "#e6f5c9",
        ),
    ]

    for (x, y), title, body, color in boxes:
        patch = FancyBboxPatch(
            (x, y),
            0.18,
            0.22,
            boxstyle="round,pad=0.018,rounding_size=0.018",
            linewidth=1.1,
            edgecolor="#404040",
            facecolor=color,
        )
        ax.add_patch(patch)
        ax.text(x + 0.09, y + 0.155, title, ha="center", va="center", weight="bold", fontsize=10)
        ax.text(x + 0.09, y + 0.075, body, ha="center", va="center", fontsize=8.5)

    arrow_pairs = [
        ((0.23, 0.69), (0.28, 0.69)),
        ((0.46, 0.69), (0.51, 0.69)),
        ((0.69, 0.69), (0.74, 0.69)),
        ((0.37, 0.58), (0.37, 0.40)),
        ((0.60, 0.58), (0.60, 0.40)),
        ((0.46, 0.29), (0.51, 0.29)),
    ]
    for start, end in arrow_pairs:
        ax.add_patch(
            FancyArrowPatch(
                start,
                end,
                arrowstyle="-|>",
                mutation_scale=14,
                linewidth=1.1,
                color="#404040",
            )
        )

    ax.text(
        0.05,
        0.93,
        "MAPT/Tau zero-shot atlas workflow",
        ha="left",
        va="center",
        fontsize=15,
        weight="bold",
    )
    ax.text(
        0.05,
        0.875,
        "Clinical labels and external predictors are added after label-free atlas scoring.",
        ha="left",
        va="center",
        fontsize=9,
        color="#333333",
    )
    ax.text(
        0.74,
        0.22,
        (
            "Interpretation rule:\nreport QC and coverage;\n"
            "do not overclaim tiny\nClinVar AUROC/AUPRC."
        ),
        ha="left",
        va="center",
        fontsize=8.5,
        color="#333333",
    )

    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300)
    plt.close(fig)
    return out


def create_vus_lollipop_plot(
    vus_rows: list[dict[str, str]],
    out: Path,
    limit: int = 15,
    score_column: str = "pathogenic_score_mean",
) -> Path:
    import matplotlib.pyplot as plt  # type: ignore
    from matplotlib.lines import Line2D  # type: ignore

    rows = sorted(vus_rows[:limit], key=lambda row: int(row["position"]))
    fig, ax = plt.subplots(figsize=(11.5, 5.2))

    region_y = -0.9
    for region, start, end in REGIONS:
        ax.axvspan(
            start,
            end,
            ymin=0.02,
            ymax=0.12,
            color=REGION_COLORS.get(region, "#dddddd"),
            alpha=0.85,
        )
        label = region.replace("microtubule_repeat_", "").replace("_", " ")
        ax.text((start + end) / 2, region_y, label, ha="center", va="center", fontsize=7)

    max_score = max(float(row[score_column]) for row in rows) if rows else 1.0
    for index, row in enumerate(rows):
        position = int(row["position"])
        score = float(row[score_column])
        region = row.get("tau_region", "")
        color = REGION_COLORS.get(region, "#666666")
        ax.vlines(position, 0, score, color="#555555", linewidth=1.0, alpha=0.75)
        ax.scatter(position, score, s=54, color=color, edgecolor="#222222", linewidth=0.5, zorder=3)
        label_y = score + 0.42 + (index % 2) * 0.42
        ax.text(position, label_y, row.get("variant_id", ""), ha="center", va="bottom", fontsize=8)

    ax.set_xlim(1, 441)
    ax.set_ylim(-1.25, max_score + 1.7)
    ax.set_xlabel("Tau-F amino-acid position")
    ax.set_ylabel("ESM-1v ensemble pathogenic score")
    ax.set_title(f"Top {len(rows)} ClinVar VUS candidates across Tau-F")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#dddddd", linewidth=0.6, alpha=0.8)

    handles = [
        Line2D(
            [0], [0], marker="o", color="none", markerfacecolor=color, label=region, markersize=7
        )
        for region, color in REGION_COLORS.items()
        if any(row.get("tau_region") == region for row in rows)
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=7, frameon=False, ncols=2)

    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300)
    plt.close(fig)
    return out


def write_top_vus_markdown(vus_rows: list[dict[str, str]], out: Path, limit: int = 10) -> Path:
    out.parent.mkdir(parents=True, exist_ok=True)
    rows = vus_rows[:limit]
    lines = [
        "# Top ClinVar VUS Candidates",
        "",
        "Ranked by ESM-1v ensemble `pathogenic_score_mean` under strict Tau-F coordinate QC.",
        "",
        "| Rank | Variant | Score | Region | Mechanism | ClinVar review |",
        "|---:|---|---:|---|---|---|",
    ]
    for index, row in enumerate(rows, start=1):
        lines.append(
            "| "
            f"{index} | {row.get('variant_id', '')} | "
            f"{float(row.get('pathogenic_score_mean', 'nan')):.2f} | "
            f"{row.get('tau_region', '')} | {row.get('mechanism_summary', '')} | "
            f"{row.get('clinvar_review_status', '')} |"
        )
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


def write_results_summary(
    domain_rows: list[dict[str, str]],
    comparison_rows: list[dict[str, str]],
    vus_rows: list[dict[str, str]],
    scatter_r: float,
    scatter_n: int,
    out: Path,
) -> Path:
    out.parent.mkdir(parents=True, exist_ok=True)
    top_domains = sorted(domain_rows, key=lambda row: float(row["mean_score"]), reverse=True)
    esm_auroc = _metric_value(comparison_rows, "AUROC", "esm1v_ensemble")
    esm_auprc = _metric_value(comparison_rows, "AUPRC", "esm1v_ensemble")
    heuristic_auroc = _metric_value(comparison_rows, "AUROC", "tau_heuristic_v1")

    lines = [
        "# ESM-1v Ensemble Result Summary",
        "",
        "## Strict ClinVar Benchmark",
        "",
        (
            "The strict P/LP vs B/LB benchmark remains very small. "
            f"ESM-1v ensemble AUROC={esm_auroc}, AUPRC={esm_auprc}. "
            f"The transparent Tau heuristic AUROC={heuristic_auroc}. "
            "These values should be treated as smoke checks, not clinical performance claims."
        ),
        "",
        "## Domain Signal",
        "",
        (
            "Mean ESM-1v pathogenic score is highest in the microtubule repeat "
            "regions and lowest in the N-terminal projection domain."
        ),
        "",
        "| Region | Mean Score | Top 10% Fraction |",
        "|---|---:|---:|",
    ]
    for row in top_domains:
        lines.append(
            f"| {row['tau_region']} | {float(row['mean_score']):.2f} | "
            f"{float(row['top_10pct_global_fraction']):.3f} |"
        )
    lines.extend(
        [
            "",
            "## ESM Versus Heuristic Baseline",
            "",
            (
                f"Across {scatter_n} variants, ESM-1v ensemble and the transparent Tau "
                f"heuristic have Pearson r={scatter_r:.2f}. This comparison helps test "
                "whether the language model adds signal beyond explicit Tau-domain rules."
            ),
            "",
            "## Top VUS",
            "",
        ]
    )
    for index, row in enumerate(vus_rows[:5], start=1):
        lines.append(
            f"{index}. {row.get('variant_id', '')}: "
            f"score={float(row.get('pathogenic_score_mean', 'nan')):.2f}; "
            f"{row.get('mechanism_summary', '')}"
        )
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out