"""Generate manuscript-oriented summary figures and markdown from result tables."""

from __future__ import annotations

from pathlib import Path
from statistics import mean


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
            f"{index} | {row.get('variant_id', '')} | {float(row.get('pathogenic_score_mean', 'nan')):.2f} | "
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
        f"The strict P/LP vs B/LB benchmark remains very small. ESM-1v ensemble AUROC={esm_auroc}, AUPRC={esm_auprc}. The transparent Tau heuristic AUROC={heuristic_auroc}. These values should be treated as smoke checks, not clinical performance claims.",
        "",
        "## Domain Signal",
        "",
        "Mean ESM-1v pathogenic score is highest in the microtubule repeat regions and lowest in the N-terminal projection domain.",
        "",
        "| Region | Mean Score | Top 10% Fraction |",
        "|---|---:|---:|",
    ]
    for row in top_domains:
        lines.append(
            f"| {row['tau_region']} | {float(row['mean_score']):.2f} | {float(row['top_10pct_global_fraction']):.3f} |"
        )
    lines.extend(
        [
            "",
            "## ESM Versus Heuristic Baseline",
            "",
            f"Across {scatter_n} variants, ESM-1v ensemble and the transparent Tau heuristic have Pearson r={scatter_r:.2f}. This comparison helps test whether the language model adds signal beyond explicit Tau-domain rules.",
            "",
            "## Top VUS",
            "",
        ]
    )
    for index, row in enumerate(vus_rows[:5], start=1):
        lines.append(
            f"{index}. {row.get('variant_id', '')}: score={float(row.get('pathogenic_score_mean', 'nan')):.2f}; {row.get('mechanism_summary', '')}"
        )
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out
