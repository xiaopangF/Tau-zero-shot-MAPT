"""Publication-oriented plotting helpers."""

from __future__ import annotations

from pathlib import Path

from .constants import AMINO_ACIDS


def create_basic_figures(
    score_rows: list[dict[str, str]],
    label_rows: list[dict[str, str]],
    score_column: str,
    outdir: str | Path,
) -> list[Path]:
    try:
        import matplotlib.pyplot as plt  # type: ignore
        import numpy as np  # type: ignore
    except ImportError as exc:
        raise RuntimeError('Figure generation requires matplotlib and numpy. Install .[analysis].') from exc

    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    label_by_variant = {row["variant_id"]: row.get("clinvar_label", "") for row in label_rows}

    values_by_label: dict[str, list[float]] = {"P_LP": [], "B_LB": [], "VUS": [], "conflicting": []}
    positions: list[int] = []
    scores: list[float] = []
    matrix = np.full((len(AMINO_ACIDS), 441), np.nan)
    aa_to_index = {aa: idx for idx, aa in enumerate(AMINO_ACIDS)}

    for row in score_rows:
        value = row.get(score_column, "")
        if value == "":
            continue
        score = float(value)
        position = int(row["position"])
        mut_aa = row.get("mut_aa", "")
        positions.append(position)
        scores.append(score)
        if mut_aa in aa_to_index and 1 <= position <= 441:
            matrix[aa_to_index[mut_aa], position - 1] = score
        label = label_by_variant.get(row["variant_id"], "")
        if label in values_by_label:
            values_by_label[label].append(score)

    paths: list[Path] = []

    fig, ax = plt.subplots(figsize=(12, 4))
    im = ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap="coolwarm")
    ax.set_xlabel("MAPT/Tau-F position")
    ax.set_ylabel("Mutant amino acid")
    ax.set_yticks(range(len(AMINO_ACIDS)))
    ax.set_yticklabels(AMINO_ACIDS)
    ax.set_title("MAPT missense zero-shot score heatmap")
    fig.colorbar(im, ax=ax, label=score_column)
    fig.tight_layout()
    path = out / "missense_heatmap.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    paths.append(path)

    fig, ax = plt.subplots(figsize=(10, 3))
    ax.scatter(positions, scores, s=5, alpha=0.45, linewidths=0)
    ax.set_xlabel("MAPT/Tau-F position")
    ax.set_ylabel(score_column)
    ax.set_title("MAPT missense zero-shot scores by position")
    fig.tight_layout()
    path = out / "score_by_position.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    paths.append(path)

    nonempty = [(label, values) for label, values in values_by_label.items() if values]
    if nonempty:
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.boxplot([values for _, values in nonempty], tick_labels=[label for label, _ in nonempty])
        ax.set_ylabel(score_column)
        ax.set_title("Zero-shot score distribution by ClinVar label")
        fig.tight_layout()
        path = out / "clinvar_score_distribution.png"
        fig.savefig(path, dpi=300)
        plt.close(fig)
        paths.append(path)

    return paths
