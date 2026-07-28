"""Benchmark utilities for label-free MAPT variant scores."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class BinaryExample:
    variant_id: str
    score: float
    label: int


def label_to_binary(label: str) -> int | None:
    if label == "P_LP":
        return 1
    if label == "B_LB":
        return 0
    return None


def make_binary_examples(
    score_rows: list[dict[str, str]],
    label_rows: list[dict[str, str]],
    score_column: str,
) -> list[BinaryExample]:
    labels = {row["variant_id"]: label_to_binary(row.get("clinvar_label", "")) for row in label_rows}
    examples: list[BinaryExample] = []
    for row in score_rows:
        variant_id = row["variant_id"]
        label = labels.get(variant_id)
        if label is None:
            continue
        value = row.get(score_column, "")
        if value == "":
            continue
        examples.append(BinaryExample(variant_id=variant_id, score=float(value), label=label))
    return examples


def auroc(examples: list[BinaryExample]) -> float:
    positives = [example for example in examples if example.label == 1]
    negatives = [example for example in examples if example.label == 0]
    if not positives or not negatives:
        return float("nan")

    wins = 0.0
    total = len(positives) * len(negatives)
    for pos in positives:
        for neg in negatives:
            if pos.score > neg.score:
                wins += 1.0
            elif pos.score == neg.score:
                wins += 0.5
    return wins / total


def average_precision(examples: list[BinaryExample]) -> float:
    positives = sum(example.label for example in examples)
    if positives == 0:
        return float("nan")

    ranked = sorted(examples, key=lambda example: example.score, reverse=True)
    precision_sum = 0.0
    tp = 0
    for index, example in enumerate(ranked, start=1):
        if example.label == 1:
            tp += 1
            precision_sum += tp / index
    return precision_sum / positives


def top_k_enrichment(examples: list[BinaryExample], percentile: float) -> float:
    if not examples:
        return float("nan")
    ranked = sorted(examples, key=lambda example: example.score, reverse=True)
    cutoff = max(1, round(len(ranked) * percentile))
    top = ranked[:cutoff]
    baseline = sum(example.label for example in examples) / len(examples)
    observed = sum(example.label for example in top) / len(top)
    if baseline == 0:
        return float("nan")
    return observed / baseline


def metrics_rows(examples: list[BinaryExample]) -> list[dict[str, object]]:
    n_pos = sum(example.label for example in examples)
    n_neg = len(examples) - n_pos
    return [
        {"metric": "n_examples", "value": len(examples)},
        {"metric": "n_pathogenic", "value": n_pos},
        {"metric": "n_benign", "value": n_neg},
        {"metric": "AUROC", "value": auroc(examples)},
        {"metric": "AUPRC", "value": average_precision(examples)},
        {"metric": "top_1pct_enrichment", "value": top_k_enrichment(examples, 0.01)},
        {"metric": "top_5pct_enrichment", "value": top_k_enrichment(examples, 0.05)},
        {"metric": "top_10pct_enrichment", "value": top_k_enrichment(examples, 0.10)},
    ]

