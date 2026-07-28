"""Compare multiple model score files against the same benchmark labels."""

from __future__ import annotations

from .evaluate import make_binary_examples, metrics_rows


def compare_model_rows(
    model_specs: list[tuple[str, list[dict[str, str]], str]],
    label_rows: list[dict[str, str]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for model_name, score_rows, score_column in model_specs:
        examples = make_binary_examples(score_rows, label_rows, score_column)
        for metric_row in metrics_rows(examples):
            rows.append(
                {
                    "model": model_name,
                    "score_column": score_column,
                    "metric": metric_row["metric"],
                    "value": metric_row["value"],
                }
            )
    return rows
