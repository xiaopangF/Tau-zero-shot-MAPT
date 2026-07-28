from mapt_zero_shot.evaluate import auroc, average_precision, make_binary_examples


def test_metrics_perfect_ranking():
    examples = make_binary_examples(
        [
            {"variant_id": "A1C", "score": "2.0"},
            {"variant_id": "A1D", "score": "1.5"},
            {"variant_id": "A1E", "score": "-1.0"},
        ],
        [
            {"variant_id": "A1C", "clinvar_label": "P_LP"},
            {"variant_id": "A1D", "clinvar_label": "P_LP"},
            {"variant_id": "A1E", "clinvar_label": "B_LB"},
        ],
        "score",
    )
    assert auroc(examples) == 1.0
    assert average_precision(examples) == 1.0

