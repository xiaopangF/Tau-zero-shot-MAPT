from mapt_zero_shot.compare import compare_model_rows


def test_compare_model_rows():
    labels = [
        {"variant_id": "P301L", "clinvar_label": "P_LP"},
        {"variant_id": "M1A", "clinvar_label": "B_LB"},
    ]
    rows = compare_model_rows(
        [
            ("good", [{"variant_id": "P301L", "s": "2"}, {"variant_id": "M1A", "s": "0"}], "s"),
            ("bad", [{"variant_id": "P301L", "s": "0"}, {"variant_id": "M1A", "s": "2"}], "s"),
        ],
        labels,
    )
    aurocs = {row["model"]: row["value"] for row in rows if row["metric"] == "AUROC"}
    assert aurocs["good"] == 1.0
    assert aurocs["bad"] == 0.0
