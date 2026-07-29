from mapt_zero_shot.concordance import (
    concordance_rows,
    concordance_summary_rows,
    top_concordance_rows,
)


def test_concordance_rows_classify_top_fraction_overlap():
    annotations = [
        {"variant_id": "A1C", "position": "1"},
        {"variant_id": "A2C", "position": "2"},
        {"variant_id": "A3C", "position": "3"},
        {"variant_id": "A4C", "position": "4"},
    ]
    esm = [
        {"variant_id": "A1C", "esm": "10"},
        {"variant_id": "A2C", "esm": "9"},
        {"variant_id": "A3C", "esm": "1"},
        {"variant_id": "A4C", "esm": "0"},
    ]
    heuristic = [
        {"variant_id": "A1C", "h": "10"},
        {"variant_id": "A2C", "h": "1"},
        {"variant_id": "A3C", "h": "9"},
        {"variant_id": "A4C", "h": "0"},
    ]
    alpha = [
        {"variant_id": "A1C", "a": "10"},
        {"variant_id": "A2C", "a": "1"},
        {"variant_id": "A3C", "a": "0"},
    ]

    rows = concordance_rows(
        annotations,
        esm,
        heuristic,
        alpha,
        esm_column="esm",
        heuristic_column="h",
        alphamissense_column="a",
        top_fraction=0.5,
    )
    by_variant = {row["variant_id"]: row for row in rows}
    assert by_variant["A1C"]["concordance_category"] == "all_three_high"
    assert by_variant["A2C"]["concordance_category"] == "esm_alphamissense_high"
    assert by_variant["A3C"]["concordance_category"] == "heuristic_only_high"
    assert by_variant["A4C"]["concordance_category"] == "no_model_high"
    assert by_variant["A4C"]["alphamissense_has_score"] is False

    summary = concordance_summary_rows(rows)
    keyed = {(row["section"], row["category"]): row["count"] for row in summary}
    assert keyed[("coverage", "total_variants")] == 4
    assert keyed[("coverage", "alphamissense_scored")] == 3
    assert keyed[("concordance_category", "all_three_high")] == 1

    top = top_concordance_rows(rows, limit=2)
    assert [row["variant_id"] for row in top] == ["A1C", "A2C"]