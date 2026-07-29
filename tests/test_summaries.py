from mapt_zero_shot.summaries import (
    alphamissense_summary_rows,
    clinvar_summary_rows,
    domain_summary_rows,
)


def test_clinvar_summary_counts_labels_and_rejections():
    rows = clinvar_summary_rows(
        [
            {"variant_id": "P301L", "clinvar_label": "P_LP"},
            {"variant_id": "A2S", "clinvar_label": "VUS"},
            {"variant_id": "M1A", "clinvar_label": ""},
        ],
        [
            {"reject_reason": "outside_reference_range", "clinvar_label": "P_LP"},
            {"reject_reason": "outside_reference_range", "clinvar_label": "VUS"},
            {"reject_reason": "reference_wt_mismatch", "clinvar_label": "B_LB"},
        ],
    )
    keyed = {(row["section"], row["category"]): row["count"] for row in rows}
    assert keyed[("accepted", "total_annotated_variants")] == 2
    assert keyed[("accepted_label", "P_LP")] == 1
    assert keyed[("rejected_reason", "outside_reference_range")] == 2


def test_domain_summary_counts_top_fraction():
    score_rows = [
        {"variant_id": "A1C", "pathogenic_score": "10"},
        {"variant_id": "A2C", "pathogenic_score": "1"},
        {"variant_id": "A3C", "pathogenic_score": "5"},
    ]
    annotation_rows = [
        {"variant_id": "A1C", "tau_region": "r1"},
        {"variant_id": "A2C", "tau_region": "r1"},
        {"variant_id": "A3C", "tau_region": "r2"},
    ]
    rows = domain_summary_rows(score_rows, annotation_rows, "pathogenic_score")
    by_region = {row["tau_region"]: row for row in rows}
    assert by_region["r1"]["n_variants"] == 2
    assert by_region["r1"]["max_score"] == 10.0
    assert by_region["r2"]["median_score"] == 5.0


def test_alphamissense_summary_counts_coverage_and_rejections():
    rows = alphamissense_summary_rows(
        [
            {
                "variant_id": "A2S",
                "position": "2",
                "tau_region": "N_terminal_projection",
                "alphamissense_score": "0.13",
                "alphamissense_class": "likely_benign",
            },
            {
                "variant_id": "A2T",
                "position": "2",
                "tau_region": "N_terminal_projection",
                "alphamissense_score": "0.14",
                "alphamissense_class": "likely_benign",
            },
            {
                "variant_id": "P301L",
                "position": "301",
                "tau_region": "microtubule_repeat_R2_exon10",
            },
        ],
        [
            {"reject_reason": "reference_wt_mismatch"},
            {"reject_reason": "outside_reference_range"},
            {"reject_reason": "outside_reference_range"},
        ],
    )
    keyed = {(row["section"], row["category"]): row["count"] for row in rows}
    assert keyed[("coverage", "total_variants")] == 3
    assert keyed[("coverage", "scored_variants")] == 2
    assert keyed[("coverage", "scored_positions")] == 1
    assert keyed[("scored_region", "N_terminal_projection")] == 2
    assert keyed[("scored_class", "likely_benign")] == 2
    assert keyed[("rejected_reason", "outside_reference_range")] == 2