from mapt_zero_shot.calibration import positive_control_rows, positive_control_summary_rows


def test_positive_controls_receive_global_rank_and_cutoff_flags():
    score_rows = [
        {"variant_id": "A1C", "pathogenic_score_mean": "10"},
        {"variant_id": "P301L", "pathogenic_score_mean": "8"},
        {"variant_id": "G272V", "pathogenic_score_mean": "3"},
        {"variant_id": "R406W", "pathogenic_score_mean": "1"},
    ]
    annotation_rows = [
        {"variant_id": "P301L", "position": "301", "tau_region": "R2"},
        {"variant_id": "G272V", "position": "272", "tau_region": "R1"},
        {"variant_id": "R406W", "position": "406", "tau_region": "tail"},
    ]
    controls = (
        {"variant_id": "P301L", "evidence_class": "pathogenic", "evidence_note": "test"},
    )

    rows = positive_control_rows(score_rows, annotation_rows, controls=controls)

    assert rows[0]["rank"] == 2
    assert rows[0]["atlas_size"] == 4
    assert rows[0]["top_10pct"] == "False"
    assert rows[0]["tau_region"] == "R2"


def test_positive_control_summary_is_explicit_about_top_decile_count():
    rows = [
        {
            "found_in_atlas": "True",
            "atlas_size": 100,
            "top_fraction": 0.03,
            "top_1pct": "False",
            "top_5pct": "True",
            "top_10pct": "True",
        },
        {
            "found_in_atlas": "True",
            "atlas_size": 100,
            "top_fraction": 0.20,
            "top_1pct": "False",
            "top_5pct": "False",
            "top_10pct": "False",
        },
    ]

    summary = {row["metric"]: row["value"] for row in positive_control_summary_rows(rows)}

    assert summary["controls_found"] == 2
    assert summary["top_5pct_count"] == 1
    assert summary["top_10pct_count"] == 1
