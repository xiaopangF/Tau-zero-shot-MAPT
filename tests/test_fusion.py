from mapt_zero_shot.fusion import composite_summary_rows, composite_zero_shot_rows


def test_composite_zero_shot_rows_standardizes_and_ranks_shared_variants():
    esm_rows = [
        {"variant_id": "A1C", "protein_change": "p.A1C", "position": "1", "wt_aa": "A", "mut_aa": "C", "pathogenic_score_mean": "1"},
        {"variant_id": "P301L", "protein_change": "p.P301L", "position": "301", "wt_aa": "P", "mut_aa": "L", "pathogenic_score_mean": "3"},
    ]
    physchem_rows = [
        {"variant_id": "A1C", "physchem_zero_shot_score": "1", "tau_region": "N_terminal_projection"},
        {"variant_id": "P301L", "physchem_zero_shot_score": "3", "tau_region": "microtubule_repeat_R2_exon10"},
    ]

    ranked = composite_zero_shot_rows(esm_rows, physchem_rows, esm_weight=3, physchem_weight=1)

    assert ranked[0]["variant_id"] == "P301L"
    assert ranked[0]["zero_shot_composite_rank"] == 1
    assert ranked[0]["esm_weight"] == 0.75
    assert ranked[0]["physchem_weight"] == 0.25
    assert ranked[0]["tau_region"] == "microtubule_repeat_R2_exon10"


def test_composite_summary_rows_reports_controls_without_training_signal():
    ranked = [
        {"variant_id": "P301L", "zero_shot_composite_rank": 1},
        {"variant_id": "A1C", "zero_shot_composite_rank": 2},
    ]

    summary = {
        row["metric"]: row["value"]
        for row in composite_summary_rows(ranked, controls=("P301L",), top_fraction=0.5)
    }

    assert summary["model_type"] == "zero_shot_esm_physchem_fusion"
    assert summary["uses_gold_standard_labels_for_scoring"] == "False"
    assert summary["gold_standard_in_top_fraction"] == 1
    assert summary["rank_P301L"] == 1