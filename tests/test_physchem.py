from mapt_zero_shot.physchem import (
    PhyschemWeights,
    feature_row,
    feature_rows,
    grid_search_weights,
    ranked_physchem_rows,
    ranked_zero_shot_physchem_rows,
    summary_rows,
    zero_shot_physchem_rows,
    zero_shot_summary_rows,
)


def test_feature_row_computes_three_requested_deltas():
    row = {"variant_id": "P301L", "position": "301", "wt_aa": "P", "mut_aa": "L"}

    scored = feature_row(row)

    assert scored["hydrophobicity_delta"] == 5.4
    assert scored["beta_sheet_delta"] == 0.75
    assert scored["net_charge_delta"] == 0.0


def test_grid_search_prefers_weights_that_rank_controls_high():
    variants = [
        {"variant_id": "A1C", "position": "1", "wt_aa": "A", "mut_aa": "C"},
        {"variant_id": "P301L", "position": "301", "wt_aa": "P", "mut_aa": "L"},
        {"variant_id": "A2S", "position": "2", "wt_aa": "A", "mut_aa": "S"},
    ]
    rows = feature_rows(variants)

    result = grid_search_weights(
        rows,
        controls=("P301L",),
        grid_min=-1,
        grid_max=1,
        grid_step=1,
        top_fraction=1 / 3,
    )

    assert result.controls_in_top_fraction == 1
    assert result.control_ranks["P301L"] == 1


def test_ranked_physchem_rows_and_summary_include_gold_standard_metrics():
    variants = [
        {"variant_id": "P301L", "position": "301", "wt_aa": "P", "mut_aa": "L"},
        {"variant_id": "A2S", "position": "2", "wt_aa": "A", "mut_aa": "S"},
    ]
    rows = feature_rows(variants)
    weights = PhyschemWeights(
        hydrophobicity_delta=1.0,
        beta_sheet_delta=0.0,
        net_charge_delta=0.0,
    )

    ranked = ranked_physchem_rows(rows, weights, controls=("P301L",))
    result = grid_search_weights(
        rows,
        controls=("P301L",),
        grid_min=1,
        grid_max=1,
        grid_step=1,
        top_fraction=0.5,
    )
    summary = {row["metric"]: row["value"] for row in summary_rows(result, 1, 1, 1, 0.5)}

    assert ranked[0]["variant_id"] == "P301L"
    assert ranked[0]["gold_standard_control"] == "True"
    assert summary["gold_standard_mean_rank"] == 1
    assert summary["rank_P301L"] == 1


def test_zero_shot_model_uses_absolute_standardized_perturbation_and_region_prior():
    variants = [
        {"variant_id": "A1C", "position": "1", "wt_aa": "A", "mut_aa": "C"},
        {"variant_id": "P301L", "position": "301", "wt_aa": "P", "mut_aa": "L"},
        {"variant_id": "R406W", "position": "406", "wt_aa": "R", "mut_aa": "W"},
    ]

    scored = zero_shot_physchem_rows(feature_rows(variants))
    by_id = {str(row["variant_id"]): row for row in scored}

    assert by_id["P301L"]["tau_region"] == "microtubule_repeat_R2_exon10"
    assert by_id["P301L"]["region_multiplier"] == 1.45
    assert by_id["P301L"]["physchem_perturbation_score"] >= 0
    assert by_id["R406W"]["net_charge_delta"] == -1.0


def test_zero_shot_ranking_only_uses_controls_for_validation():
    variants = [
        {"variant_id": "A1C", "position": "1", "wt_aa": "A", "mut_aa": "C"},
        {"variant_id": "P301L", "position": "301", "wt_aa": "P", "mut_aa": "L"},
        {"variant_id": "R406W", "position": "406", "wt_aa": "R", "mut_aa": "W"},
    ]
    scored = zero_shot_physchem_rows(feature_rows(variants))
    ranked = ranked_zero_shot_physchem_rows(scored, controls=("P301L", "R406W"))
    summary = {
        row["metric"]: row["value"]
        for row in zero_shot_summary_rows(ranked, controls=("P301L", "R406W"), top_fraction=1 / 3)
    }

    assert {row["variant_id"] for row in ranked} == {"A1C", "P301L", "R406W"}
    assert summary["model_type"] == "zero_shot_physchem_perturbation"
    assert summary["uses_gold_standard_labels_for_scoring"] == "False"
    assert summary["gold_standard_count"] == 2