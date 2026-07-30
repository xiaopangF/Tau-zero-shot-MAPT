from mapt_zero_shot.physchem import (
    PhyschemWeights,
    feature_row,
    feature_rows,
    grid_search_weights,
    ranked_physchem_rows,
    summary_rows,
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
