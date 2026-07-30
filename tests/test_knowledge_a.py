from mapt_zero_shot.knowledge_a import (
    directed_beta_disruption_score,
    microtubule_interface_score,
    ptm_microenvironment_score,
    score_scheme_a_rows,
    scheme_a_gold_rows,
    scheme_a_summary_rows,
    splicing_proximity_score,
)
from mapt_zero_shot.physchem import feature_row, feature_rows


def test_scheme_a_directed_beta_score_requires_r_region_and_beta_gain():
    p301l = feature_row({"variant_id": "P301L", "position": "301", "wt_aa": "P", "mut_aa": "L"})
    n100v = feature_row({"variant_id": "N100V", "position": "100", "wt_aa": "N", "mut_aa": "V"})

    assert directed_beta_disruption_score(p301l) == 2.5
    assert directed_beta_disruption_score(n100v) == 0.0


def test_scheme_a_ptm_score_stacks_charge_and_hydrophobicity_rules():
    row = feature_row({"variant_id": "R202W", "position": "202", "wt_aa": "R", "mut_aa": "W"})

    assert ptm_microenvironment_score(row) == 3.5


def test_scheme_a_splicing_score_uses_distance_to_285():
    assert splicing_proximity_score({"position": "285"}) == 3.0
    assert splicing_proximity_score({"position": "279"}) == 2.0
    assert splicing_proximity_score({"position": "260"}) == 1.0
    assert splicing_proximity_score({"position": "259"}) == 0.0


def test_scheme_a_microtubule_interface_score_uses_fixed_intervals_once():
    assert microtubule_interface_score({"position": "255"}) == 1.5
    assert microtubule_interface_score({"position": "285"}) == 1.5
    assert microtubule_interface_score({"position": "337"}) == 0.0
    assert microtubule_interface_score({"position": "346"}) == 0.0


def test_score_scheme_a_rows_outputs_final_rank_and_gold_details():
    variants = [
        {"variant_id": "A1C", "protein_change": "p.A1C", "position": "1", "wt_aa": "A", "mut_aa": "C"},
        {"variant_id": "P301L", "protein_change": "p.P301L", "position": "301", "wt_aa": "P", "mut_aa": "L"},
        {"variant_id": "R202W", "protein_change": "p.R202W", "position": "202", "wt_aa": "R", "mut_aa": "W"},
    ]
    rows = feature_rows(variants)

    ranked = score_scheme_a_rows(rows, controls=("P301L",))
    gold = scheme_a_gold_rows(ranked, controls=("P301L",))
    summary = {
        row["metric"]: row["value"]
        for row in scheme_a_summary_rows(ranked, controls=("P301L",), top_fraction=1 / 3)
    }

    assert len(ranked) == 3
    assert gold[0]["variant_id"] == "P301L"
    assert gold[0]["directed_beta_score"] == 2.5
    assert "scheme_a_final_score" in gold[0]
    assert "scheme_a_rank" in gold[0]
    assert summary["model_type"] == "scheme_a_knowledge_driven_with_interface"
    assert summary["uses_gold_standard_labels_for_scoring"] == "False"