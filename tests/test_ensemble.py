from mapt_zero_shot.ensemble import ensemble_score_rows


def test_ensemble_mean_and_std():
    rows = ensemble_score_rows(
        [
            [{"variant_id": "P301L", "protein_change": "p.P301L", "position": "301", "wt_aa": "P", "mut_aa": "L", "model": "m1", "pathogenic_score": "2"}],
            [{"variant_id": "P301L", "protein_change": "p.P301L", "position": "301", "wt_aa": "P", "mut_aa": "L", "model": "m2", "pathogenic_score": "4"}],
        ]
    )
    assert len(rows) == 1
    assert rows[0]["pathogenic_score_mean"] == 3.0
    assert rows[0]["n_models"] == 2
    assert rows[0]["models"] == "m1;m2"
