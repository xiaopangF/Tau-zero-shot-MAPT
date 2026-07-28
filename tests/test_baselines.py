from mapt_zero_shot.baselines import heuristic_score_row


def test_heuristic_scores_tau_sensitive_features():
    row = heuristic_score_row(
        {
            "variant_id": "P301L",
            "tau_region": "microtubule_repeat_R2_exon10",
            "tau_motif": "",
            "near_known_pathogenic_hotspot_3aa": "True",
            "near_ptm_site_3aa": "False",
            "charge_change": "none",
            "special_residue_change": "P_loss",
        }
    )
    assert row["heuristic_score"] == 6.25
    assert row["model"] == "tau_heuristic_v1"
    assert "hotspot_3aa=2" in row["heuristic_components"]


def test_heuristic_scores_low_for_unremarkable_n_terminal_variant():
    row = heuristic_score_row(
        {
            "variant_id": "M1A",
            "tau_region": "N_terminal_projection",
            "tau_motif": "",
            "near_known_pathogenic_hotspot_3aa": "False",
            "near_ptm_site_3aa": "False",
            "charge_change": "none",
            "special_residue_change": "",
        }
    )
    assert row["heuristic_score"] == 0.25
