from mapt_zero_shot.annotations import annotate_variant_row, region_for_position


def test_region_annotations():
    assert region_for_position(301) == "microtubule_repeat_R2_exon10"
    assert region_for_position(406) == "C_terminal_tail"


def test_variant_annotation_fields():
    row = annotate_variant_row({"variant_id": "P301L", "position": 301, "wt_aa": "P", "mut_aa": "L"})
    assert row["near_known_pathogenic_hotspot_3aa"] == "True"
    assert row["special_residue_change"] == "P_loss"

