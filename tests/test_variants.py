from mapt_zero_shot.constants import MAPT_441_SEQUENCE
from mapt_zero_shot.variants import generate_missense_variants, parse_protein_change


def test_reference_length():
    assert len(MAPT_441_SEQUENCE) == 441


def test_generate_all_missense_variants():
    variants = generate_missense_variants()
    assert len(variants) == 441 * 19
    assert variants[0].variant_id == "M1A"
    assert all(variant.wt_aa != variant.mut_aa for variant in variants)


def test_parse_protein_change_one_letter():
    parsed = parse_protein_change("p.R406W")
    assert parsed is not None
    assert parsed.variant_id == "R406W"


def test_parse_protein_change_three_letter():
    parsed = parse_protein_change("p.Pro301Leu")
    assert parsed is not None
    assert parsed.variant_id == "P301L"
