import gzip

from mapt_zero_shot.clinvar import load_mapt_clinvar_with_qc


def test_clinvar_reference_wt_qc(tmp_path):
    path = tmp_path / "variant_summary.txt.gz"
    header = ["GeneSymbol", "VariationID", "Name", "ClinicalSignificance", "ReviewStatus"]
    rows = [
        ["MAPT", "1", "NM_test(MAPT):c.902C>T (p.Pro301Leu)", "Pathogenic", "criteria provided"],
        ["MAPT", "2", "NM_test(MAPT):c.914A>G (p.Gln305Arg)", "Benign", "criteria provided"],
    ]
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(row) + "\n")

    result = load_mapt_clinvar_with_qc(str(path))
    assert "P301L" in result.accepted_by_variant
    assert "Q305R" not in result.accepted_by_variant
    assert result.rejected_rows[0]["reject_reason"] == "reference_wt_mismatch"
    assert result.rejected_rows[0]["reference_wt_aa"] == "S"
