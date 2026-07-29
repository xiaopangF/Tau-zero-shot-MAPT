import gzip

from mapt_zero_shot.external_scores import load_alphamissense_with_qc


def test_alphamissense_reader_skips_comments_and_checks_reference(tmp_path):
    path = tmp_path / "alpha.tsv.gz"
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        handle.write("# Copyright notice\n")
        handle.write("# License notice\n")
        handle.write(
            "#CHROM\tPOS\tuniprot_id\ttranscript_id\tprotein_variant\t"
            "am_pathogenicity\tam_class\n"
        )
        handle.write("1\t1\tP10636\tTX1\tp.Pro301Leu\t0.99\tlikely_pathogenic\n")
        handle.write("1\t2\tP10636\tTX1\tp.Gln305Arg\t0.12\tlikely_benign\n")

    result = load_alphamissense_with_qc(str(path), uniprot_id="P10636")
    assert result.accepted_by_variant["P301L"]["alphamissense_score"] == "0.99"
    assert "Q305R" not in result.accepted_by_variant
    assert result.rejected_rows[0]["reject_reason"] == "reference_wt_mismatch"
    assert result.rejected_rows[0]["reference_wt_aa"] == "S"
