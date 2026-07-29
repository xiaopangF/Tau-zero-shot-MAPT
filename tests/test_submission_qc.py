from mapt_zero_shot.submission_qc import submission_qc_rows


def test_submission_qc_checks_key_numbers_and_required_files(tmp_path):
    required = tmp_path / "exists.txt"
    required.write_text("ok", encoding="utf-8")
    rows = submission_qc_rows(
        "This draft mentions 2 variants, 1 accepted, 3 rejected, 4 scored, 5 and 6.",
        variants=[{"variant_id": "A1C"}, {"variant_id": "A2C"}],
        clinvar_summary=[
            {"section": "accepted", "category": "total_annotated_variants", "count": "1"},
            {"section": "rejected", "category": "total_rejected_rows", "count": "3"},
        ],
        alphamissense_summary=[
            {"section": "coverage", "category": "scored_variants", "count": "4"}
        ],
        concordance_summary=[
            {"section": "concordance_category", "category": "esm_heuristic_high", "count": "5"},
            {"section": "concordance_category", "category": "esm_only_high", "count": "6"},
            {"section": "concordance_category", "category": "heuristic_only_high", "count": "7"},
        ],
        required_files=[required, tmp_path / "missing.txt"],
    )
    keyed = {row["check"]: row["status"] for row in rows}
    assert keyed["variant_count"] == "pass"
    assert keyed["concordance_heuristic_only_high"] == "fail"
    assert keyed[f"file_exists:{required.as_posix()}"] == "pass"
    assert keyed[f"file_exists:{(tmp_path / 'missing.txt').as_posix()}"] == "fail"