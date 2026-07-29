from pathlib import Path


def test_submission_package_script_lists_core_outputs():
    script = Path("scripts/prepare_submission_package.ps1").read_text(encoding="utf-8")
    assert "supplementary_table_1_full_missense_atlas.tsv" in script
    assert "supplementary_table_6a_model_concordance_all.tsv" in script
    assert "figure_1_workflow_schematic.png" in script
    assert "manifest.tsv" in script