from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_latex_reading_version_includes_core_assets_and_citations():
    tex = (ROOT / "manuscript" / "latex" / "main.tex").read_text(encoding="utf-8")

    assert r"\bibliography{../references}" in tex
    assert "figure_workflow_schematic.png" in tex
    assert "missense_heatmap.png" in tex
    assert "figure_vus_lollipop.png" in tex
    assert r"\citep{meier2021language}" in tex


def test_latex_build_script_supports_local_miktex_and_bibtex():
    script = (ROOT / "scripts" / "build_latex_pdf.ps1").read_text(encoding="utf-8")

    assert "MiKTeX" in script
    assert "pdflatex" in script
    assert "bibtex main" in script
    assert "main.pdf" in script
