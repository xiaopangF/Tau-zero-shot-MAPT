$ErrorActionPreference = "Stop"
$env:PYTHONPATH = "src"

if (Test-Path ".\.venv\Scripts\python.exe") {
    $Python = ".\.venv\Scripts\python.exe"
} else {
    $Python = "python"
}

New-Item -ItemType Directory -Force results\manuscript_assets | Out-Null

& $Python -m mapt_zero_shot.cli calibrate-positive-controls `
    --scores results\scores\mapt_esm1v_ensemble.tsv `
    --annotations data\processed\mapt_all_missense_variants.tsv `
    --out results\mapt_esm1v_positive_controls.tsv `
    --summary-out results\mapt_esm1v_positive_control_summary.tsv `
    --figure results\manuscript_assets\figure_positive_control_calibration.png
