param(
    [string]$Device = "auto",
    [switch]$AllowCpu
)

$ErrorActionPreference = "Stop"
$env:PYTHONPATH = "src"

$hasCuda = python -c "import torch; print('1' if torch.cuda.is_available() else '0')"
if ($Device -eq "auto") {
    if ($hasCuda.Trim() -eq "1") {
        $ResolvedDevice = "cuda"
    } else {
        $ResolvedDevice = "cpu"
    }
} else {
    $ResolvedDevice = $Device
}

if ($ResolvedDevice -eq "cpu" -and -not $AllowCpu) {
    throw "CUDA was not detected. ESM-1v is large; rerun with -AllowCpu only if you intentionally want a long CPU job."
}

New-Item -ItemType Directory -Force results/scores | Out-Null

$models = @(
    "esm1v_t33_650M_UR90S_1",
    "esm1v_t33_650M_UR90S_2",
    "esm1v_t33_650M_UR90S_3",
    "esm1v_t33_650M_UR90S_4",
    "esm1v_t33_650M_UR90S_5"
)

$scoreFiles = @()
for ($i = 0; $i -lt $models.Length; $i++) {
    $model = $models[$i]
    $out = "results/scores/mapt_esm1v_$($i + 1).tsv"
    $scoreFiles += $out
    if (Test-Path $out) {
        Write-Host "Skipping existing $out"
        continue
    }
    Write-Host "Scoring $model on $ResolvedDevice..."
    python -m mapt_zero_shot.cli score-esm --model $model --device $ResolvedDevice --out $out
}

python -m mapt_zero_shot.cli ensemble --scores $scoreFiles --score-column pathogenic_score --out results/scores/mapt_esm1v_ensemble.tsv

if (Test-Path data/processed/mapt_clinvar_benchmark.tsv) {
    python -m mapt_zero_shot.cli evaluate --scores results/scores/mapt_esm1v_ensemble.tsv --labels data/processed/mapt_clinvar_benchmark.tsv --score-column pathogenic_score_mean --out results/mapt_esm1v_ensemble_metrics.tsv
    python -m mapt_zero_shot.cli prioritize --scores results/scores/mapt_esm1v_ensemble.tsv --annotations data/processed/mapt_clinvar_benchmark.tsv --score-column pathogenic_score_mean --label VUS --limit 50 --out results/mapt_esm1v_top50_vus_priority.tsv
    python -m mapt_zero_shot.cli figures --scores results/scores/mapt_esm1v_ensemble.tsv --labels data/processed/mapt_clinvar_benchmark.tsv --score-column pathogenic_score_mean --outdir results/figures/esm1v_ensemble
}
