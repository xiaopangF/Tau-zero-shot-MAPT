param()

$ErrorActionPreference = "Stop"
$PSNativeCommandUseErrorActionPreference = $true

$checkpointDir = Join-Path $env:USERPROFILE ".cache\torch\hub\checkpoints"
New-Item -ItemType Directory -Force $checkpointDir | Out-Null

$models = @(
    "esm1v_t33_650M_UR90S_1",
    "esm1v_t33_650M_UR90S_2",
    "esm1v_t33_650M_UR90S_3",
    "esm1v_t33_650M_UR90S_4",
    "esm1v_t33_650M_UR90S_5"
)

foreach ($model in $models) {
    $url = "https://dl.fbaipublicfiles.com/fair-esm/models/$model.pt"
    $out = Join-Path $checkpointDir "$model.pt"
    $head = curl.exe -L -I $url
    $contentLengthLine = $head | Where-Object { $_ -match '^Content-Length:' } | Select-Object -Last 1
    if (-not $contentLengthLine) {
        throw "Could not determine Content-Length for $url"
    }
    $expected = [int64](($contentLengthLine -split ':', 2)[1].Trim())
    $actual = if (Test-Path $out) { (Get-Item $out).Length } else { 0 }
    if ($actual -eq $expected) {
        Write-Host "$model already complete ($actual bytes)."
        continue
    }
    if ($actual -gt $expected) {
        Remove-Item -LiteralPath $out
    }
    Write-Host "Downloading $model to $out"
    curl.exe -L --retry 20 --retry-delay 5 --retry-all-errors -C - -o $out $url
    $actual = (Get-Item $out).Length
    if ($actual -ne $expected) {
        throw "$model download size mismatch: expected $expected, got $actual"
    }
}
