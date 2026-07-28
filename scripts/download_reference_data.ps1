param(
    [switch]$AlphaMissense
)

$ErrorActionPreference = "Stop"
New-Item -ItemType Directory -Force data/raw, data/external | Out-Null

$clinvarUrl = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
$clinvarOut = "data/raw/variant_summary.txt.gz"
Write-Host "Downloading ClinVar variant summary..."
Invoke-WebRequest -Uri $clinvarUrl -OutFile $clinvarOut
Write-Host "Wrote $clinvarOut"

if ($AlphaMissense) {
    $alphaUrl = "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz"
    $alphaOut = "data/external/AlphaMissense_hg38.tsv.gz"
    Write-Host "Downloading AlphaMissense hg38 scores. This file is large."
    Invoke-WebRequest -Uri $alphaUrl -OutFile $alphaOut
    Write-Host "Wrote $alphaOut"
}
