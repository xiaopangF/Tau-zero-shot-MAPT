$ErrorActionPreference = "Stop"

$texDir = "manuscript\latex"
$main = "main.tex"
$mainPdf = "main.pdf"
$miktexBin = Join-Path $env:LOCALAPPDATA "Programs\MiKTeX\miktex\bin\x64"

if (Test-Path $miktexBin) {
    $env:PATH = "$miktexBin;$env:PATH"
}

if (-not (Test-Path (Join-Path $texDir $main))) {
    throw "Could not find $texDir\$main"
}

Push-Location $texDir
try {
    if (Get-Command pdflatex -ErrorAction SilentlyContinue) {
        pdflatex -interaction=nonstopmode -halt-on-error $main
        bibtex main
        pdflatex -interaction=nonstopmode -halt-on-error $main
        pdflatex -interaction=nonstopmode -halt-on-error $main
    } elseif (Get-Command tectonic -ErrorAction SilentlyContinue) {
        tectonic $main
    } elseif (Get-Command latexmk -ErrorAction SilentlyContinue) {
        latexmk -pdf -interaction=nonstopmode -halt-on-error $main
    } else {
        throw "No LaTeX engine found. Install MiKTeX, TeX Live, TinyTeX, or Tectonic first."
    }

    if (-not (Test-Path $mainPdf)) {
        throw "LaTeX finished without creating $texDir\$mainPdf"
    }
    Write-Host "Wrote $texDir\$mainPdf"
} finally {
    Pop-Location
}