# LaTeX reading version

This directory contains a cleaner LaTeX version of the MAPT/Tau manuscript for
reading and sharing.

Main file:

- `main.tex`

Build script from the repository root:

```powershell
.\scripts\build_latex_pdf.ps1
```

The script needs one LaTeX engine installed, such as MiKTeX, TeX Live, TinyTeX,
or Tectonic. The current Windows machine did not have `xelatex`, `pdflatex`,
`latexmk`, or `pandoc` installed when this folder was created.

The LaTeX file references figures from:

- `results/manuscript_assets/`
- `results/figures/esm1v_ensemble/`

Those generated figure files are local outputs and are not committed to GitHub.