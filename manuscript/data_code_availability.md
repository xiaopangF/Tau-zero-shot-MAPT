# Data and Code Availability Draft

## Code availability

All source code for the MAPT/Tau zero-shot atlas workflow is available at:

`https://github.com/xiaopangF/Tau-zero-shot-MAPT`

The repository contains the command-line workflow for Tau-F missense variant
enumeration, Tau-specific annotation, ESM scoring, ESM-1v ensemble averaging,
ClinVar import and coordinate QC, AlphaMissense import and coordinate QC,
heuristic baseline scoring, model comparison, concordance analysis,
prioritization, and manuscript figure generation.

## Data availability

Large raw and generated data files are not committed directly to GitHub. They can
be regenerated from public resources and the repository workflow.

Public input resources used in the current workflow include:

- NCBI ClinVar `variant_summary.txt.gz`.
- AlphaMissense public hg38 score table: `AlphaMissense_hg38.tsv.gz`.
- Public ESM model checkpoints downloaded through the `fair-esm` / PyTorch model
  loading workflow.

Generated local files include:

- `data/processed/mapt_all_missense_variants.tsv`
- `data/processed/mapt_clinvar_benchmark.tsv`
- `data/processed/mapt_clinvar_rejected.tsv`
- `data/processed/mapt_alphamissense.tsv`
- `data/processed/mapt_alphamissense_rejected.tsv`
- `results/scores/mapt_esm1v_ensemble.tsv`
- `results/mapt_esm1v_ensemble_domain_summary.tsv`
- `results/mapt_esm1v_top50_vus_priority.tsv`
- `results/mapt_model_concordance.tsv`
- `results/mapt_model_concordance_summary.tsv`
- `results/mapt_model_concordance_top.tsv`
- `results/manuscript_assets/`

For journal submission, these generated outputs should be deposited as
supplementary files or in a data repository such as Zenodo, Figshare, or an
institutional repository. The GitHub repository should remain the source for code
and reproducibility instructions, while large generated tables should be provided
through the manuscript's data-deposition record.

## Reproducibility summary

A new user can reproduce the major workflow by:

1. Creating the Python environment described in the README.
2. Downloading ClinVar and AlphaMissense public input files.
3. Enumerating the 8379 Tau-F missense variants.
4. Importing ClinVar and AlphaMissense with strict coordinate QC.
5. Running or reusing the ESM-1v ensemble scoring workflow.
6. Generating summaries, figures, VUS rankings, and concordance tables.

The main manuscript claim does not depend on training a clinical classifier.
Clinical labels are used only for QC, smoke-check benchmarking, and
interpretation.

## Current limitation for final publication

Before journal submission, the generated result tables and manuscript figures
should be archived with a permanent DOI. The current GitHub repository alone is
not ideal for hosting large generated TSV files because the project intentionally
ignores `data/raw`, `data/external`, `data/processed`, and `results` outputs.