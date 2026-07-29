# Data directory

Large or externally downloaded files are intentionally not committed.

Recommended raw inputs:

- `data/raw/variant_summary.txt.gz`: NCBI ClinVar variant summary.
- `data/external/AlphaMissense_hg38.tsv.gz`: AlphaMissense public score table.
- Optional: gnomAD MAPT region export.
- Optional: UniProt MAPT TSV export.

Generated files:

- `data/processed/mapt_all_missense_variants.tsv`
- `data/processed/mapt_clinvar_benchmark.tsv`
- `data/processed/mapt_alphamissense.tsv`
- `data/processed/mapt_alphamissense_rejected.tsv`
