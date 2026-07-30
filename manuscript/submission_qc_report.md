# Submission QC Report

This report checks whether the manuscript draft is internally consistent with
current generated result tables and whether key submission files exist.

## Summary

- pass: 14
- warn: 0
- fail: 0

## Checks

| Check | Status | Detail |
|---|---|---|
| variant_count | pass | expected manuscript to contain 8379 |
| clinvar_accepted | pass | expected manuscript to contain 33 |
| clinvar_rejected | pass | expected manuscript to contain 1119 |
| alphamissense_scored | pass | expected manuscript to contain 877 |
| concordance_esm_heuristic_high | pass | expected manuscript to contain 209 |
| concordance_esm_only_high | pass | expected manuscript to contain 629 |
| concordance_heuristic_only_high | pass | expected manuscript to contain 629 |
| risky_phrase:accurately predicts pathogenicity | pass | occurrences=0; risky_contexts=0 |
| risky_phrase:clinically validates | pass | occurrences=0; risky_contexts=0 |
| risky_phrase:clinical diagnostic model | pass | occurrences=3; risky_contexts=0 |
| risky_phrase:diagnostic test | pass | occurrences=1; risky_contexts=0 |
| risky_phrase:reclassified pathogenic | pass | occurrences=0; risky_contexts=0 |
| file_exists:manuscript/latex/main.pdf | pass | present |
| file_exists:results/manuscript_assets/figure_positive_control_calibration.png | pass | present |
