# Supplementary Tables Draft

## Supplementary Table 1. Full Tau-F missense atlas

Proposed file name: `supplementary_table_1_full_missense_atlas.tsv`

Current source file:

- `data/processed/mapt_all_missense_variants.tsv`

Description: complete list of all 8379 possible single amino-acid substitutions
in the 441-aa Tau-F / hTau40 reference sequence. Columns include variant ID,
protein change, position, reference amino acid, mutant amino acid, Tau region,
aggregation motif annotation, PTM proximity, pathogenic hotspot proximity,
charge-change class, and special-residue-change class.

Use in manuscript: supports the full atlas claim.

## Supplementary Table 2. ESM-1v ensemble scores

Proposed file name: `supplementary_table_2_esm1v_ensemble_scores.tsv`

Current source file:

- `results/scores/mapt_esm1v_ensemble.tsv`

Description: ESM-1v five-checkpoint ensemble scores for all 8379 Tau-F missense
variants. Columns include per-variant metadata, number of models, model names,
`pathogenic_score_mean`, and `pathogenic_score_std`.

Use in manuscript: supports all zero-shot score analyses.

## Supplementary Table 3. ClinVar coordinate QC

Proposed file names:

- `supplementary_table_3a_clinvar_accepted.tsv`
- `supplementary_table_3b_clinvar_rejected.tsv`
- `supplementary_table_3c_clinvar_qc_summary.tsv`

Current source files:

- `data/processed/mapt_clinvar_benchmark.tsv`
- `data/processed/mapt_clinvar_rejected.tsv`
- `results/mapt_clinvar_qc_summary.tsv`

Description: accepted ClinVar MAPT rows under strict 441-aa Tau-F reference
matching, rejected rows, and summary counts. Rejected rows are classified as
non-parseable missense, outside the Tau-F reference range, or reference amino-acid
mismatch.

Use in manuscript: supports the statement that direct clinical benchmarking is
limited after strict coordinate QC.

## Supplementary Table 4. AlphaMissense coordinate QC

Proposed file names:

- `supplementary_table_4a_alphamissense_accepted.tsv`
- `supplementary_table_4b_alphamissense_rejected.tsv`
- `supplementary_table_4c_alphamissense_qc_summary.tsv`

Current source files:

- `data/processed/mapt_alphamissense.tsv`
- `data/processed/mapt_alphamissense_rejected.tsv`
- `results/mapt_alphamissense_qc_summary.tsv`

Description: AlphaMissense scores accepted under strict Tau-F coordinate QC,
rejected AlphaMissense MAPT rows, and coverage summary. The current run accepted
877 scored Tau-F variants and rejected rows that were outside the 441-aa range or
had reference-residue mismatches.

Use in manuscript: supports AlphaMissense coverage and coordinate-mismatch
claims.

## Supplementary Table 5. Top ClinVar VUS priority list

Proposed file name: `supplementary_table_5_top_vus_priority.tsv`

Current source file:

- `results/mapt_esm1v_top50_vus_priority.tsv`

Description: top ClinVar VUS candidates ranked by ESM-1v ensemble
`pathogenic_score_mean`, with Tau region, mechanism summary, ClinVar review
status, and score.

Use in manuscript: supports the VUS prioritization section and Figure 5.

## Supplementary Table 6. Model concordance and discordance

Proposed file names:

- `supplementary_table_6a_model_concordance_all.tsv`
- `supplementary_table_6b_model_concordance_summary.tsv`
- `supplementary_table_6c_model_concordance_top.tsv`

Current source files:

- `results/mapt_model_concordance.tsv`
- `results/mapt_model_concordance_summary.tsv`
- `results/mapt_model_concordance_top.tsv`

Description: agreement and disagreement across ESM-1v, the transparent Tau
heuristic baseline, and AlphaMissense using the top-decile high-priority rule.
The table marks whether each variant is high-priority for each method and assigns
a concordance category.

Use in manuscript: supports the model concordance/discordance Results section and
Discussion.

## Notes for submission packaging

Large generated TSV files are not committed to GitHub. For journal submission,
copy the current source files to the proposed supplementary-table filenames and
upload them as supplementary data or deposit them in an external data repository.