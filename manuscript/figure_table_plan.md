# Figure and Table Plan

## Main figures

### Figure 1: Study workflow and coordinate QC

Purpose: show the pipeline from Tau-F sequence selection to all missense variants,
zero-shot scoring, ClinVar/AlphaMissense QC, and VUS prioritization.

Current source files:

- `data/processed/mapt_all_missense_variants.tsv`
- `data/processed/mapt_clinvar_benchmark.tsv`
- `data/processed/mapt_clinvar_rejected.tsv`
- `data/processed/mapt_alphamissense.tsv`
- `data/processed/mapt_alphamissense_rejected.tsv`

Status: generated locally as `results/manuscript_assets/figure_workflow_schematic.png`.

### Figure 2: ESM-1v Tau-F missense atlas

Purpose: show zero-shot pathogenic-score distribution across Tau-F positions and
mutant amino acids.

Current source files:

- `results/scores/mapt_esm1v_ensemble.tsv`
- `results/figures/esm1v_ensemble/missense_heatmap.png`
- `results/figures/esm1v_ensemble/score_by_position.png`

Status: generated locally; likely needs publication styling.

### Figure 3: Domain-level ESM-1v signal

Purpose: show that repeat regions have higher mean ESM-1v scores than the
N-terminal projection domain.

Current source files:

- `results/mapt_esm1v_ensemble_domain_summary.tsv`
- `results/manuscript_assets/figure_domain_summary.png`

Status: generated locally.

### Figure 4: Model comparison and baseline control

Purpose: show the strict ClinVar benchmark as a smoke check and compare ESM-1v
with the transparent Tau heuristic baseline.

Current source files:

- `results/mapt_model_comparison.tsv`
- `results/manuscript_assets/figure_model_comparison.png`
- `results/manuscript_assets/figure_esm1v_vs_heuristic.png`

Important note: benchmark n is too small for clinical performance claims.

### Figure 5: Prioritized ClinVar VUS candidates

Purpose: present top VUS candidates and their mechanistic annotations.

Current source files:

- `results/mapt_esm1v_top50_vus_priority.tsv`
- `results/manuscript_assets/top_vus_candidates.md`

Status: generated locally as `results/manuscript_assets/figure_vus_lollipop.png`; table also exists.

## Supplementary tables

### Supplementary Table 1: Full Tau-F missense atlas

Source: `data/processed/mapt_all_missense_variants.tsv`

### Supplementary Table 2: ESM-1v ensemble scores

Source: `results/scores/mapt_esm1v_ensemble.tsv`

### Supplementary Table 3: ClinVar accepted and rejected coordinate QC

Sources:

- `data/processed/mapt_clinvar_benchmark.tsv`
- `data/processed/mapt_clinvar_rejected.tsv`
- `results/mapt_clinvar_qc_summary.tsv`

### Supplementary Table 4: AlphaMissense coverage and rejected rows

Sources:

- `data/processed/mapt_alphamissense.tsv`
- `data/processed/mapt_alphamissense_rejected.tsv`
- `results/mapt_alphamissense_qc_summary.tsv`

### Supplementary Table 5: Top VUS priority list

Source: `results/mapt_esm1v_top50_vus_priority.tsv`

### Supplementary Table 6: Model concordance and discordance

Sources:

- `results/mapt_model_concordance.tsv`
- `results/mapt_model_concordance_summary.tsv`
- `results/mapt_model_concordance_top.tsv`

## Immediate gaps

- Add a publication-style workflow schematic.
- Add a VUS lollipop or position plot.
- Add a concordance/discordance table for ESM-1v, heuristic, and AlphaMissense
  where AlphaMissense coverage exists.
- Add Methods text for coordinate systems and strict reference-residue matching.