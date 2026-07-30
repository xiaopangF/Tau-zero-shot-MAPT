# Analysis Notes

## Current smoke workflow status

These notes describe the local lightweight ESM-2 smoke run, not the final
manuscript-scale ESM-1v result.

## ClinVar coordinate QC

Strict 441-aa Tau-F coordinate QC accepted 33 annotated missense variants and
rejected 1119 MAPT ClinVar rows.

Accepted labels:

- B/LB: 2
- P/LP: 1
- VUS: 26
- conflicting: 4

Rejected reasons:

- no_parseable_missense: 717
- outside_reference_range: 288
- reference_wt_mismatch: 114

Interpretation: direct ClinVar benchmarking is small under strict Tau-F
coordinates, so manuscript claims should emphasize atlas generation, QC,
mechanistic annotation, and VUS prioritization rather than AUROC alone.

## Lightweight ESM-2 domain signal

Using `esm2_t6_8M_UR50D` and `pathogenic_score`, the microtubule repeat regions
show higher average scores than the projection domains. This is a smoke-model
observation and should be re-evaluated with ESM-1v ensemble before manuscript use.

Top mean scores by region in the smoke run:

- R2 / exon 10: 3.43
- R3: 2.93
- R4: 2.87
- R1: 2.58
- proline-rich region: 1.58
- N-terminal projection: 0.96
- C-terminal tail: 0.88

## Heuristic baseline smoke result

The transparent `tau_heuristic_v1` baseline has been added to control for obvious
Tau-specific rules. On the tiny strict ClinVar benchmark it gives AUROC 0.75 and
AUPRC 1.0, while the lightweight ESM-2 smoke model performs poorly. Because the
benchmark has only three P/LP vs B/LB examples, this should be treated as a
sanity check and reviewer-facing motivation for including a simple baseline, not
as a performance conclusion.

The important manuscript question becomes: does the ESM-1v ensemble prioritize
known/VUS MAPT variants beyond this transparent Tau heuristic baseline?

## ESM-1v ensemble local run

The full ESM-1v five-checkpoint ensemble has been run locally on the RTX 5060 GPU.
Generated local files include:

- `results/scores/mapt_esm1v_1.tsv`
- `results/scores/mapt_esm1v_2.tsv`
- `results/scores/mapt_esm1v_3.tsv`
- `results/scores/mapt_esm1v_4.tsv`
- `results/scores/mapt_esm1v_5.tsv`
- `results/scores/mapt_esm1v_ensemble.tsv`
- `results/mapt_esm1v_ensemble_metrics.tsv`
- `results/mapt_esm1v_ensemble_domain_summary.tsv`
- `results/mapt_esm1v_top50_vus_priority.tsv`
- `results/figures/esm1v_ensemble/`

Strict ClinVar binary benchmark remains too small for performance claims:

- n_examples: 3
- n_pathogenic: 1
- n_benign: 2
- AUROC: 0.5
- AUPRC: 0.5

Domain-level ESM-1v ensemble signal is stronger in microtubule repeat regions
than the N-terminal projection domain. Mean `pathogenic_score_mean` by region:

- R4: 13.22
- R3: 13.21
- R2 / exon 10: 12.96
- R1: 12.64
- C-terminal tail: 12.03
- proline-rich region: 9.83
- N-terminal projection: 3.08

Interpretation: the project now has manuscript-scale zero-shot scores, but
clinical-label evaluation is still limited by strict MAPT/Tau-F coordinate QC.
Next manuscript work should emphasize atlas/domain/VUS analyses and add external
baselines such as AlphaMissense.

## Manuscript assets generated

Local manuscript assets were generated in `results/manuscript_assets/`:

- `figure_domain_summary.png`
- `figure_model_comparison.png`
- `figure_esm1v_vs_heuristic.png`
- `results_summary.md`
- `top_vus_candidates.md`

The ESM-1v ensemble and Tau heuristic have Pearson r=0.58 across all 8379
variants, suggesting moderate overlap between learned zero-shot signal and the
transparent Tau-rule baseline. This supports including both in the manuscript.

## AlphaMissense external baseline QC

The public `AlphaMissense_hg38.tsv.gz` table was imported with `--uniprot-id
P10636` and strict 441-aa Tau-F reference-residue QC.

Coverage in the current local run:

- total Tau-F missense variants: 8379
- AlphaMissense-scored variants: 877
- AlphaMissense-scored positions: 149
- likely benign: 798
- ambiguous: 67
- likely pathogenic: 12

Rejected AlphaMissense MAPT rows:

- outside_reference_range: 2048
- reference_wt_mismatch: 1884

Interpretation: AlphaMissense is useful as an external reference and QC example,
but its direct coverage of the 441-aa Tau-F atlas is limited. In the strict
ClinVar binary benchmark it scores only 2 variants, so the AlphaMissense AUROC
and AUPRC from `results/mapt_model_comparison.tsv` should not be used as a
manuscript performance claim. The stronger manuscript use is to report coverage,
coordinate mismatch, and selected concordance/discordance examples.

## Positive-control calibration

The five-control calibration was run with `scripts/run_positive_control_calibration.ps1`.
All five established MAPT missense controls were found in the 8379-variant atlas, but
none entered the ESM-1v top 10%. The best-ranked control was G272V at rank 2097/8379
(top 25.03%). This is now treated as a major calibration limitation rather than a
hidden result.