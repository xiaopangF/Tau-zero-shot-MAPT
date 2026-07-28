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
