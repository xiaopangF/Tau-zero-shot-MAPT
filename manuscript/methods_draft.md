# Methods Draft

## Study design

This study builds a label-free missense variant atlas for MAPT/Tau. In simple
terms, we first listed every possible single amino-acid change in Tau-F, then
used pretrained protein language models to score those changes without using
clinical labels for training or tuning. Clinical databases were added only after
scoring, as evaluation and interpretation layers.

The intended use of the atlas is variant prioritization. It is not designed as a
standalone clinical diagnostic test.

## Tau reference sequence and coordinate system

All analyses used the 441-amino-acid adult CNS Tau-F / hTau40 / 2N4R isoform as
the reference sequence. This coordinate system corresponds to the commonly used
441-aa Tau sequence and is tracked in the project as the Tau-F / P10636-8 atlas
coordinate system.

This choice is important because MAPT has multiple isoforms. Some external
resources use different transcript or protein coordinates. To avoid mixing
coordinate systems, all clinical and external predictor rows were checked against
the 441-aa Tau-F reference residue before being merged into the atlas.

## Missense variant enumeration

For each of the 441 Tau-F positions, all 19 possible non-reference amino-acid
substitutions were generated. This produced:

`441 positions x 19 substitutions = 8379 missense variants`

Each variant was assigned a compact identifier such as `P301L`, meaning the
reference amino acid proline at position 301 is changed to leucine.

## Tau-specific variant annotation

Each variant was annotated with interpretable Tau biology features. These
included:

- Tau region: N-terminal projection domain, proline-rich region, repeat regions
  R1-R4, or C-terminal tail.
- Aggregation motifs, including PHF6-star/VQIINK and PHF6/VQIVYK.
- Proximity to selected Tau post-translational modification sites.
- Proximity to known pathogenic MAPT hotspot positions.
- Charge-changing substitutions.
- Changes involving special residues such as proline, glycine, or cysteine.

These annotations were not used to train the ESM model. They were used for
interpretation, prioritization, and the transparent heuristic baseline.

## Zero-shot ESM scoring

Variants were scored with pretrained ESM protein language models using masked
marginal scoring. For each Tau-F position, the model estimated the log probability
of the reference amino acid and the mutant amino acid in the same sequence
context.

The raw ESM log-likelihood ratio was calculated as:

`esm_llr = mutant_log_probability - wildtype_log_probability`

Lower `esm_llr` means the mutant amino acid is less compatible with the sequence
patterns learned by the model. For easier prioritization, the project uses:

`pathogenic_score = -esm_llr`

Therefore, higher `pathogenic_score` means the variant is more suspicious under
the zero-shot model.

The manuscript-scale run used five ESM-1v checkpoints. Scores from the five
models were averaged to produce `pathogenic_score_mean`; between-model variation
was stored as `pathogenic_score_std`.

## ClinVar import and strict coordinate QC

ClinVar annotations were imported after zero-shot scoring. A ClinVar MAPT row was
accepted only if it described a parseable missense variant and the reported
reference amino acid matched the 441-aa Tau-F reference sequence at that position.

Rows were rejected if:

- no parseable missense protein change could be extracted;
- the reported position was outside the 441-aa Tau-F reference range;
- the reported reference amino acid did not match Tau-F.

This strict QC prevents variants from other MAPT isoforms or coordinate systems
from being incorrectly treated as Tau-F variants.

ClinVar labels were grouped into four interpretation classes:

- `P_LP`: pathogenic or likely pathogenic;
- `B_LB`: benign or likely benign;
- `VUS`: uncertain significance;
- `conflicting`: conflicting classifications of pathogenicity.

Only `P_LP` and `B_LB` were used for binary AUROC/AUPRC smoke checks.

## AlphaMissense import and external predictor QC

The public AlphaMissense hg38 score table was imported as an external reference
predictor. Rows were filtered to MAPT/UniProt `P10636` and then subjected to the
same strict Tau-F reference-residue QC.

AlphaMissense rows were accepted only when the parsed protein substitution
matched the 441-aa Tau-F reference residue at the mapped position. Rows outside
the 1-441 range or with reference-residue mismatches were written to a rejected
QC table.

Because AlphaMissense coverage under strict Tau-F coordinates was incomplete,
AlphaMissense was treated as an external comparison and QC layer, not as the main
benchmark for manuscript performance claims.

## Transparent Tau heuristic baseline

A simple rule-based Tau heuristic baseline was created to test whether ESM scores
only reflected obvious Tau-domain features. The heuristic assigns points for
features such as repeat-domain location, aggregation motif proximity, known
hotspot proximity, PTM proximity, charge changes, and special-residue changes.

This baseline is not a learned model. Its purpose is interpretability: if ESM and
the heuristic agree, the ESM score may be partly explained by known Tau biology;
if they disagree, the variant may be useful for follow-up analysis.

## Model evaluation and prioritization

Binary evaluation used only strict ClinVar `P_LP` and `B_LB` variants. Metrics
included AUROC, AUPRC, and top-percentile enrichment. Because the strict binary
ClinVar set contained only a small number of examples, these metrics were treated
as pipeline checks rather than clinical performance estimates.

For VUS prioritization, ClinVar VUS variants were ranked by ESM-1v ensemble
`pathogenic_score_mean`. Ranked variants were reported with their Tau region,
mechanistic annotations, ClinVar review status, and score.

## Domain-level summary

To summarize regional behavior, ESM-1v ensemble scores were grouped by Tau-F
region. For each region, the mean, median, interquartile range, maximum score,
and fraction of variants in the global top score percentiles were calculated.

This analysis asks whether the zero-shot model assigns systematically different
sensitivity to different Tau regions.

## Software and reproducibility

The workflow is implemented as a reproducible Python command-line package. Core
steps include variant enumeration, ClinVar import, AlphaMissense import, ESM
scoring, ensemble averaging, benchmarking, prioritization, summary-table
generation, and manuscript figure generation.

Generated large data files and result tables are not committed to GitHub. The
repository stores the code and instructions needed to regenerate them.