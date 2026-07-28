# Manuscript Outline

## Title

Zero-shot prioritization and mechanistic annotation of MAPT/Tau missense variants

## Abstract

1. MAPT/Tau variants cause inherited tauopathies, but many missense variants
   remain difficult to interpret.
2. We generated a full missense atlas for canonical 441-aa Tau and scored each
   substitution in 2N4R Tau-F/P10636-8 using pretrained protein sequence models without clinical label
   training.
3. We benchmarked label-free scores against ClinVar/UniProt annotations and
   annotated high-priority variants by Tau region, microtubule repeats,
   aggregation motifs, PTM proximity, charge change, and pathogenic hotspot
   proximity.
4. The atlas prioritizes known pathogenic variants and proposes ranked VUS
   candidates for follow-up.

## Main sections

- Introduction: MAPT genetics, Tau biology, limitations of generic predictors.
- Methods: reference sequence, variant enumeration, zero-shot scoring, clinical
  label handling, Tau-specific annotations, metrics.
- Results: full atlas, benchmark performance, domain-level behavior, model
  concordance/discordance, VUS prioritization.
- Discussion: strengths, Tau-specific failure modes, clinical limitations, and
  experimental follow-up.

