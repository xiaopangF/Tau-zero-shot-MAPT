# Discussion Draft

## Principal findings

This study produced a complete zero-shot missense atlas for the 441-aa Tau-F /
hTau40 isoform. The atlas covers all 8379 possible single amino-acid substitutions
and combines ESM-1v ensemble scores with Tau-specific mechanistic annotations,
positive-control calibration, ClinVar coordinate QC, supplementary AlphaMissense QC,
and ranked VUS prioritization.

The main finding is not that the model is clinically validated. The main finding is
that a reproducible, label-free workflow can generate a full Tau-F prioritization map
while explicitly controlling for MAPT isoform-coordinate mismatch.

## Biological interpretation

The ESM-1v ensemble assigned higher average scores to the microtubule repeat regions
than to the N-terminal projection domain. This is biologically plausible because the
repeat regions are central to Tau microtubule binding and aggregation-related biology.
However, the five-control calibration showed that known MAPT missense controls were
not concentrated in the ESM-1v top decile. The regional signal should therefore be
interpreted as a prioritization pattern, not proof that every high-scoring substitution
is pathogenic.

## Calibration against established MAPT pathogenic controls

Five established MAPT missense controls, G272V, P301L, V337M, R406W, and N279K, were
all present in the atlas, but none was placed in the ESM-1v top decile. G272V was the
highest-ranked control at 2097/8379, corresponding to the top 25.03% of scores; the
remaining controls fell between the top 36.71% and top 46.19%. The present ESM-1v
score should therefore not be described as a calibrated MAPT pathogenicity probability.
It is better interpreted as a sequence-compatibility prioritization signal whose
usefulness depends on independent Tau biology and variant-level evidence.

## Value of strict coordinate QC

Strict ClinVar and AlphaMissense matching showed that many external MAPT rows cannot
be safely mapped to the 441-aa Tau-F reference. AlphaMissense covered only 877 of
8379 atlas variants, so it is retained as a supplementary reference rather than a
primary full-atlas comparator. This does not mean that ClinVar or AlphaMissense are
low-quality resources; it means that MAPT isoform biology makes naive merging risky.

## Interpretation of the ClinVar benchmark

The strict binary ClinVar benchmark contains only three P/LP versus B/LB examples.
Therefore, AUROC and AUPRC values should be treated as smoke checks rather than
clinical-performance estimates. Current VUS rankings are also hypothesis-generating,
not clinical interpretations.

## VUS prioritization and a falsifiable hypothesis

G440E is a concrete hypothesis-generating case. Although it lies in the C-terminal
tail rather than a canonical microtubule-binding repeat, its high score could reflect
sequence constraints related to Tau conformational dynamics, membrane association, or
phosphorylation accessibility. This should be tested in cellular or biochemical assays
rather than treated as a conclusion.

## Limitations and future work

ESM-1v does not directly model splicing, isoform expression, post-translational
modification state, neuronal cell context, aggregation kinetics, or patient-level
genetic background. The negative positive-control calibration indicates that a generic
protein language model may not rank every known Tau disease mutation highly.
AlphaMissense coverage is incomplete, and all current VUS candidates require
independent evidence and experimental validation.

The next analyses should add population-frequency filtering, conservation and disorder
features, mechanistic regression, and functional assays. In particular, we will test
whether ESM scores are associated with conservation, hydrophobicity, disorder, charge
change, or aggregation-related features.

## Conclusion

This work provides a reproducible zero-shot Tau-F missense atlas with strict coordinate
QC and interpretable prioritization layers. Its most defensible use is as a transparent
starting point for MAPT/Tau variant review, hypothesis generation, and follow-up
prioritization, not as a standalone clinical diagnostic model.

## Plain-language summary

For a beginner, the discussion means:

1. We made a complete Tau mutation map.
2. The AI model thinks Tau repeat regions are especially sensitive, but known disease mutations were not ranked high enough.
3. So we cannot claim this is a clinical diagnostic model.
4. AlphaMissense covers only a small part of our Tau map, so it is supplementary.
5. The next scientific step is to explain the scores using conservation, disorder, hydrophobicity, and population frequency.