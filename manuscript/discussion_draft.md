# Discussion Draft

## Principal findings

This study produced a complete zero-shot missense atlas for the 441-aa Tau-F /
hTau40 isoform. The atlas covers all 8379 possible single amino-acid
substitutions and combines ESM-1v ensemble scores with Tau-specific mechanistic
annotations, ClinVar coordinate QC, AlphaMissense coordinate QC, and ranked VUS
prioritization.

The main finding is not that the model is clinically validated. The main finding
is that a reproducible, label-free workflow can generate a full Tau-F
prioritization map while explicitly controlling for MAPT isoform-coordinate
mismatch. This distinction is important because MAPT has multiple isoforms and
external variant resources do not always align cleanly to the 441-aa Tau-F
reference sequence.

## Biological interpretation

The ESM-1v ensemble assigned higher average pathogenicity-prioritization scores
to the microtubule repeat regions than to the N-terminal projection domain. This
pattern is biologically plausible because the repeat regions are central to Tau
microtubule binding and aggregation-related biology. The result suggests that the
protein language model captures broad regional constraint in Tau-F, even though
it was not trained specifically on MAPT clinical labels.

This domain-level signal should be interpreted as a prioritization pattern. It
does not prove that every high-scoring repeat-region substitution is pathogenic.
Rather, it identifies regions and substitutions where the model sees stronger
sequence incompatibility and where follow-up biological review is most justified.

## Value of strict coordinate QC

A major practical lesson from this project is that coordinate QC is not a minor
technical detail. Strict ClinVar matching accepted only a small number of
annotated Tau-F missense variants and rejected many MAPT rows that could not be
safely mapped to the 441-aa reference. AlphaMissense showed a similar issue: only
877 of 8379 Tau-F missense variants received accepted AlphaMissense scores after
strict matching, while thousands of MAPT rows were rejected because they were
outside the Tau-F range or had reference-residue mismatches.

This does not mean ClinVar or AlphaMissense are low-quality resources. It means
that MAPT isoform biology makes naive merging risky. For a Tau-specific atlas,
reference-residue checking is necessary to avoid mixing data from different
coordinate systems.

## Interpretation of the ClinVar benchmark

The strict binary ClinVar benchmark contains only three P/LP versus B/LB examples.
Therefore, AUROC, AUPRC, and enrichment values should be treated as smoke checks.
They can show that the evaluation code runs and that score direction is handled
consistently, but they cannot support a strong clinical-performance claim.

For this reason, the manuscript should avoid language such as "accurately
predicts pathogenicity" or "clinically validates MAPT variants." Safer language
is "prioritizes variants," "generates a zero-shot atlas," and "provides a
reproducible framework for follow-up."

## VUS prioritization as the most useful current clinical layer

Most strictly mapped ClinVar variants were VUS. This makes VUS prioritization a
more appropriate current use case than binary classification. The ranked VUS list
highlights variants such as G440E, G333A, K290Q, T377A, and D34Y for further
review. These variants should be treated as candidates for literature review,
population-frequency checking, segregation analysis, functional assays, or expert
curation.

Importantly, prioritization is not reclassification. A high ESM-1v score means a
variant is worth attention, not that it is clinically pathogenic.

## Relationship to the heuristic baseline

The transparent Tau heuristic baseline had moderate correlation with ESM-1v
scores across all variants. This indicates that the language model and explicit
Tau-domain rules overlap but are not identical. This is useful in two ways. First,
it helps interpret high ESM scores that occur in known sensitive regions or
motifs. Second, it identifies discordant variants where the model assigns high
priority even though simple Tau rules do not, or vice versa.

A future version of the manuscript should include a dedicated concordance and
discordance table. This table would help reviewers see where the protein language
model adds information beyond obvious handcrafted features.

## Limitations

This project has several important limitations.

First, the clinical benchmark is very small after strict Tau-F coordinate QC.
This prevents strong claims about clinical accuracy.

Second, ESM-1v is a sequence model. It does not directly model all disease
mechanisms, including splicing effects, isoform expression, post-translational
modification state, neuronal cell context, aggregation kinetics, or patient-level
genetic background.

Third, AlphaMissense coverage of the 441-aa Tau-F atlas is incomplete under strict
reference matching. AlphaMissense is therefore useful as an external reference
and QC example, but not as a full head-to-head benchmark for all Tau-F variants.

Fourth, the current VUS rankings have not been experimentally validated. They are
hypothesis-generating candidates, not clinical interpretations.

## Future work

The next computational step is to add a concordance and discordance analysis
across ESM-1v, the Tau heuristic baseline, and AlphaMissense where coverage
exists. This would identify variants consistently prioritized by multiple methods
and variants where methods disagree.

The next biological step is to review the top VUS candidates against additional
evidence, including population frequency, literature reports, segregation data,
functional studies, and domain-specific Tau biology. If experimental work is
possible, high-priority variants could be tested in assays related to Tau
microtubule binding, aggregation propensity, phosphorylation, or cellular toxicity.

## Conclusion

This work provides a reproducible zero-shot Tau-F missense atlas with strict
coordinate QC and interpretable prioritization layers. The strongest claim is not
clinical diagnosis. The strongest claim is that the atlas gives researchers a
transparent starting point for MAPT/Tau variant review, hypothesis generation,
and follow-up prioritization.

## Plain-language summary

For a beginner, the discussion means:

1. We made a complete Tau mutation map.
2. The AI model thinks Tau repeat regions are especially sensitive.
3. The clinical database has too few clean examples to prove accuracy.
4. AlphaMissense is useful, but its coordinates do not fully match our Tau map.
5. The safest use is to rank uncertain variants for follow-up, not to diagnose
   patients.