# Results Draft

## A complete Tau-F missense atlas

We generated a complete single amino-acid substitution atlas for the 441-residue
2N4R Tau-F / hTau40 reference sequence. The atlas contains 8379 missense
variants, corresponding to 19 possible substitutions at each reference position.
Each variant was annotated with Tau region, repeat-domain membership, aggregation
motif proximity, post-translational modification proximity, charge-change class,
special-residue changes, and proximity to known pathogenic MAPT hotspots.

This design keeps the scoring layer label-free. Clinical annotations are added
only after the full variant list and zero-shot model scores have been generated.

## Strict clinical-label coordinate QC limits direct benchmarking

ClinVar import was performed with strict 441-aa Tau-F coordinate checking. A
ClinVar missense row was accepted only when the reported reference amino acid
matched the Tau-F reference residue at the same position. This retained 33
annotated variants and rejected 1119 MAPT rows that did not cleanly map to the
Tau-F atlas.

Accepted ClinVar labels were sparse: 1 pathogenic/likely pathogenic variant, 2
benign/likely benign variants, 26 variants of uncertain significance, and 4
variants with conflicting classifications. Therefore, the direct P/LP versus B/LB
benchmark contains only 3 examples. All AUROC and AUPRC values from this strict
benchmark should be interpreted as pipeline sanity checks rather than estimates
of clinical performance.

## ESM-1v ensemble scores highlight Tau repeat regions

We scored all 8379 variants with a five-checkpoint ESM-1v ensemble and used the
mean negative log-likelihood ratio as the pathogenicity-prioritization score. At
the domain level, the highest average scores were observed in the microtubule
repeat regions, while the N-terminal projection domain had the lowest average
score.

Mean ESM-1v ensemble pathogenic scores by region were:

| Region | Mean score | Global top 10 percent fraction |
|---|---:|---:|
| microtubule_repeat_R4 | 13.22 | 0.262 |
| microtubule_repeat_R3 | 13.21 | 0.255 |
| microtubule_repeat_R2_exon10 | 12.96 | 0.248 |
| microtubule_repeat_R1 | 12.64 | 0.211 |
| C_terminal_tail | 12.03 | 0.114 |
| proline_rich_region | 9.83 | 0.057 |
| N_terminal_projection | 3.08 | 0.000 |

This domain pattern is consistent with the biological importance of the repeat
region, but it should be presented as a prioritization signal rather than a
clinical validation result.

## Comparison with a transparent Tau heuristic baseline

To test whether the language-model scores merely reproduce obvious Tau-specific
rules, we built a transparent heuristic baseline from region, motif, hotspot,
PTM, charge-change, and special-residue features. Across all 8379 variants, the
ESM-1v ensemble and the Tau heuristic baseline had Pearson r=0.58, indicating
moderate overlap.

This result supports using the heuristic baseline as an interpretability control.
The two methods are related but not identical, which creates an opportunity to
study variants where ESM-1v and explicit Tau-domain rules disagree.

## ClinVar VUS prioritization

Because most strictly mapped ClinVar variants were VUS, the most useful clinical
analysis in the current data is prioritization rather than binary classification.
The top ESM-1v-ranked ClinVar VUS candidates were:

| Rank | Variant | Score | Region | Mechanistic annotation |
|---:|---|---:|---|---|
| 1 | G440E | 10.37 | C_terminal_tail | charge_gain_or_loss; G_loss |
| 2 | G333A | 8.73 | microtubule_repeat_R3 | G_loss |
| 3 | K290Q | 8.48 | microtubule_repeat_R2_exon10 | charge_gain_or_loss |
| 4 | T377A | 7.02 | C_terminal_tail | region signal |
| 5 | D34Y | 4.12 | N_terminal_projection | charge_gain_or_loss |

These variants are candidates for literature review, database cross-checking,
and possible experimental prioritization. They should not be described as
reclassified pathogenic variants without independent evidence.

## AlphaMissense provides an external QC and coverage comparison

We imported the public AlphaMissense hg38 table as an external reference model
using strict Tau-F reference-residue QC. AlphaMissense scores were available for
877 of 8379 Tau-F missense variants, covering 149 positions. Among scored
variants, 798 were AlphaMissense likely benign, 67 ambiguous, and 12 likely
pathogenic.

The importer rejected 3932 MAPT AlphaMissense rows after filtering for P10636:
2048 were outside the 441-aa Tau-F reference range, and 1884 had a reference
amino acid mismatch. This confirms that AlphaMissense is useful as an external
reference and coordinate-QC example, but its direct coverage of the 441-aa Tau-F
atlas is limited.

In the strict ClinVar binary benchmark, AlphaMissense scored only 2 examples.
Therefore, its apparent AUROC/AUPRC in the local model-comparison table should
not be reported as a performance claim.

## Current interpretation

The current evidence supports a manuscript centered on atlas generation,
coordinate QC, domain-level zero-shot behavior, and VUS prioritization. The
strict clinical benchmark is too small to support strong claims about clinical
classification accuracy. The most defensible claim is that the project provides a
reproducible Tau-F missense prioritization atlas with transparent QC and multiple
external/contextual interpretation layers.