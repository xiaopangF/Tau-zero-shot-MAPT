# Figure Legends Draft

## Figure 1. MAPT/Tau zero-shot atlas workflow

Workflow for constructing the Tau-F missense variant atlas. The 441-aa Tau-F /
hTau40 / 2N4R reference sequence was used to enumerate all 8379 possible single
amino-acid substitutions. Variants were scored with pretrained ESM models without
using clinical labels for model training or tuning. ClinVar and AlphaMissense
were added afterward as QC, evaluation, and interpretation layers. Strict
reference-residue matching was used to prevent mixing Tau-F coordinates with
other MAPT isoform coordinates. Final outputs include domain-level score
summaries and ranked ClinVar VUS candidates.

Plain-language meaning: this figure shows the recipe for the whole project.

## Figure 2. ESM-1v missense score atlas across Tau-F

Heatmap and position-level view of ESM-1v ensemble scores for all 8379 Tau-F
missense variants. Each position corresponds to one amino-acid site in Tau-F, and
each possible mutant amino acid receives a zero-shot score. Higher
`pathogenic_score_mean` indicates that the mutant amino acid is less compatible
with the sequence patterns learned by the pretrained model. The figure provides a
complete visual atlas of model-prioritized substitutions across the Tau-F
protein.

Plain-language meaning: this figure is the mutation risk map.

## Figure 3. Domain-level ESM-1v score patterns

Mean ESM-1v ensemble pathogenic scores were summarized by Tau-F region. The
microtubule repeat regions R1-R4 showed higher average scores than the
N-terminal projection domain. The plot also shows the fraction of variants in
each region that fall into the global top 10 percent of ESM-1v scores. This
regional pattern suggests that the zero-shot model assigns greater sensitivity to
Tau repeat regions, which are central to microtubule binding and Tau aggregation
biology.

Plain-language meaning: this figure asks which parts of Tau look most sensitive
to mutation.

## Figure 4. Model comparison and heuristic baseline control

Model scores were compared against the strict ClinVar binary benchmark using only
variants labeled pathogenic/likely pathogenic or benign/likely benign after
Tau-F coordinate QC. The benchmark contains very few examples and should be
interpreted as a smoke check rather than a clinical performance estimate. The
figure also includes comparison with a transparent Tau heuristic baseline and a
scatter plot comparing ESM-1v ensemble scores with heuristic scores across all
8379 variants.

Plain-language meaning: this figure checks whether the model behaves reasonably,
but it does not prove clinical accuracy.

## Figure 5. Prioritized ClinVar VUS candidates across Tau-F

Lollipop plot of the top-ranked ClinVar variants of uncertain significance
(VUS), ordered by ESM-1v ensemble `pathogenic_score_mean`. Each point marks the
Tau-F position of a VUS candidate, and the height of the point indicates its
zero-shot score. Tau-F regions are shown along the protein axis. Candidate
variants are labeled by protein-change identifier and should be interpreted as
priorities for follow-up review, not as clinical reclassifications.

Plain-language meaning: this figure shows the uncertain variants that our model
thinks deserve the most attention.

## Supplementary Figure or Table Notes

### Supplementary Table 1. Full Tau-F missense atlas

Complete list of all 8379 single amino-acid substitutions in Tau-F, including
variant identifiers, protein changes, positions, reference amino acids, mutant
amino acids, and Tau-specific annotations.

### Supplementary Table 2. ESM-1v ensemble scores

Zero-shot scores for all 8379 Tau-F missense variants. The table includes scores
from the five-checkpoint ESM-1v ensemble, including the mean pathogenic score and
between-model standard deviation.

### Supplementary Table 3. ClinVar coordinate QC

Accepted and rejected ClinVar MAPT rows under strict 441-aa Tau-F coordinate
matching. Rejected rows document whether the issue was non-parseable missense
annotation, an out-of-range position, or reference-residue mismatch.

### Supplementary Table 4. AlphaMissense coordinate QC

Accepted and rejected AlphaMissense MAPT rows under strict Tau-F coordinate
matching. This table documents AlphaMissense coverage of the Tau-F atlas and
explains why AlphaMissense is used as an external reference rather than a full
main-model benchmark.

### Supplementary Table 5. Top ClinVar VUS candidates

Ranked list of ClinVar VUS candidates prioritized by ESM-1v ensemble score,
including region, mechanism annotation, review status, and score.