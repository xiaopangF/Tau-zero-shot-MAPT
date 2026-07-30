# MAPT/Tau Zero-Shot Variant Atlas

This project builds a reproducible, label-free prioritization atlas for MAPT/Tau
missense variants. The v1 scope is a computational bioinformatics project:
enumerate all single amino-acid substitutions in canonical 441-aa human Tau,
score them with zero-shot protein language models, benchmark against clinical
annotations, and add Tau-specific mechanistic annotations.

## Project outputs

- `data/processed/mapt_all_missense_variants.tsv`: all 441 x 19 MAPT missense variants.
- `results/scores/mapt_esm_scores.tsv`: zero-shot ESM scores.
- `data/processed/mapt_clinvar_benchmark.tsv`: variants merged with ClinVar labels.
- `results/mapt_metrics.tsv`: AUROC/AUPRC and enrichment metrics.
- `results/figures/`: publication-oriented figures.

## Setup

Use Python 3.10-3.12 for the full GPU workflow because PyTorch and fair-esm
support often lags newer Python versions.

```powershell
conda env create -f environment.yml
conda activate mapt-zero-shot
pip install -e ".[analysis,esm]"
```

The core enumeration and TSV utilities use only the Python standard library.

## Minimal workflow

Generate the full missense atlas. Coordinates are for 2N4R Tau-F / hTau40 (`P10636-8`), not UniProt canonical Big-tau (`P10636-1`).

```powershell
python -m mapt_zero_shot.cli enumerate `
  --out data/processed/mapt_all_missense_variants.tsv
```

Merge ClinVar annotations after downloading `variant_summary.txt.gz` from NCBI:

```powershell
python -m mapt_zero_shot.cli import-clinvar `
  --variants data/processed/mapt_all_missense_variants.tsv `
  --clinvar data/raw/variant_summary.txt.gz `
  --out data/processed/mapt_clinvar_benchmark.tsv
```

Score all missense variants with an ESM model:

```powershell
python -m mapt_zero_shot.cli score-esm `
  --model esm1v_t33_650M_UR90S_1 `
  --out results/scores/mapt_esm1v_scores.tsv
```

Evaluate scores against ClinVar labels:

```powershell
python -m mapt_zero_shot.cli evaluate `
  --scores results/scores/mapt_esm1v_scores.tsv `
  --labels data/processed/mapt_clinvar_benchmark.tsv `
  --score-column pathogenic_score `
  --out results/mapt_metrics.tsv
```

Create figures:

```powershell
python -m mapt_zero_shot.cli figures `
  --scores results/scores/mapt_esm1v_scores.tsv `
  --labels data/processed/mapt_clinvar_benchmark.tsv `
  --outdir results/figures
```

## Data policy

Clinical labels are never used to train or tune the zero-shot scores. ClinVar,
UniProt disease annotations, and population frequency data are evaluation and
interpretation layers only.

## Repository privacy note

Manuscript drafts, LaTeX sources, generated PDFs, and submission-package files are
kept local and are intentionally not tracked in GitHub. This public repository keeps
the reproducible analysis code and workflow documentation.
## Reproducible smoke workflow used in this workspace

The current workspace has already run the lightweight ESM-2 model
`esm2_t6_8M_UR50D` as an end-to-end smoke workflow. This is not the final model
for the paper, but it verifies enumeration, ESM scoring, ClinVar import,
prioritization, and figure generation.

```powershell
.\scripts\download_reference_data.ps1
.\scripts\run_smoke_workflow.ps1
```

Current generated files include:

- `data/processed/mapt_all_missense_variants.tsv` with 8379 variants.
- `data/raw/variant_summary.txt.gz` downloaded from NCBI ClinVar.
- `data/processed/mapt_clinvar_benchmark.tsv` with MAPT ClinVar annotations.
- `results/scores/mapt_esm2_t6_8M_scores.tsv` with 8379 ESM LLR and pathogenic_score values.
- `results/mapt_esm2_t6_8M_metrics.tsv` as a smoke benchmark table.
- `results/mapt_top50_vus_priority.tsv` for ClinVar VUS prioritization.
- `results/figures/missense_heatmap.png`.
- `results/figures/score_by_position.png`.
- `results/figures/clinvar_score_distribution.png`.

For manuscript-scale runs, replace the smoke model with ESM-1v or a larger ESM-2
checkpoint:

```powershell
python -m mapt_zero_shot.cli score-esm `
  --model esm1v_t33_650M_UR90S_1 `
  --out results/scores/mapt_esm1v_scores.tsv
```

ClinVar currently provides only a small number of MAPT missense P/LP and B/LB
benchmark examples after strict 441-aa coordinate matching. Treat the lightweight
metrics as a pipeline check, not a manuscript-level performance claim.

## Score direction

The ESM output contains two related score fields:

- `esm_llr = mutant_logp - wildtype_logp`; lower values mean the mutant residue is less compatible with the pretrained sequence model.
- `pathogenic_score = -esm_llr`; higher values are used for pathogenicity prioritization, ClinVar benchmarking, figures, and VUS ranking.

Use `pathogenic_score` for all main ranking analyses.

## ESM-1v ensemble

Run each ESM-1v checkpoint separately, then average them:

```powershell
python -m mapt_zero_shot.cli score-esm --model esm1v_t33_650M_UR90S_1 --out results/scores/mapt_esm1v_1.tsv
python -m mapt_zero_shot.cli score-esm --model esm1v_t33_650M_UR90S_2 --out results/scores/mapt_esm1v_2.tsv
python -m mapt_zero_shot.cli score-esm --model esm1v_t33_650M_UR90S_3 --out results/scores/mapt_esm1v_3.tsv
python -m mapt_zero_shot.cli score-esm --model esm1v_t33_650M_UR90S_4 --out results/scores/mapt_esm1v_4.tsv
python -m mapt_zero_shot.cli score-esm --model esm1v_t33_650M_UR90S_5 --out results/scores/mapt_esm1v_5.tsv
python -m mapt_zero_shot.cli ensemble `
  --scores results/scores/mapt_esm1v_1.tsv results/scores/mapt_esm1v_2.tsv results/scores/mapt_esm1v_3.tsv results/scores/mapt_esm1v_4.tsv results/scores/mapt_esm1v_5.tsv `
  --score-column pathogenic_score `
  --out results/scores/mapt_esm1v_ensemble.tsv
```

## ClinVar coordinate QC

`import-clinvar` keeps only protein substitutions whose wild-type residue matches
this atlas' 441-aa Tau-F reference sequence. Rejected MAPT ClinVar rows can be
written with `--rejected-out` and are categorized as:

- `no_parseable_missense`
- `outside_reference_range`
- `reference_wt_mismatch`

This protects the benchmark from mixing Big-tau or other MAPT isoform coordinates
into the 441-aa Tau-F atlas.

## Formal ESM-1v run script

Use the guarded script for the full ESM-1v ensemble:

```powershell
.\scripts\run_esm1v_ensemble.ps1
```

The script refuses to run ESM-1v on CPU by default. If you intentionally want a
long CPU run, pass `-AllowCpu`.

## Automated tests

The repository includes GitHub Actions CI in `.github/workflows/tests.yml`. It
installs the package with the lightweight test dependencies and runs:

```powershell
python -m pytest -q
```

## Heuristic baseline

A transparent `tau_heuristic_v1` baseline is available through:

```powershell
python -m mapt_zero_shot.cli score-heuristic `
  --annotations data/processed/mapt_all_missense_variants.tsv `
  --out results/scores/mapt_heuristic_scores.tsv
```

This is not a learned model. It assigns points for Tau region, aggregation motif,
known pathogenic hotspot proximity, PTM proximity, charge change, and Pro/Gly/Cys
changes. Its purpose is to test whether ESM scores add information beyond obvious
Tau-specific rules.

## Model comparison table

Use `compare-models` to put multiple methods in one benchmark table. Use relative
paths in `name:path:score_column` specs on Windows.

```powershell
python -m mapt_zero_shot.cli compare-models `
  --labels data/processed/mapt_clinvar_benchmark.tsv `
  --model esm2_t6_8M:results/scores/mapt_esm2_t6_8M_scores.tsv:pathogenic_score `
  --model tau_heuristic_v1:results/scores/mapt_heuristic_scores.tsv:heuristic_score `
  --out results/mapt_model_comparison.tsv
```

## CUDA environment used locally

The workstation has an NVIDIA GeForce RTX 5060 Laptop GPU. The default Python
3.14 environment had CPU-only PyTorch (`torch 2.13.0+cpu`), so CUDA was not
visible there. A project-local Python 3.12 virtual environment was created at
`.venv` with CUDA PyTorch:

```text
torch 2.11.0+cu128
cuda_available True
```

Use `.\scripts\check_cuda.ps1` to verify the active project environment.

For robust ESM-1v downloads, use:

```powershell
.\scripts\download_esm1v_checkpoints.ps1
```

This checks remote file sizes and prevents half-downloaded checkpoint caches from
being mistaken for valid weights.

## Publication-oriented summary assets

After running the ESM-1v ensemble, the workflow can generate summary plots and
ranked tables under `results/`. Generated result files are ignored by Git and should
be archived separately when needed.
## Zero-shot physicochemical model

The recommended physicochemical model is `score-physchem-zero-shot`:

```powershell
python -m mapt_zero_shot.cli score-physchem-zero-shot `
  --variants data/processed/mapt_all_missense_variants.tsv `
  --out results/scores/mapt_physchem_zero_shot_ranked.tsv `
  --summary-out results/mapt_physchem_zero_shot_summary.tsv
```

It computes Kyte-Doolittle hydrophobicity delta, Chou-Fasman beta-sheet
propensity delta, and net-charge delta. Each feature is standardized across the
complete unlabeled 8379-variant atlas, converted to an absolute perturbation
magnitude, and combined with fixed Tau priors for repeat-region membership,
PHF6/PHF6* aggregation motifs, and proximity to common Tau PTM sites. The model
does not use ClinVar labels or the five established pathogenic controls to choose
weights. Those controls appear only in the summary as an external calibration
check.

Use `fuse-zero-shot` to combine an ESM score table with the physicochemical table
after unlabeled z-score normalization. The default is ESM-dominant (`0.75`) with
a smaller physicochemical contribution (`0.25`):

```powershell
python -m mapt_zero_shot.cli fuse-zero-shot `
  --esm-scores results/scores/mapt_esm1v_ensemble.tsv `
  --physchem-scores results/scores/mapt_physchem_zero_shot_ranked.tsv `
  --out results/scores/mapt_zero_shot_composite_ranked.tsv `
  --summary-out results/mapt_zero_shot_composite_summary.tsv
```

`score-physchem-grid` is retained as a legacy reproducibility command. It searches
weights using the five controls and therefore is label-informed; it must not be
presented as the main zero-shot model in the paper.
## Scheme A knowledge-driven model

`score-knowledge-a` applies fixed, label-free rules to an existing
physicochemical feature table. It uses directed beta-sheet disruption in the Tau
repeat region, PTM-neighborhood charge/hydrophobicity disturbance, proximity to
the exon-10 splicing-sensitive center residue 285, and a small baseline
physicochemical z-score term.

```powershell
python -m mapt_zero_shot.cli score-knowledge-a `
  --physchem-features results/scores/mapt_physchem_zero_shot_ranked.tsv `
  --out results/scores/mapt_scheme_a_ranked.tsv `
  --gold-out results/mapt_scheme_a_gold_controls.tsv `
  --summary-out results/mapt_scheme_a_summary.tsv
```

The five established controls are not used to tune any Scheme A parameter. They
are reported only as a validation table with baseline physicochemical z-score,
directed beta score, PTM disturbance score, splicing score, final score, and rank.
## AlphaMissense external baseline

AlphaMissense can be imported as an external reference model after downloading the
public hg38 table:

```powershell
.\scripts\download_reference_data.ps1 -AlphaMissense
python -m mapt_zero_shot.cli import-alphamissense `
  --variants data/processed/mapt_all_missense_variants.tsv `
  --alphamissense data/external/AlphaMissense_hg38.tsv.gz `
  --uniprot-id P10636 `
  --out data/processed/mapt_alphamissense.tsv `
  --rejected-out data/processed/mapt_alphamissense_rejected.tsv
python -m mapt_zero_shot.cli summarize-alphamissense `
  --alpha data/processed/mapt_alphamissense.tsv `
  --rejected data/processed/mapt_alphamissense_rejected.tsv `
  --out results/mapt_alphamissense_qc_summary.tsv
```

The importer applies strict Tau-F reference-residue QC. This is important because
AlphaMissense MAPT entries are indexed by Ensembl/UniProt coordinates that do not
fully match the 441-aa Tau-F atlas. Rows whose wild-type residue does not match
Tau-F, or whose position is outside 1-441, are written to the rejected table.

AlphaMissense can also be included in model comparison:

```powershell
python -m mapt_zero_shot.cli compare-models `
  --labels data/processed/mapt_clinvar_benchmark.tsv `
  --model esm2_t6_8M:results/scores/mapt_esm2_t6_8M_scores.tsv:pathogenic_score `
  --model tau_heuristic_v1:results/scores/mapt_heuristic_scores.tsv:heuristic_score `
  --model esm1v_ensemble:results/scores/mapt_esm1v_ensemble.tsv:pathogenic_score_mean `
  --model alphamissense:data/processed/mapt_alphamissense.tsv:alphamissense_score `
  --out results/mapt_model_comparison.tsv
```

In the current local run, AlphaMissense covers 877 of 8379 Tau-F missense
variants and only 2 strict ClinVar binary benchmark variants, so its AUROC/AUPRC
should be described as an underpowered sanity check rather than a performance
claim.
