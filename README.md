# MAPT/Tau Zero-Shot Variant Atlas

This project builds a reproducible, label-free prioritization atlas for MAPT/Tau
missense variants. The v1 scope is a computational bioinformatics paper:
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

## Manuscript claim

Recommended wording:

> We generated a zero-shot MAPT/Tau missense variant atlas using pretrained
> protein sequence models and evaluated whether label-free scores prioritize
> known pathogenic MAPT variants and clinically uncertain missense variants.


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

Use `pathogenic_score` for all main manuscript ranking analyses.

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

Use the guarded script for the manuscript-scale ESM-1v ensemble:

```powershell
.\scripts\run_esm1v_ensemble.ps1
```

The script refuses to run ESM-1v on CPU by default. If you intentionally want a
long CPU run, pass `-AllowCpu`.
