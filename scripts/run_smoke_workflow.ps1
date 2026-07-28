$ErrorActionPreference = "Stop"
$env:PYTHONPATH = "src"

python -m mapt_zero_shot.cli enumerate --out data/processed/mapt_all_missense_variants.tsv
python -m mapt_zero_shot.cli score-esm --model esm2_t6_8M_UR50D --out results/scores/mapt_esm2_t6_8M_scores.tsv

if (Test-Path data/raw/variant_summary.txt.gz) {
    python -m mapt_zero_shot.cli import-clinvar --variants data/processed/mapt_all_missense_variants.tsv --clinvar data/raw/variant_summary.txt.gz --out data/processed/mapt_clinvar_benchmark.tsv --rejected-out data/processed/mapt_clinvar_rejected.tsv
    python -m mapt_zero_shot.cli evaluate --scores results/scores/mapt_esm2_t6_8M_scores.tsv --labels data/processed/mapt_clinvar_benchmark.tsv --score-column pathogenic_score --out results/mapt_esm2_t6_8M_metrics.tsv
    python -m mapt_zero_shot.cli summarize-clinvar --benchmark data/processed/mapt_clinvar_benchmark.tsv --rejected data/processed/mapt_clinvar_rejected.tsv --out results/mapt_clinvar_qc_summary.tsv
    python -m mapt_zero_shot.cli summarize-domains --scores results/scores/mapt_esm2_t6_8M_scores.tsv --annotations data/processed/mapt_all_missense_variants.tsv --score-column pathogenic_score --out results/mapt_esm2_t6_8M_domain_summary.tsv
    python -m mapt_zero_shot.cli prioritize --scores results/scores/mapt_esm2_t6_8M_scores.tsv --annotations data/processed/mapt_clinvar_benchmark.tsv --score-column pathogenic_score --label VUS --limit 50 --out results/mapt_top50_vus_priority.tsv
    python -m mapt_zero_shot.cli figures --scores results/scores/mapt_esm2_t6_8M_scores.tsv --labels data/processed/mapt_clinvar_benchmark.tsv --score-column pathogenic_score --outdir results/figures
} else {
    python -m mapt_zero_shot.cli summarize-domains --scores results/scores/mapt_esm2_t6_8M_scores.tsv --annotations data/processed/mapt_all_missense_variants.tsv --score-column pathogenic_score --out results/mapt_esm2_t6_8M_domain_summary.tsv
    python -m mapt_zero_shot.cli prioritize --scores results/scores/mapt_esm2_t6_8M_scores.tsv --annotations data/processed/mapt_all_missense_variants.tsv --score-column pathogenic_score --include-unlabeled --limit 50 --out results/mapt_top50_unlabeled_priority.tsv
    python -m mapt_zero_shot.cli figures --scores results/scores/mapt_esm2_t6_8M_scores.tsv --labels data/processed/mapt_all_missense_variants.tsv --score-column pathogenic_score --outdir results/figures
}
