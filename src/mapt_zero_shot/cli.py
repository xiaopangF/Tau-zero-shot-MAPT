"""Command-line interface for the MAPT zero-shot atlas workflow."""

from __future__ import annotations

import argparse
from pathlib import Path

from .annotations import annotate_variant_row
from .clinvar import load_mapt_clinvar_with_qc
from .ensemble import ensemble_score_rows
from .evaluate import make_binary_examples, metrics_rows
from .external_scores import load_alphamissense
from .figures import create_basic_figures
from .io import merge_fieldnames, read_tsv, write_tsv
from .prioritize import prioritize_rows
from .score_esm import masked_marginal_scores
from .variants import generate_missense_variants


def cmd_enumerate(args: argparse.Namespace) -> None:
    rows = [annotate_variant_row(variant.as_row()) for variant in generate_missense_variants()]
    fieldnames = [
        "variant_id",
        "protein_change",
        "position",
        "wt_aa",
        "mut_aa",
        "tau_region",
        "tau_motif",
        "near_ptm_site_3aa",
        "nearest_ptm_distance",
        "near_known_pathogenic_hotspot_3aa",
        "nearest_hotspot_distance",
        "charge_change",
        "special_residue_change",
    ]
    write_tsv(args.out, rows, fieldnames)
    print(f"Wrote {len(rows)} variants to {args.out}")


def cmd_import_clinvar(args: argparse.Namespace) -> None:
    variants = read_tsv(args.variants)
    result = load_mapt_clinvar_with_qc(args.clinvar)
    clinvar_by_id = result.accepted_by_variant
    rows = []
    for row in variants:
        merged = dict(row)
        merged.update(clinvar_by_id.get(row["variant_id"], {}))
        rows.append(merged)
    fieldnames = merge_fieldnames(rows)
    write_tsv(args.out, rows, fieldnames)
    if args.rejected_out:
        rejected_fieldnames = merge_fieldnames(result.rejected_rows)
        write_tsv(args.rejected_out, result.rejected_rows, rejected_fieldnames)
        print(f"Wrote {len(result.rejected_rows)} rejected ClinVar rows to {args.rejected_out}")
    print(f"Wrote {len(rows)} annotated variants to {args.out}")


def cmd_import_alphamissense(args: argparse.Namespace) -> None:
    variants = read_tsv(args.variants)
    alpha_by_id = load_alphamissense(
        args.alphamissense,
        transcript_id=args.transcript_id,
        uniprot_id=args.uniprot_id,
        position_offset=args.position_offset,
    )
    rows = []
    for row in variants:
        merged = dict(row)
        merged.update(alpha_by_id.get(row["variant_id"], {}))
        rows.append(merged)
    fieldnames = merge_fieldnames(rows)
    write_tsv(args.out, rows, fieldnames)
    print(f"Wrote {len(rows)} variants with AlphaMissense fields to {args.out}")


def cmd_score_esm(args: argparse.Namespace) -> None:
    rows = masked_marginal_scores(model_name=args.model, device=args.device)
    fieldnames = [
        "variant_id",
        "protein_change",
        "position",
        "wt_aa",
        "mut_aa",
        "model",
        "esm_wt_logp",
        "esm_mut_logp",
        "esm_llr",
        "pathogenic_score",
    ]
    write_tsv(args.out, rows, fieldnames)
    print(f"Wrote {len(rows)} ESM scores to {args.out}")


def cmd_ensemble(args: argparse.Namespace) -> None:
    tables = [read_tsv(path) for path in args.scores]
    rows = ensemble_score_rows(tables, args.score_column)
    mean_column = f"{args.score_column}_mean"
    std_column = f"{args.score_column}_std"
    fieldnames = [
        "variant_id",
        "protein_change",
        "position",
        "wt_aa",
        "mut_aa",
        "ensemble_score_column",
        "n_models",
        "models",
        mean_column,
        std_column,
    ]
    write_tsv(args.out, rows, fieldnames)
    print(f"Wrote {len(rows)} ensemble scores to {args.out}")


def cmd_evaluate(args: argparse.Namespace) -> None:
    scores = read_tsv(args.scores)
    labels = read_tsv(args.labels)
    examples = make_binary_examples(scores, labels, args.score_column)
    rows = metrics_rows(examples)
    write_tsv(args.out, rows, ["metric", "value"])
    print(f"Wrote metrics for {len(examples)} benchmark variants to {args.out}")


def cmd_prioritize(args: argparse.Namespace) -> None:
    scores = read_tsv(args.scores)
    annotations = read_tsv(args.annotations)
    labels = set(args.label)
    rows = prioritize_rows(scores, annotations, args.score_column, labels, args.include_unlabeled)
    if args.limit is not None:
        rows = rows[: args.limit]
    fieldnames = merge_fieldnames(rows)
    write_tsv(args.out, rows, fieldnames)
    print(f"Wrote {len(rows)} prioritized variants to {args.out}")


def cmd_figures(args: argparse.Namespace) -> None:
    scores = read_tsv(args.scores)
    labels = read_tsv(args.labels)
    paths = create_basic_figures(scores, labels, args.score_column, args.outdir)
    for path in paths:
        print(f"Wrote {path}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    enumerate_parser = sub.add_parser("enumerate", help="Generate all MAPT missense variants.")
    enumerate_parser.add_argument("--out", required=True, type=Path)
    enumerate_parser.set_defaults(func=cmd_enumerate)

    clinvar_parser = sub.add_parser("import-clinvar", help="Merge ClinVar labels into variants.")
    clinvar_parser.add_argument("--variants", required=True)
    clinvar_parser.add_argument("--clinvar", required=True)
    clinvar_parser.add_argument("--out", required=True, type=Path)
    clinvar_parser.add_argument("--rejected-out", default=None, type=Path)
    clinvar_parser.set_defaults(func=cmd_import_clinvar)

    alpha_parser = sub.add_parser("import-alphamissense", help="Merge AlphaMissense scores.")
    alpha_parser.add_argument("--variants", required=True)
    alpha_parser.add_argument("--alphamissense", required=True)
    alpha_parser.add_argument("--out", required=True, type=Path)
    alpha_parser.add_argument("--transcript-id", default=None)
    alpha_parser.add_argument("--uniprot-id", default=None)
    alpha_parser.add_argument("--position-offset", type=int, default=0)
    alpha_parser.set_defaults(func=cmd_import_alphamissense)

    score_parser = sub.add_parser("score-esm", help="Score all variants using fair-esm.")
    score_parser.add_argument("--model", default="esm1v_t33_650M_UR90S_1")
    score_parser.add_argument("--device", default=None)
    score_parser.add_argument("--out", required=True, type=Path)
    score_parser.set_defaults(func=cmd_score_esm)

    ensemble_parser = sub.add_parser("ensemble", help="Average multiple score TSV files.")
    ensemble_parser.add_argument("--scores", nargs="+", required=True)
    ensemble_parser.add_argument("--score-column", default="pathogenic_score")
    ensemble_parser.add_argument("--out", required=True, type=Path)
    ensemble_parser.set_defaults(func=cmd_ensemble)

    eval_parser = sub.add_parser("evaluate", help="Evaluate scores against ClinVar P/LP vs B/LB.")
    eval_parser.add_argument("--scores", required=True)
    eval_parser.add_argument("--labels", required=True)
    eval_parser.add_argument("--score-column", default="pathogenic_score")
    eval_parser.add_argument("--out", required=True, type=Path)
    eval_parser.set_defaults(func=cmd_evaluate)

    priority_parser = sub.add_parser("prioritize", help="Rank VUS or unlabeled variants.")
    priority_parser.add_argument("--scores", required=True)
    priority_parser.add_argument("--annotations", required=True)
    priority_parser.add_argument("--score-column", default="pathogenic_score")
    priority_parser.add_argument("--label", action="append", default=["VUS"])
    priority_parser.add_argument("--include-unlabeled", action="store_true")
    priority_parser.add_argument("--limit", type=int, default=None)
    priority_parser.add_argument("--out", required=True, type=Path)
    priority_parser.set_defaults(func=cmd_prioritize)

    fig_parser = sub.add_parser("figures", help="Create basic publication figures.")
    fig_parser.add_argument("--scores", required=True)
    fig_parser.add_argument("--labels", required=True)
    fig_parser.add_argument("--score-column", default="pathogenic_score")
    fig_parser.add_argument("--outdir", required=True, type=Path)
    fig_parser.set_defaults(func=cmd_figures)

    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()

