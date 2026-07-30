"""Command-line interface for the MAPT zero-shot atlas workflow."""

from __future__ import annotations

import argparse
from math import ceil
from pathlib import Path

from .annotations import annotate_variant_row
from .baselines import heuristic_score_rows
from .calibration import (
    create_positive_control_plot,
    positive_control_rows,
    positive_control_summary_rows,
)
from .compare import compare_model_rows
from .concordance import (
    concordance_rows,
    concordance_summary_rows,
    top_concordance_rows,
)
from .clinvar import load_mapt_clinvar_with_qc
from .ensemble import ensemble_score_rows
from .evaluate import make_binary_examples, metrics_rows
from .external_scores import load_alphamissense_with_qc
from .figures import create_basic_figures
from .fusion import composite_summary_rows, composite_zero_shot_rows
from .io import merge_fieldnames, read_tsv, write_tsv
from .knowledge_a import score_scheme_a_rows, scheme_a_gold_rows, scheme_a_summary_rows
from .physchem import (
    feature_rows,
    grid_search_weights,
    ranked_physchem_rows,
    ranked_zero_shot_physchem_rows,
    summary_rows as physchem_summary_rows,
    zero_shot_physchem_rows,
    zero_shot_summary_rows,
)
from .manuscript_assets import (
    create_domain_summary_plot,
    create_model_comparison_plot,
    create_vus_lollipop_plot,
    create_score_scatter_plot,
    create_workflow_schematic,
    write_results_summary,
    write_top_vus_markdown,
)
from .prioritize import prioritize_rows
from .score_esm import masked_marginal_scores
from .submission_qc import submission_qc_rows, write_submission_qc_markdown
from .summaries import alphamissense_summary_rows, clinvar_summary_rows, domain_summary_rows
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
    result = load_alphamissense_with_qc(
        args.alphamissense,
        transcript_id=args.transcript_id,
        uniprot_id=args.uniprot_id,
        position_offset=args.position_offset,
    )
    alpha_by_id = result.accepted_by_variant
    rows = []
    for row in variants:
        merged = dict(row)
        merged.update(alpha_by_id.get(row["variant_id"], {}))
        rows.append(merged)
    fieldnames = merge_fieldnames(rows)
    write_tsv(args.out, rows, fieldnames)
    if args.rejected_out:
        rejected_fieldnames = merge_fieldnames(result.rejected_rows)
        write_tsv(args.rejected_out, result.rejected_rows, rejected_fieldnames)
        print(
            f"Wrote {len(result.rejected_rows)} rejected AlphaMissense rows to "
            f"{args.rejected_out}"
        )
    print(f"Wrote {len(rows)} variants with AlphaMissense fields to {args.out}")


def cmd_score_physchem_zero_shot(args: argparse.Namespace) -> None:
    variants = read_tsv(args.variants)
    rows = feature_rows(variants)
    scored = zero_shot_physchem_rows(rows)
    ranked = ranked_zero_shot_physchem_rows(scored, controls=tuple(args.control))
    fieldnames = [
        "variant_id",
        "protein_change",
        "position",
        "wt_aa",
        "mut_aa",
        "tau_region",
        "tau_motif",
        "near_ptm_site_3aa",
        "motif_prior",
        "ptm_prior",
        "mechanistic_prior_score",
        "hydrophobicity_delta",
        "beta_sheet_delta",
        "net_charge_delta",
        "standardized_hydrophobicity_delta",
        "standardized_beta_sheet_delta",
        "standardized_net_charge_delta",
        "region_multiplier",
        "physchem_perturbation_score",
        "physchem_zero_shot_score",
        "physchem_rank",
        "physchem_top_fraction",
        "physchem_top_percent",
        "gold_standard_control",
    ]
    write_tsv(args.out, ranked, fieldnames)
    write_tsv(
        args.summary_out,
        zero_shot_summary_rows(ranked, controls=tuple(args.control), top_fraction=args.top_fraction),
        ["metric", "value"],
    )
    control_rows = [row for row in ranked if row["gold_standard_control"] == "True"]
    mean_rank = sum(int(row["physchem_rank"]) for row in control_rows) / len(control_rows)
    top_cutoff = ceil(len(ranked) * args.top_fraction)
    controls_in_top = sum(int(row["physchem_rank"]) <= top_cutoff for row in control_rows)
    print("Model: zero-shot physicochemical perturbation (no label-informed tuning)")
    print(
        f"Gold-standard controls in top {args.top_fraction:.2%}: "
        f"{controls_in_top}/{len(control_rows)}; mean rank={mean_rank:.2f}"
    )
    print(f"Wrote ranked variants to {args.out}")

def cmd_score_physchem_grid(args: argparse.Namespace) -> None:
    variants = read_tsv(args.variants)
    rows = feature_rows(variants)
    result = grid_search_weights(
        rows,
        controls=tuple(args.control),
        grid_min=args.grid_min,
        grid_max=args.grid_max,
        grid_step=args.grid_step,
        top_fraction=args.top_fraction,
    )
    ranked = ranked_physchem_rows(rows, result.weights, controls=tuple(args.control))
    fieldnames = [
        "variant_id",
        "protein_change",
        "position",
        "wt_aa",
        "mut_aa",
        "hydrophobicity_delta",
        "beta_sheet_delta",
        "net_charge_delta",
        "physchem_score",
        "physchem_rank",
        "physchem_top_fraction",
        "physchem_top_percent",
        "gold_standard_control",
    ]
    write_tsv(args.out, ranked, fieldnames)
    write_tsv(
        args.summary_out,
        physchem_summary_rows(
            result, args.grid_min, args.grid_max, args.grid_step, args.top_fraction
        ),
        ["metric", "value"],
    )
    print(
        "Best weights: "
        f"hydrophobicity={result.weights.hydrophobicity_delta}, "
        f"beta_sheet={result.weights.beta_sheet_delta}, "
        f"net_charge={result.weights.net_charge_delta}"
    )
    print(
        f"Gold-standard controls in top {args.top_fraction:.2%}: "
        f"{result.controls_in_top_fraction}/{len(result.control_ranks)}; "
        f"mean rank={result.mean_control_rank:.2f}"
    )
    print(f"Wrote ranked variants to {args.out}")

def cmd_fuse_zero_shot(args: argparse.Namespace) -> None:
    esm_rows = read_tsv(args.esm_scores)
    physchem_rows = read_tsv(args.physchem_scores)
    ranked = composite_zero_shot_rows(
        esm_rows,
        physchem_rows,
        esm_column=args.esm_column,
        physchem_column=args.physchem_column,
        esm_weight=args.esm_weight,
        physchem_weight=args.physchem_weight,
    )
    control_set = set(args.control)
    for row in ranked:
        row["gold_standard_control"] = str(row["variant_id"] in control_set)
    fieldnames = [
        "variant_id",
        "protein_change",
        "position",
        "wt_aa",
        "mut_aa",
        "tau_region",
        "tau_motif",
        "near_ptm_site_3aa",
        "esm_raw_score",
        "physchem_raw_score",
        "esm_standardized_score",
        "physchem_standardized_score",
        "esm_weight",
        "physchem_weight",
        "zero_shot_composite_score",
        "zero_shot_composite_rank",
        "zero_shot_composite_top_fraction",
        "zero_shot_composite_top_percent",
        "gold_standard_control",
    ]
    write_tsv(args.out, ranked, fieldnames)
    write_tsv(
        args.summary_out,
        composite_summary_rows(ranked, controls=tuple(args.control), top_fraction=args.top_fraction),
        ["metric", "value"],
    )
    control_rows = [row for row in ranked if row["gold_standard_control"] == "True"]
    mean_rank = sum(int(row["zero_shot_composite_rank"]) for row in control_rows) / len(control_rows)
    top_cutoff = ceil(len(ranked) * args.top_fraction)
    controls_in_top = sum(
        int(row["zero_shot_composite_rank"]) <= top_cutoff for row in control_rows
    )
    print("Model: zero-shot ESM + physicochemical fusion (no label-informed tuning)")
    print(
        f"Gold-standard controls in top {args.top_fraction:.2%}: "
        f"{controls_in_top}/{len(control_rows)}; mean rank={mean_rank:.2f}"
    )
    print(f"Wrote ranked variants to {args.out}")

def cmd_score_knowledge_a(args: argparse.Namespace) -> None:
    feature_table = read_tsv(args.physchem_features)
    ranked = score_scheme_a_rows(feature_table, controls=tuple(args.control))
    fieldnames = [
        "variant_id",
        "protein_change",
        "position",
        "wt_aa",
        "mut_aa",
        "tau_region",
        "hydrophobicity_delta",
        "beta_sheet_delta",
        "net_charge_delta",
        "baseline_physchem_z",
        "directed_beta_score",
        "ptm_microenvironment_score",
        "splicing_proximity_score",
        "microtubule_interface_disturbance_strength",
        "microtubule_interface_score",
        "scheme_a_final_score",
        "scheme_a_rank",
        "scheme_a_top_fraction",
        "scheme_a_top_percent",
        "gold_standard_control",
    ]
    gold_fieldnames = [
        "variant_id",
        "baseline_physchem_z",
        "directed_beta_score",
        "ptm_microenvironment_score",
        "splicing_proximity_score",
        "microtubule_interface_disturbance_strength",
        "microtubule_interface_score",
        "scheme_a_final_score",
        "scheme_a_rank",
    ]
    write_tsv(args.out, ranked, fieldnames)
    gold_rows = scheme_a_gold_rows(ranked, controls=tuple(args.control))
    write_tsv(args.gold_out, gold_rows, gold_fieldnames)
    write_tsv(
        args.summary_out,
        scheme_a_summary_rows(ranked, controls=tuple(args.control), top_fraction=args.top_fraction),
        ["metric", "value"],
    )
    top_cutoff = ceil(len(ranked) * args.top_fraction)
    controls_in_top = sum(int(row["scheme_a_rank"]) <= top_cutoff for row in gold_rows)
    mean_rank = sum(int(row["scheme_a_rank"]) for row in gold_rows) / len(gold_rows)
    print("Model: Scheme A knowledge-driven score (fixed parameters, no label tuning)")
    print(
        f"Gold-standard controls in top {args.top_fraction:.2%}: "
        f"{controls_in_top}/{len(gold_rows)}; mean rank={mean_rank:.2f}"
    )
    print(f"Wrote ranked variants to {args.out}")
    print(f"Wrote gold-standard detail table to {args.gold_out}")

def cmd_score_heuristic(args: argparse.Namespace) -> None:
    annotations = read_tsv(args.annotations)
    rows = heuristic_score_rows(annotations)
    fieldnames = merge_fieldnames(rows)
    write_tsv(args.out, rows, fieldnames)
    print(f"Wrote {len(rows)} heuristic baseline scores to {args.out}")

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


def _parse_model_spec(spec: str) -> tuple[str, str, str]:
    parts = spec.split(":")
    if len(parts) != 3:
        raise ValueError("Model specs must use name:path:score_column format")
    return parts[0], parts[1], parts[2]


def cmd_compare_models(args: argparse.Namespace) -> None:
    labels = read_tsv(args.labels)
    specs = []
    for spec in args.model:
        model_name, path, score_column = _parse_model_spec(spec)
        specs.append((model_name, read_tsv(path), score_column))
    rows = compare_model_rows(specs, labels)
    write_tsv(args.out, rows, ["model", "score_column", "metric", "value"])
    print(f"Wrote {len(rows)} model comparison rows to {args.out}")

def cmd_summarize_clinvar(args: argparse.Namespace) -> None:
    benchmark = read_tsv(args.benchmark)
    rejected = read_tsv(args.rejected) if args.rejected else None
    rows = clinvar_summary_rows(benchmark, rejected)
    write_tsv(args.out, rows, ["section", "category", "count"])
    print(f"Wrote {len(rows)} ClinVar summary rows to {args.out}")


def cmd_summarize_alphamissense(args: argparse.Namespace) -> None:
    alpha_rows = read_tsv(args.alpha)
    rejected_rows = read_tsv(args.rejected) if args.rejected else None
    rows = alphamissense_summary_rows(alpha_rows, rejected_rows)
    write_tsv(args.out, rows, ["section", "category", "count"])
    print(f"Wrote {len(rows)} AlphaMissense summary rows to {args.out}")

def cmd_calibrate_positive_controls(args: argparse.Namespace) -> None:
    scores = read_tsv(args.scores)
    annotations = read_tsv(args.annotations) if args.annotations else []
    rows = positive_control_rows(scores, annotations, score_column=args.score_column)
    write_tsv(
        args.out,
        rows,
        [
            "variant_id",
            "protein_change",
            "position",
            "tau_region",
            "evidence_class",
            "evidence_note",
            "found_in_atlas",
            "score",
            "score_std",
            "atlas_size",
            "rank",
            "top_fraction",
            "top_percent",
            "top_1pct",
            "top_5pct",
            "top_10pct",
            "calibration_interpretation",
        ],
    )
    write_tsv(args.summary_out, positive_control_summary_rows(rows), ["metric", "value"])
    if args.figure:
        create_positive_control_plot(rows, args.figure)
    print(f"Wrote {len(rows)} positive-control rows to {args.out}")


def cmd_summarize_domains(args: argparse.Namespace) -> None:
    scores = read_tsv(args.scores)
    annotations = read_tsv(args.annotations)
    rows = domain_summary_rows(scores, annotations, args.score_column)
    fieldnames = [
        "tau_region",
        "n_variants",
        "mean_score",
        "median_score",
        "q25_score",
        "q75_score",
        "max_score",
        "top_5pct_global_count",
        "top_5pct_global_fraction",
        "top_10pct_global_count",
        "top_10pct_global_fraction",
    ]
    write_tsv(args.out, rows, fieldnames)
    print(f"Wrote {len(rows)} domain summary rows to {args.out}")


def cmd_manuscript_assets(args: argparse.Namespace) -> None:
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    domain_rows = read_tsv(args.domain_summary)
    comparison_rows = read_tsv(args.model_comparison)
    vus_rows = read_tsv(args.vus_priority)
    esm_rows = read_tsv(args.esm_scores)
    heuristic_rows = read_tsv(args.heuristic_scores)
    clinvar_summary_rows_input = read_tsv(args.clinvar_summary) if args.clinvar_summary else None
    alphamissense_summary_rows_input = (
        read_tsv(args.alphamissense_summary) if args.alphamissense_summary else None
    )

    create_workflow_schematic(
        outdir / "figure_workflow_schematic.png",
        atlas_variant_count=len(esm_rows),
        clinvar_summary_rows=clinvar_summary_rows_input,
        alphamissense_summary_rows=alphamissense_summary_rows_input,
    )
    create_domain_summary_plot(domain_rows, outdir / "figure_domain_summary.png")
    create_model_comparison_plot(comparison_rows, outdir / "figure_model_comparison.png")
    _, scatter_r, scatter_n = create_score_scatter_plot(
        esm_rows, heuristic_rows, outdir / "figure_esm1v_vs_heuristic.png"
    )
    create_vus_lollipop_plot(vus_rows, outdir / "figure_vus_lollipop.png", args.top_n)
    write_top_vus_markdown(vus_rows, outdir / "top_vus_candidates.md", args.top_n)
    write_results_summary(
        domain_rows,
        comparison_rows,
        vus_rows,
        scatter_r,
        scatter_n,
        outdir / "results_summary.md",
    )
    print(f"Wrote manuscript assets to {outdir}")


def cmd_concordance(args: argparse.Namespace) -> None:
    annotations = read_tsv(args.annotations)
    esm_rows = read_tsv(args.esm_scores)
    heuristic_rows = read_tsv(args.heuristic_scores)
    alphamissense_rows = read_tsv(args.alphamissense_scores)
    rows = concordance_rows(
        annotations,
        esm_rows,
        heuristic_rows,
        alphamissense_rows,
        esm_column=args.esm_column,
        heuristic_column=args.heuristic_column,
        alphamissense_column=args.alphamissense_column,
        top_fraction=args.top_fraction,
    )
    fieldnames = [
        "variant_id",
        "protein_change",
        "position",
        "wt_aa",
        "mut_aa",
        "tau_region",
        "tau_motif",
        "charge_change",
        "special_residue_change",
        "esm_score",
        "heuristic_score",
        "alphamissense_score",
        "esm_top_decile",
        "heuristic_top_decile",
        "alphamissense_top_decile",
        "alphamissense_has_score",
        "n_models_high",
        "concordance_category",
    ]
    write_tsv(args.out, rows, fieldnames)
    write_tsv(
        args.summary_out,
        concordance_summary_rows(rows),
        ["section", "category", "count"],
    )
    top_rows = top_concordance_rows(rows, category=args.category, limit=args.limit)
    write_tsv(args.top_out, top_rows, fieldnames)
    print(f"Wrote {len(rows)} concordance rows to {args.out}")
    print(f"Wrote {len(top_rows)} top concordance rows to {args.top_out}")

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
    alpha_parser.add_argument("--rejected-out", default=None, type=Path)
    alpha_parser.add_argument("--transcript-id", default=None)
    alpha_parser.add_argument("--uniprot-id", default=None)
    alpha_parser.add_argument("--position-offset", type=int, default=0)
    alpha_parser.set_defaults(func=cmd_import_alphamissense)

    physchem_zero_shot_parser = sub.add_parser(
        "score-physchem-zero-shot",
        help="Score physicochemical disruption with a fixed label-free model.",
    )
    physchem_zero_shot_parser.add_argument("--variants", required=True)
    physchem_zero_shot_parser.add_argument("--out", required=True, type=Path)
    physchem_zero_shot_parser.add_argument("--summary-out", required=True, type=Path)
    physchem_zero_shot_parser.add_argument("--top-fraction", type=float, default=0.01)
    physchem_zero_shot_parser.add_argument(
        "--control",
        action="append",
        default=["G272V", "P301L", "V337M", "R406W", "N279K"],
        help="Gold-standard control used for validation only. May be repeated.",
    )
    physchem_zero_shot_parser.set_defaults(func=cmd_score_physchem_zero_shot)
    physchem_parser = sub.add_parser(
        "score-physchem-grid",
        help="Grid-search hydrophobicity, beta-sheet, and charge weights.",
    )
    physchem_parser.add_argument("--variants", required=True)
    physchem_parser.add_argument("--out", required=True, type=Path)
    physchem_parser.add_argument("--summary-out", required=True, type=Path)
    physchem_parser.add_argument("--grid-min", type=float, default=-5.0)
    physchem_parser.add_argument("--grid-max", type=float, default=5.0)
    physchem_parser.add_argument("--grid-step", type=float, default=0.5)
    physchem_parser.add_argument("--top-fraction", type=float, default=0.01)
    physchem_parser.add_argument(
        "--control",
        action="append",
        default=["G272V", "P301L", "V337M", "R406W", "N279K"],
        help="Gold-standard control variant ID. May be repeated.",
    )
    physchem_parser.set_defaults(func=cmd_score_physchem_grid)
    fusion_parser = sub.add_parser(
        "fuse-zero-shot",
        help="Fuse ESM and physicochemical zero-shot scores after unlabeled normalization.",
    )
    fusion_parser.add_argument("--esm-scores", required=True)
    fusion_parser.add_argument("--physchem-scores", required=True)
    fusion_parser.add_argument("--esm-column", default="pathogenic_score_mean")
    fusion_parser.add_argument("--physchem-column", default="physchem_zero_shot_score")
    fusion_parser.add_argument("--esm-weight", type=float, default=0.75)
    fusion_parser.add_argument("--physchem-weight", type=float, default=0.25)
    fusion_parser.add_argument("--top-fraction", type=float, default=0.01)
    fusion_parser.add_argument("--out", required=True, type=Path)
    fusion_parser.add_argument("--summary-out", required=True, type=Path)
    fusion_parser.add_argument(
        "--control",
        action="append",
        default=["G272V", "P301L", "V337M", "R406W", "N279K"],
        help="Gold-standard control used for validation only. May be repeated.",
    )
    fusion_parser.set_defaults(func=cmd_fuse_zero_shot)
    knowledge_a_parser = sub.add_parser(
        "score-knowledge-a",
        help="Score variants with fixed Scheme A knowledge-driven rules.",
    )
    knowledge_a_parser.add_argument("--physchem-features", required=True)
    knowledge_a_parser.add_argument("--out", required=True, type=Path)
    knowledge_a_parser.add_argument("--gold-out", required=True, type=Path)
    knowledge_a_parser.add_argument("--summary-out", required=True, type=Path)
    knowledge_a_parser.add_argument("--top-fraction", type=float, default=0.01)
    knowledge_a_parser.add_argument(
        "--control",
        action="append",
        default=["G272V", "P301L", "V337M", "R406W", "N279K"],
        help="Gold-standard control used for validation only. May be repeated.",
    )
    knowledge_a_parser.set_defaults(func=cmd_score_knowledge_a)
    heuristic_parser = sub.add_parser(
        "score-heuristic", help="Score variants with a transparent Tau heuristic baseline."
    )
    heuristic_parser.add_argument("--annotations", required=True)
    heuristic_parser.add_argument("--out", required=True, type=Path)
    heuristic_parser.set_defaults(func=cmd_score_heuristic)
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

    compare_parser = sub.add_parser(
        "compare-models", help="Compare multiple score files against ClinVar labels."
    )
    compare_parser.add_argument("--labels", required=True)
    compare_parser.add_argument(
        "--model", action="append", required=True, help="Use name:path:score_column format."
    )
    compare_parser.add_argument("--out", required=True, type=Path)
    compare_parser.set_defaults(func=cmd_compare_models)
    clinvar_summary_parser = sub.add_parser(
        "summarize-clinvar", help="Summarize accepted and rejected ClinVar rows."
    )
    clinvar_summary_parser.add_argument("--benchmark", required=True)
    clinvar_summary_parser.add_argument("--rejected", default=None)
    clinvar_summary_parser.add_argument("--out", required=True, type=Path)
    clinvar_summary_parser.set_defaults(func=cmd_summarize_clinvar)
    alpha_summary_parser = sub.add_parser(
        "summarize-alphamissense", help="Summarize AlphaMissense coverage and QC rejection reasons."
    )
    alpha_summary_parser.add_argument("--alpha", required=True)
    alpha_summary_parser.add_argument("--rejected", default=None)
    alpha_summary_parser.add_argument("--out", required=True, type=Path)
    alpha_summary_parser.set_defaults(func=cmd_summarize_alphamissense)

    calibration_parser = sub.add_parser(
        "calibrate-positive-controls",
        help="Rank established MAPT controls within the complete ESM atlas.",
    )
    calibration_parser.add_argument("--scores", required=True)
    calibration_parser.add_argument("--annotations", default=None)
    calibration_parser.add_argument("--score-column", default="pathogenic_score_mean")
    calibration_parser.add_argument("--out", required=True, type=Path)
    calibration_parser.add_argument("--summary-out", required=True, type=Path)
    calibration_parser.add_argument("--figure", default=None, type=Path)
    calibration_parser.set_defaults(func=cmd_calibrate_positive_controls)

    domain_summary_parser = sub.add_parser(
        "summarize-domains", help="Summarize score distributions by Tau region."
    )
    domain_summary_parser.add_argument("--scores", required=True)
    domain_summary_parser.add_argument("--annotations", required=True)
    domain_summary_parser.add_argument("--score-column", default="pathogenic_score")
    domain_summary_parser.add_argument("--out", required=True, type=Path)
    domain_summary_parser.set_defaults(func=cmd_summarize_domains)

    assets_parser = sub.add_parser(
        "manuscript-assets", help="Create manuscript summary figures and markdown."
    )
    assets_parser.add_argument("--domain-summary", required=True)
    assets_parser.add_argument("--model-comparison", required=True)
    assets_parser.add_argument("--vus-priority", required=True)
    assets_parser.add_argument("--esm-scores", required=True)
    assets_parser.add_argument("--heuristic-scores", required=True)
    assets_parser.add_argument("--clinvar-summary", default=None)
    assets_parser.add_argument("--alphamissense-summary", default=None)
    assets_parser.add_argument("--top-n", type=int, default=10)
    assets_parser.add_argument("--outdir", required=True, type=Path)
    assets_parser.set_defaults(func=cmd_manuscript_assets)
    concordance_parser = sub.add_parser(
        "concordance", help="Analyze agreement and disagreement across scoring methods."
    )
    concordance_parser.add_argument("--annotations", required=True)
    concordance_parser.add_argument("--esm-scores", required=True)
    concordance_parser.add_argument("--heuristic-scores", required=True)
    concordance_parser.add_argument("--alphamissense-scores", required=True)
    concordance_parser.add_argument("--esm-column", default="pathogenic_score_mean")
    concordance_parser.add_argument("--heuristic-column", default="heuristic_score")
    concordance_parser.add_argument("--alphamissense-column", default="alphamissense_score")
    concordance_parser.add_argument("--top-fraction", type=float, default=0.10)
    concordance_parser.add_argument("--category", default=None)
    concordance_parser.add_argument("--limit", type=int, default=50)
    concordance_parser.add_argument("--out", required=True, type=Path)
    concordance_parser.add_argument("--summary-out", required=True, type=Path)
    concordance_parser.add_argument("--top-out", required=True, type=Path)
    concordance_parser.set_defaults(func=cmd_concordance)
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
