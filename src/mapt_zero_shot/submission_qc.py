"""Submission-readiness checks for manuscript drafts."""

from __future__ import annotations

from pathlib import Path


NEGATION_MARKERS = ("not", "avoid", "should not", "rather than", "without")


RISKY_PHRASES = (
    "accurately predicts pathogenicity",
    "clinically validates",
    "clinical diagnostic model",
    "diagnostic test",
    "reclassified pathogenic",
)


def _summary_lookup(rows: list[dict[str, str]], section: str, category: str) -> str:
    for row in rows:
        if row.get("section") == section and row.get("category") == category:
            return row.get("count", "")
    return ""


def _check_number(text: str, label: str, value: object) -> dict[str, object]:
    value_text = str(value)
    return {
        "check": label,
        "status": "pass" if value_text in text else "fail",
        "detail": f"expected manuscript to contain {value_text}",
    }


def _check_file(path: Path) -> dict[str, object]:
    return {
        "check": f"file_exists:{path.as_posix()}",
        "status": "pass" if path.exists() else "fail",
        "detail": "present" if path.exists() else "missing",
    }


def submission_qc_rows(
    manuscript_text: str,
    variants: list[dict[str, str]],
    clinvar_summary: list[dict[str, str]],
    alphamissense_summary: list[dict[str, str]],
    concordance_summary: list[dict[str, str]],
    required_files: list[Path] | None = None,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    required_files = required_files or []

    rows.append(_check_number(manuscript_text, "variant_count", len(variants)))
    rows.append(
        _check_number(
            manuscript_text,
            "clinvar_accepted",
            _summary_lookup(clinvar_summary, "accepted", "total_annotated_variants"),
        )
    )
    rows.append(
        _check_number(
            manuscript_text,
            "clinvar_rejected",
            _summary_lookup(clinvar_summary, "rejected", "total_rejected_rows"),
        )
    )
    rows.append(
        _check_number(
            manuscript_text,
            "alphamissense_scored",
            _summary_lookup(alphamissense_summary, "coverage", "scored_variants"),
        )
    )
    rows.append(
        _check_number(
            manuscript_text,
            "concordance_esm_heuristic_high",
            _summary_lookup(concordance_summary, "concordance_category", "esm_heuristic_high"),
        )
    )
    rows.append(
        _check_number(
            manuscript_text,
            "concordance_esm_only_high",
            _summary_lookup(concordance_summary, "concordance_category", "esm_only_high"),
        )
    )
    rows.append(
        _check_number(
            manuscript_text,
            "concordance_heuristic_only_high",
            _summary_lookup(concordance_summary, "concordance_category", "heuristic_only_high"),
        )
    )

    lowered = manuscript_text.lower()
    for phrase in RISKY_PHRASES:
        phrase_lower = phrase.lower()
        count = lowered.count(phrase_lower)
        risky_count = 0
        start = 0
        while True:
            index = lowered.find(phrase_lower, start)
            if index == -1:
                break
            context = lowered[max(0, index - 80) : index + len(phrase_lower) + 80]
            if not any(marker in context for marker in NEGATION_MARKERS):
                risky_count += 1
            start = index + len(phrase_lower)
        rows.append(
            {
                "check": f"risky_phrase:{phrase}",
                "status": "warn" if risky_count else "pass",
                "detail": f"occurrences={count}; risky_contexts={risky_count}",
            }
        )

    for path in required_files:
        rows.append(_check_file(path))
    return rows


def write_submission_qc_markdown(rows: list[dict[str, object]], out: Path) -> Path:
    out.parent.mkdir(parents=True, exist_ok=True)
    counts = {"pass": 0, "warn": 0, "fail": 0}
    for row in rows:
        counts[str(row["status"])] = counts.get(str(row["status"]), 0) + 1

    lines = [
        "# Submission QC Report",
        "",
        "This report checks whether the manuscript draft is internally consistent with",
        "current generated result tables and whether key submission files exist.",
        "",
        "## Summary",
        "",
        f"- pass: {counts.get('pass', 0)}",
        f"- warn: {counts.get('warn', 0)}",
        f"- fail: {counts.get('fail', 0)}",
        "",
        "## Checks",
        "",
        "| Check | Status | Detail |",
        "|---|---|---|",
    ]
    for row in rows:
        lines.append(f"| {row['check']} | {row['status']} | {row['detail']} |")
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out