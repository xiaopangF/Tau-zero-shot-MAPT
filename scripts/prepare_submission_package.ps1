$ErrorActionPreference = "Stop"

$packageRoot = "submission_package"
$manuscriptDir = Join-Path $packageRoot "manuscript"
$figuresDir = Join-Path $packageRoot "figures"
$suppDir = Join-Path $packageRoot "supplementary_tables"
$refsDir = Join-Path $packageRoot "references"

New-Item -ItemType Directory -Force $manuscriptDir, $figuresDir, $suppDir, $refsDir | Out-Null

$items = @(
    @{ Source = "manuscript/manuscript_draft.md"; Target = "$manuscriptDir/manuscript_draft.md"; Type = "manuscript"; Description = "Combined manuscript draft" },
    @{ Source = "manuscript/figure_legends.md"; Target = "$manuscriptDir/figure_legends.md"; Type = "manuscript"; Description = "Figure legends draft" },
    @{ Source = "manuscript/data_code_availability.md"; Target = "$manuscriptDir/data_code_availability.md"; Type = "manuscript"; Description = "Data and code availability draft" },
    @{ Source = "manuscript/supplementary_tables.md"; Target = "$manuscriptDir/supplementary_tables.md"; Type = "manuscript"; Description = "Supplementary table descriptions" },
    @{ Source = "manuscript/submission_qc_report.md"; Target = "$manuscriptDir/submission_qc_report.md"; Type = "qc"; Description = "Submission QC report" },
    @{ Source = "manuscript/references.bib"; Target = "$refsDir/references.bib"; Type = "references"; Description = "BibTeX bibliography" },

    @{ Source = "results/manuscript_assets/figure_workflow_schematic.png"; Target = "$figuresDir/figure_1_workflow_schematic.png"; Type = "figure"; Description = "Figure 1 workflow schematic" },
    @{ Source = "results/figures/esm1v_ensemble/missense_heatmap.png"; Target = "$figuresDir/figure_2a_esm1v_heatmap.png"; Type = "figure"; Description = "Figure 2A ESM-1v heatmap" },
    @{ Source = "results/figures/esm1v_ensemble/score_by_position.png"; Target = "$figuresDir/figure_2b_score_by_position.png"; Type = "figure"; Description = "Figure 2B score by position" },
    @{ Source = "results/manuscript_assets/figure_domain_summary.png"; Target = "$figuresDir/figure_3_domain_summary.png"; Type = "figure"; Description = "Figure 3 domain summary" },
    @{ Source = "results/manuscript_assets/figure_model_comparison.png"; Target = "$figuresDir/figure_4a_model_comparison.png"; Type = "figure"; Description = "Figure 4A model comparison" },
    @{ Source = "results/manuscript_assets/figure_esm1v_vs_heuristic.png"; Target = "$figuresDir/figure_4b_esm1v_vs_heuristic.png"; Type = "figure"; Description = "Figure 4B ESM-1v vs heuristic" },
    @{ Source = "results/manuscript_assets/figure_vus_lollipop.png"; Target = "$figuresDir/figure_5_vus_lollipop.png"; Type = "figure"; Description = "Figure 5 VUS lollipop" },

    @{ Source = "data/processed/mapt_all_missense_variants.tsv"; Target = "$suppDir/supplementary_table_1_full_missense_atlas.tsv"; Type = "supplementary_table"; Description = "Full Tau-F missense atlas" },
    @{ Source = "results/scores/mapt_esm1v_ensemble.tsv"; Target = "$suppDir/supplementary_table_2_esm1v_ensemble_scores.tsv"; Type = "supplementary_table"; Description = "ESM-1v ensemble scores" },
    @{ Source = "data/processed/mapt_clinvar_benchmark.tsv"; Target = "$suppDir/supplementary_table_3a_clinvar_accepted.tsv"; Type = "supplementary_table"; Description = "ClinVar accepted rows" },
    @{ Source = "data/processed/mapt_clinvar_rejected.tsv"; Target = "$suppDir/supplementary_table_3b_clinvar_rejected.tsv"; Type = "supplementary_table"; Description = "ClinVar rejected rows" },
    @{ Source = "results/mapt_clinvar_qc_summary.tsv"; Target = "$suppDir/supplementary_table_3c_clinvar_qc_summary.tsv"; Type = "supplementary_table"; Description = "ClinVar QC summary" },
    @{ Source = "data/processed/mapt_alphamissense.tsv"; Target = "$suppDir/supplementary_table_4a_alphamissense_accepted.tsv"; Type = "supplementary_table"; Description = "AlphaMissense accepted rows" },
    @{ Source = "data/processed/mapt_alphamissense_rejected.tsv"; Target = "$suppDir/supplementary_table_4b_alphamissense_rejected.tsv"; Type = "supplementary_table"; Description = "AlphaMissense rejected rows" },
    @{ Source = "results/mapt_alphamissense_qc_summary.tsv"; Target = "$suppDir/supplementary_table_4c_alphamissense_qc_summary.tsv"; Type = "supplementary_table"; Description = "AlphaMissense QC summary" },
    @{ Source = "results/mapt_esm1v_top50_vus_priority.tsv"; Target = "$suppDir/supplementary_table_5_top_vus_priority.tsv"; Type = "supplementary_table"; Description = "Top VUS priority list" },
    @{ Source = "results/mapt_model_concordance.tsv"; Target = "$suppDir/supplementary_table_6a_model_concordance_all.tsv"; Type = "supplementary_table"; Description = "Model concordance all variants" },
    @{ Source = "results/mapt_model_concordance_summary.tsv"; Target = "$suppDir/supplementary_table_6b_model_concordance_summary.tsv"; Type = "supplementary_table"; Description = "Model concordance summary" },
    @{ Source = "results/mapt_model_concordance_top.tsv"; Target = "$suppDir/supplementary_table_6c_model_concordance_top.tsv"; Type = "supplementary_table"; Description = "Model concordance top rows" }
)

$manifestRows = New-Object System.Collections.Generic.List[object]
foreach ($item in $items) {
    if (-not (Test-Path $item.Source)) {
        throw "Required source file is missing: $($item.Source)"
    }
    Copy-Item -LiteralPath $item.Source -Destination $item.Target -Force
    $copied = Get-Item -LiteralPath $item.Target
    $manifestRows.Add([pscustomobject]@{
        type = $item.Type
        source = $item.Source
        packaged_file = $item.Target
        bytes = $copied.Length
        description = $item.Description
    })
}

$manifestPath = Join-Path $packageRoot "manifest.tsv"
$manifestRows | Export-Csv -LiteralPath $manifestPath -Delimiter "`t" -NoTypeInformation
Write-Host "Prepared submission package at $packageRoot"
Write-Host "Wrote manifest to $manifestPath"