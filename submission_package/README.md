# Submission package

This directory is for locally assembled journal, preprint, or data-repository
materials.

Large generated files copied into this directory are intentionally ignored by
Git. Only this README is committed.

Use the packaging script from the repository root:

```powershell
.\scripts\prepare_submission_package.ps1
```

The script creates:

- `submission_package/manuscript/`
- `submission_package/figures/`
- `submission_package/supplementary_tables/`
- `submission_package/references/`
- `submission_package/manifest.tsv`

Before journal submission, upload the large supplementary tables and figures to
the journal system or to a data repository such as Zenodo, Figshare, or an
institutional archive.