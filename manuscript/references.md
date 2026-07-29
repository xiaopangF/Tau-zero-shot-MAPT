# References Notes

This file explains the core references in plain language. The machine-readable
BibTeX entries are in `manuscript/references.bib`.

## Core Tau / MAPT references

### Hutton et al. 1998

Citation key: `hutton1998tau`

Why we cite it: this is one of the key papers linking mutations in the tau gene
(MAPT) to inherited frontotemporal dementia and parkinsonism linked to chromosome
17. It supports the statement that MAPT mutations can directly cause
neurodegenerative disease.

### Spillantini et al. 2000

Citation key: `spillantini2000tau`

Why we cite it: this review explains Tau gene mutations, FTDP-17, and Tau
isoforms. It supports our discussion that Tau biology and coordinates are
complicated because multiple isoforms exist.

## Variant database reference

### Landrum et al. 2018

Citation key: `landrum2018clinvar`

Why we cite it: this is the ClinVar database paper. We cite it when explaining
that ClinVar provides clinical variant interpretations, but that submitted labels
can be sparse or uncertain.

## Protein language model references

### Rives et al. 2021

Citation key: `rives2021biological`

Why we cite it: this paper supports the idea that protein language models learn
biological information from large protein sequence datasets.

### Meier et al. 2021

Citation key: `meier2021language`

Why we cite it: this is the ESM-1v / zero-shot mutation-effect prediction paper.
It directly supports our method of using pretrained protein language models
without Tau-specific training labels.

## External missense predictor reference

### Cheng et al. 2023

Citation key: `cheng2023alphamissense`

Why we cite it: this is the AlphaMissense paper. We cite it when explaining the
external AlphaMissense baseline and why coordinate QC is needed before comparing
its scores with the Tau-F atlas.

## How to use these references later

In Markdown or Pandoc-style writing, cite papers like this:

`Tau mutations can cause inherited tauopathy [@hutton1998tau].`

When converting to a manuscript, use `references.bib` as the bibliography file.