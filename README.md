# Sezary syndrome transcript detection analysis

Analysis code supporting the research letter "Type 2-associated transcriptional
features of an expanded T-cell clonotype in Sezary syndrome."

The analysis covers CD45-enriched skin-biopsy libraries from two patients with
Sezary syndrome in GEO GSE263881. HLA08 contributed two technical reactions
from one specimen at T3; HLA19 contributed one library at T5.

## Reconstruct the analysis from public GEO files

Python 3.10 or later is required.

```
python -m pip install -r requirements.txt
python ctcl_public_audit.py
python reproduce_submission.py --cells ctcl_audit_results/cell_marker_audit.csv.gz
```

The first analysis script downloads the public processed expression matrices,
feature annotations, barcodes, and filtered TCR contigs for the three SS
libraries. It applies the specified quality-control and CD3D selection rules
and matches qualifying TRA/TRB sequences by library and barcode. It covers
only the samples and marker genes needed for the report.

The summary script writes the numerical source data for Figures 1 and 2,
Supplementary Tables S1-S4, clonotype sequences, sensitivity summaries,
cohort totals, and runtime information to results/. It uses only the Python
standard library once the cell-level input is available.

## Run from the accompanying supplementary data

The journal companion file Sezary_Data.zip contains the selected cell-level
input and previously verified numerical outputs. Extract it, then supply its
input path to the summary script:

```
python reproduce_submission.py --cells /path/to/Sezary_Data/data/cell_marker_audit.csv.gz
```

The packaged input contains only the 7,929 selected SS cells and relevant
columns. The public reconstruction exports all cells from the three libraries;
the summary script applies the selection flag in either input. Input-file
hashes may therefore differ, while the eight summary CSVs must agree.

## Render the figures

After regenerating results/, run the following with Arial installed:

```
python figure_source/Figure_1_source/generate_figure_1.py --font-dir C:/Windows/Fonts
python figure_source/Figure_2_source/generate_figure_2.py --font-dir C:/Windows/Fonts
```

Use the local Arial directory on other systems, or omit --font-dir when Arial
is discoverable by matplotlib. Font files are not distributed.

## Public source accessions

GEO series: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263881

| Specimen | Expression | TCR |
|---|---|---|
| HLA19 T5 | GSM8200498 | GSM8200499 |
| HLA08 T3 reaction 1 | GSM8200500 | GSM8200501 |
| HLA08 T3 reaction 2 | GSM8200502 | GSM8200503 |

Source publication: https://doi.org/10.1182/blood.2024024817

data/source_manifest.json preserves the original recorded public-input hashes.
data/runtime_versions.json records the original public-data audit environment.
The supplementary data are maintained separately from this code package.

## Validation and scope

The standard-library summary script was run on the focused selected-cell input.
All eight numerical CSV outputs matched the previously verified submission
outputs exactly, and the rarefaction summaries matched Tables S3 and S4.
The optional fresh GEO download and artwork-rendering routes were not rerun
during this code-package preparation. Python source files were syntax checked.

All comparisons are descriptive. The HLA08 reactions are technical repeats,
not independent patient replicates. Malignant-cell identity was not
independently confirmed. Transcript detection does not establish protein
expression, functional signaling, a validated CD4 lineage, or a Th2A identity.
The code does not reconstruct the source study's Seurat/Harmony clustering.

## Associated manuscript authors

Julia A. Yescas, Joan Guitart, Xiaolong A. Zhou, and Ayman Grada.
