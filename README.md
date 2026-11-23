# High-throughput identification of chromatin-opening transcription factors

> **Note:** This repository is under active construction. Code, data, and documentation are currently being organized and will be added incrementally. Scripts are shared for transparency and do not constitute streamlined pipelines.

## Repository layout

```text
├── figures/                         # Notebooks for generating manuscript figures
│   ├── Figure1.ipynb
│   ├── Figure2.ipynb
│   ├── Figure3.ipynb
│   └── load_packages.R              # Packages and helper functions used by the notebooks
├── figure_data/                     # Small processed files read by the notebooks (see its README)
└── processing_scripts/              # Data processing and analysis scripts
    ├── 240-TF_scRNA-seq_screen/     # Initial 240-TF scRNA-seq screen
    └── 60-TF_Multiome_screen/       # Follow-up 60-TF single-cell multiome screen
```

## What's in each directory

- **`figures/`** – Jupyter notebooks for generating the manuscript's main and supplemental figure panels.
- **`figure_data/`** – Small processed files that let parts of the notebooks run without the full datasets. Each file is described in `figure_data/README.md`.
- **`processing_scripts/240-TF_scRNA-seq_screen/`** – Scripts for processing and analyzing data from the initial pooled CRISPRa screen of 240 transcription factors using single-cell RNA-seq.
- **`processing_scripts/60-TF_Multiome_screen/`** – Scripts for processing and analyzing data from the follow-up pooled CRISPRa screen of 60 transcription factors using joint single-cell RNA-seq and ATAC-seq.

## Data availability

Raw and processed data are available from GEO (GSE353547).
