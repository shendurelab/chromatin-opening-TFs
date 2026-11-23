# 60-TF multiome screen: processing scripts

Scripts that process the 60-TF multiome (RNA + ATAC) CRISPR screen from Cell Ranger ARC output to the final Signac object. Raw data and processed objects are available from GEO (GSE353547).

## Setup

- Edit the `<...>` placeholders for paths in `preprocess-s1.load-filter.R`, `preprocess-s2.call-tf-peaks.sh` and `preprocess-s3.aggr-tf-peaks.R`.
- Run all steps from the project directory, which must contain `data/`. All other paths are relative to it.
- Software: R 4.3.2, Seurat 5.1.0, Signac 1.13.0, monocle3 1.3.5.
- Cell Ranger ARC (`cellranger-arc aggr`) output is expected in `data/count/arc_merged/outs/`, and the gRNA assignment table (`gRNA-associations_4-umi_150-reads.tsv`, from GEO) in `data/gRNA_count/`.

## Run order

1. `preprocess-s1.load-filter.R` | Load and filter cells, add gRNA assignments and cell-cycle phase
2. `preprocess-s2.call-tf-peaks.sh` and `preprocess-s2.call-tf-peaks.R` | Call peaks for each TF's gRNA cells (MACS2 via Signac)
3. `preprocess-s3.aggr-tf-peaks.R` | Combine TF and NTC peaks, compute ATAC LSI/UMAP and neighbor graph
4. `preprocess-s4.process-RNA.R` | RNA dimension reduction and clustering, ATAC TF-IDF/SVD
