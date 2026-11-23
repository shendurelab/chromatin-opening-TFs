# figure_data

Small processed files read by the notebooks in `../figures/`. Large inputs (Seurat objects, per-TF DESeq2 results) are not included; see the top-level README.

| File | Description | Used in |
|---|---|---|
| 240-TF_MoransI_gRNA_K500.rds | gRNA-level Moran's I results, 240-TF scRNA-seq screen | Figure1.ipynb |
| targets.tsv | Top 60 TFs selected for the follow-up multiome screen | Figure1.ipynb |
| fc_out_raw.csv | Average log2FC in accessibility of each TF's DG sites in every TF's gRNA cells (60 x 61 matrix, includes non-DG column) | Figure3.ipynb |
| n_da_sites.csv | Number of DG sites per TF (plus number of non-DG peaks) | Figure3.ipynb |
| dg-motif-enrichment.monaLisa.rds | monaLisa motif enrichment in each TF's DG sites | Figure3.ipynb |
| chromVAR_overall_diffMotifs.csv | chromVAR motif accessibility differences, TF gRNA cells vs NTC | Figure3.ipynb |
| peak-umap.DA-sites.csv.gz | UMAP coordinates of DG/DL sites (one row per peak) | Figure3.ipynb |
| motif_ids.tsv | JASPAR motif ID to name lookup for the FIMO files | Figure3.ipynb |
| MA0465.3_fimo-min-p-value.tsv.gz | Per-peak minimum FIMO p-value, Cdx2 motif | Figure3.ipynb |
| MA0099.4_fimo-min-p-value.tsv.gz | Per-peak minimum FIMO p-value, FOS::JUN motif | Figure3.ipynb |
| MA0814.3_fimo-min-p-value.tsv.gz | Per-peak minimum FIMO p-value, Tfap2c motif | Figure3.ipynb |
| MA0142.1_fimo-min-p-value.tsv.gz | Per-peak minimum FIMO p-value, Pou5f1::Sox2 motif | Figure3.ipynb |

DG/DL sites for each TF are listed in Table S2.
