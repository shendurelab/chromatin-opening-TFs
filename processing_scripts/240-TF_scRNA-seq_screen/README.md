# 240-TF scRNA-seq screen: processing scripts

Raw data and processed objects are available from GEO (GSE353547).

The Monocle3 object used for the analyses in the manuscript was cell-cycle regressed, generated with:

    Rscript process_cds.R --cellcycle_regression TRUE --out ./data/monocle_cds.rds
