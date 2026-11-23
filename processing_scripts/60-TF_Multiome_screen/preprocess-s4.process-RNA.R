library(plyr)
library(dplyr)
library(Seurat)
library(patchwork)
library(stringr)
library(ggplot2)
library(gridExtra)
library(stringr)
library(scales)
library(reshape2)
library(biomaRt)
library(monocle3)
library(tidyr)
library(ggpubr)
library(viridis)
library(Signac)
suppressMessages(library(EnsDb.Mmusculus.v79))
suppressMessages(library(BSgenome.Mmusculus.UCSC.mm10))
suppressMessages(library(Matrix))
suppressMessages(library(magrittr))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))

# ---------------------------------------------------------------------------------
# Run RNA dimension reduction via Monocle3 using raw counts
# and finally add those reduced dimensions back into the Signac/
# Seurat object for further clustering in Seurat.
# ---------------------------------------------------------------------------------

# Load the existing Signac/Seurat object
obj <- readRDS("./data/atac-rna.preprocessed-s3.signac.rds")

# Convert to a Monocle3 cell_data_set, using the raw RNA counts
pheno.data <- obj@meta.data
DefaultAssay(obj) <- "RNA"

# Use raw counts (obj@assays$RNA@counts), not log-normalized data
RNA.counts.data <- GetAssayData(obj, assay="RNA", layer = "counts")

# Create a data frame of features for Monocle3
RNA.feature.data <- data.frame(
  gene_short_name = rownames(RNA.counts.data),
  row.names       = rownames(RNA.counts.data)
)

# Build the Monocle3 cell_data_set
RNA.cds <- new_cell_data_set(
  expression_data = RNA.counts.data,
  cell_metadata   = pheno.data,
  gene_metadata   = RNA.feature.data
)

# Preprocessing the RNA data with Monocle3
# num_dim=100 means Monocle3 will calculate up to 100 principal components
RNA.cds <- preprocess_cds(RNA.cds, num_dim = 100)

# Optionally regress out total RNA count (stored in meta.data as nCount_RNA)
# This helps align/normalize cells by library size
print("Regressing out UMI counts...")
RNA.cds <- align_cds(
  RNA.cds,
  residual_model_formula_str = '~log_nCount_RNA'
)

# Dimensionality reduction (UMAP by default) and clustering
RNA.cds <- reduce_dimension(RNA.cds)
RNA.cds <- cluster_cells(RNA.cds)

# Extract the Monocle3 PCA and UMAP embeddings
# Monocle3 automatically stores the PCA coordinates in reducedDims(RNA.cds)$PCA
# and the UMAP coordinates in reducedDims(RNA.cds)$UMAP
pca.mat  <- reducedDims(RNA.cds)$PCA
umap.mat <- reducedDims(RNA.cds)$UMAP

# Ensuring the default assay of the object is RNA
DefaultAssay(obj) <- "RNA"

# Add the Monocle3 PCA/UMAP back to the Signac/Seurat object as new DimReduc slots
obj[['pcaRNA']] <- CreateDimReducObject(
  embeddings = pca.mat,
  key        = 'pcaRNA_',
  assay      = DefaultAssay(obj)
)

obj[['umapRNA']] <- CreateDimReducObject(
  embeddings = umap.mat,
  key        = 'umapRNA_',
  assay      = DefaultAssay(obj)
)

# -------------------------------------------------------------------------------
# Normalizing & scaling the RNA data again in Seurat
# -------------------------------------------------------------------------------

DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(obj)
obj <- ScaleData(obj, features = all.genes)

# Ensure there are no NA values in gRNA column
obj@meta.data$gRNA <- ifelse(is.na(obj@meta.data$gRNA), "", obj@meta.data$gRNA)

# Use the Monocle3 PCA ("pcaRNA") as the reduction for Seurat's neighbor graph
obj <- FindNeighbors(obj, assay = "RNA", reduction = "pcaRNA", dims = 1:20)
obj <- FindClusters(obj, assay = "RNA", resolution = 1.2)


# -------------------------------------------------------------------------------
# TFIDF normalization of the ATAC matrix with peaks that were called for each TF
# -------------------------------------------------------------------------------

DefaultAssay(obj) <- "tf.all.peaks"
obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj, min.cutoff = 'q0')
obj <- RunSVD(obj)

# Save the updated object
saveRDS(obj, file = "./data/atac-rna.preprocessed-s4.signac.rds")

sessionInfo()
