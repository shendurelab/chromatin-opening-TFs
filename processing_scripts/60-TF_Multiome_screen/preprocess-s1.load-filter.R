# Loading in required packages
suppressMessages(library(dplyr))
suppressMessages(library(plyr))
suppressMessages(library(tidyr))
suppressMessages(library(patchwork))
suppressMessages(library(stringr))

suppressMessages(library(gridExtra))
suppressMessages(library(gridGraphics))
suppressMessages(library(grid))
suppressMessages(library(scales))
suppressMessages(library(reshape2))
suppressMessages(library(biomaRt))

suppressMessages(library(Seurat))
suppressMessages(library(monocle3))
suppressMessages(library(Signac))

suppressMessages(library(gprofiler2))
suppressMessages(library(EnsDb.Mmusculus.v79))
suppressMessages(library(BSgenome.Mmusculus.UCSC.mm10))

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "mm10"

# <path to project directory> must contain:
#   data/count/arc_merged/outs/        (cellranger-arc aggr output)
#   data/gRNA_count/gRNA-associations_4-umi_150-reads.tsv   (available on GEO)
# Run this script from the project directory; output is written to ./data/
raw_path <- "<path to project directory>"

#----------------------------
# Loading in the data
#----------------------------
message("Loading dataset as Seurat object")

path <- paste0(raw_path,"/data/count/arc_merged/outs/")
counts <- Read10X(paste0(path,"filtered_feature_bc_matrix/"))
fragpath <- paste0(path,"atac_fragments.tsv.gz")

# create a Seurat object containing the RNA adata
arc <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
arc[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

arc[["percent.mt"]] <- PercentageFeatureSet(arc, pattern = "^mt-")

#----------------------------
# Filtering low quality CBCs 
#----------------------------
message("Filtering low quality cells")

arc <- subset(
  x = arc,
  subset = nCount_ATAC < 50000 &
    nCount_RNA < 20000 &
    nCount_ATAC > 800 &
    nCount_RNA > 800 &
    percent.mt < 30
)

#----------------------------
# Adding the gRNA assignments 
#----------------------------
message("Adding gRNA information")

gRNA.path <- paste0(raw_path,'/data/gRNA_count/gRNA-associations_4-umi_150-reads.tsv')
gRNA <- read.table(gRNA.path, sep = "\t",header=T)

# formatting the gRNA CBC to match the
# formatting of the meta data CBCs
replicates <- str_split_fixed(gRNA$CBC, "_", 2)
replicate <- str_split_fixed(replicates[,1], "R", 2)[,2]
cbc <- str_split_fixed(replicates[,2], "-", 2)[,1]
gRNA$CBC <- paste0(cbc,"-",replicate)

# adding a column to the meta data with the 
# cell barcode
arc@meta.data$CBC <- rownames(arc@meta.data)

# adding the target assignment and gRNA 
# assignment to the meta data by matching
# the gRNA CBC to the meta data CBC
arc@meta.data$target <- gRNA$target[match(arc@meta.data$CBC, gRNA$CBC)]
arc@meta.data$gRNA <- gRNA$gRNAs[match(arc@meta.data$CBC, gRNA$CBC)]

# adding a column to the meta data indicating 
# replicae number
arc@meta.data$replicate <- str_split_fixed(arc@meta.data$CBC, "-", 2)[,2]

#----------------------------
# Assigning cell-cycle phase 
#----------------------------
message("Adding cell-cycle phase")

DefaultAssay(arc) <- 'RNA'

seurat.cell.cycle <- arc 
seurat.cell.cycle <- NormalizeData(seurat.cell.cycle)
seurat.cell.cycle <- FindVariableFeatures(seurat.cell.cycle, selection.method = "vst")
seurat.cell.cycle <- ScaleData(seurat.cell.cycle, features = rownames(seurat.cell.cycle))

# Getting S phase and G1 marker genes
m.s.genes = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
m.g2m.genes = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name

# Performing the cell-cycle scoring with Seurat
seurat.cell.cycle <- CellCycleScoring(seurat.cell.cycle, s.features = m.s.genes,
                                      g2m.features = m.g2m.genes, set.ident = TRUE)

# Creating a data frame that has cell barcodes and their corresponing cell-cycle phase assignment
cell.cycle.df <- seurat.cell.cycle@meta.data[, 'Phase', drop=FALSE]

# Adding the cell-cycle phase assignments to the original Seurat object meta data
arc@meta.data$phase <- cell.cycle.df$Phase[match(row.names(arc@meta.data),
                                                           row.names(cell.cycle.df))]

# Lastly, adding columns to the metadata that have the log of the RNA and ATAC UMI counts
arc@meta.data$log_nCount_RNA <- log2(arc@meta.data$nCount_RNA)
arc@meta.data$log_nCount_ATAC <- log2(arc@meta.data$nCount_ATAC)

#----------------------------
# Saving the Signac object
#----------------------------
saveRDS(arc, file = "./data/atac-rna.preprocessed-s1.signac.rds")
message("Initial filtering complete")

sessionInfo()
