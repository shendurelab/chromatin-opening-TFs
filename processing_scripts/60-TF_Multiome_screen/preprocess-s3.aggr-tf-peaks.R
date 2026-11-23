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

###########################
# loading in existing
# Signac object with the
# peaks called on entire
# dataset
###########################

# reading in data
obj <- readRDS("./data/atac-rna.preprocessed-s1.signac.rds")

# defining path to the fragments file ouputed by cellranger 
raw_path <- "<path to project directory>"
path <- paste0(raw_path,"/data/count/arc_merged/outs/")
frag.path <- paste0(path,"atac_fragments.tsv.gz")

# converting NA values to whitespace in
# target and gRNA metadata columns
obj@meta.data$target[is.na(obj@meta.data$target)] <- ""
obj@meta.data$gRNA[is.na(obj@meta.data$gRNA)] <- ""

# adding a column with the number of gRNAs
# assigned to each cell
meta <- as.data.frame(obj@meta.data)
meta$n.guides <- str_count(meta$gRNA, "_") + 1
meta$n.guides[is.na(meta$n.guides)] <- 0

# Reading in a list of our 60 targets
targets <- read.csv("./data/targets.tsv",sep="\t",header=FALSE)$V1

############################
# aggregating peaks that
# were called on each TF
############################

# initializing a list to store
# peak GRange objects in for
# each TF
grlist <- list()

# iterating through each TF
for (tf in targets) {
    # reading the peaks called on that TF as a GRanges object 
    tf.gr <- readRDS(paste0("./data/peaks/",tf,"/",tf,"_peaks.gr.rds"))

    # getting a list of standard chromosomes in mouse (19 autosomes
    # and X and Y sex chromosomes
    valid_chr <- grep("^chr[0-9XY]+$", unique(as.character(seqnames(tf.gr))), value = TRUE)

    # subsetting the GRanges object to standard chromosomes
    # only
    tf.gr <- tf.gr[mcols(tf.gr)$peak_called_in == tf & 
                   as.character(seqnames(tf.gr)) %in% valid_chr]
    message("Number of ", tf, " specific peaks: ", length(tf.gr))
    # adding the peaks to the list of
    # GRanges objects	   
    grlist[[tf]] <- tf.gr
}

############################
# calling peaks in NTC cells 
# and cells that weren't
# assigned a gRNA
############################
   
print("Calling peaks on cells not assigned a gRNA as well as NTC cells...")

DefaultAssay(obj) <- "ATAC"

# adding a column to object metadata indicating
# whether a cell is a NTC cell, a cell not 
# assigned a gRNA, or other cells (cells assigned
# a targetting gRNA)
obj@meta.data$peak.group <- ifelse(obj@meta.data$target == "NTC", "NTC", 
                            ifelse(obj@meta.data$target == "", "no_guide", "other"))

# setting the identify of the object
Idents(object = obj) <- "peak.group"

# subsetting to cell not assigned a gRNA
# and to NTC cells only
sub.obj <- subset(x = obj, peak.group != "other")

# Calling peaks on the cells with CallPeaks 
ntc.peaks <- CallPeaks(
  object = sub.obj,
  group.by = "peak.group",
  outdir = paste0("./data/peaks/NTC/"),
  name = "NTC",
  cleanup = FALSE
)

# saving the GRanges object
saveRDS(ntc.peaks, file=paste0("./data/peaks/NTC/no-gRNA_peaks.gr.rds"))

# clearing some unneccessary variables
obj@meta.data$peak.group <- NULL
sub.obj <- NULL

############################
# creating a new single cell
# matrix using the peaks 
# that were called on 
# individual TFs as well as 
# NTC cells and on cells
# that weren't assigned a
# gRNA
############################

# getting a list of standard chromosomes in mouse (19 autosomes
# and X and Y sex chromosomes
valid_chr <- grep("^chr[0-9XY]+$", unique(as.character(seqnames(ntc.peaks))), value = TRUE)

# subsetting the GRanges object to standard chromosomes
# only
ntc.peaks <- ntc.peaks[as.character(seqnames(ntc.peaks)) %in% valid_chr]

# adding the NTC and unassigned peaks to our list
grlist[["unassigned"]] <- ntc.peaks

# reducing the GRanges list into a single object
gr.combined <- Reduce(f = c, x = grlist)
# reducing peaks that overlap into a single peak
gr <- reduce(x = gr.combined, with.revmap = TRUE)

# removing nonstandard chromosomes and blacklisted regions
gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
gr <- subsetByOverlaps(x = gr, ranges = blacklist_mm10, invert = TRUE)

# creating the new feature matrix
mat <- FeatureMatrix(
       fragments = Fragments(obj),
       features = gr,
       cells = colnames(obj)
       )

###########################
# creating Monocle3 object
# + and running dimensionality
# reduction with Monocle3
###########################
print("Creating a Monocle3 object and running dimensionality reduction...")
# getting metadata
pheno.data <- obj@meta.data

# binarizing the atac data
mat@x[mat@x > 0] <- 1
feature.data <- data.frame(gene_short_name = row.names(mat), row.names = row.names(mat))

# creating the Monocle3 object
peaks.cds <- new_cell_data_set(expression_data = mat,
                              cell_metadata = pheno.data, gene_metadata = feature.data)

peaks.cds <- monocle3::detect_genes(peaks.cds)
# ensuring there are no peaks
# with zero reads
peaks.cds <- peaks.cds[Matrix::rowSums(exprs(peaks.cds)) != 0,]

# Preprocessing the data with Monocle3 
peaks.cds <- detect_genes(peaks.cds)
peaks.cds <- estimate_size_factors(peaks.cds)
peaks.cds <- preprocess_cds(peaks.cds, method = "LSI")

peaks.cds <- reduce_dimension(peaks.cds, reduction_method = 'UMAP',
                              max_components = 2,
                              preprocess_method = "LSI")

peaks.cds <- cluster_cells(peaks.cds)
saveRDS(clusters(peaks.cds), file = "./data/monocle3_clusters.rds")

###########################
# converting the Monocle3
# object back to an assay
# in our Signac object
###########################
print("Converting to a Signac object...")

# getting the mm10 genome annotation
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"

# getting the count matrix from the  
counts <- exprs(peaks.cds)

# create a chromatin assay and add it to the 
# object as "tf-back.peaks"
obj[["tf.all.peaks"]] <- CreateChromatinAssay(
 counts = counts,
 sep = c(":", "-"),
 fragments = frag.path,
 annotation = annotation
)

# setting the default assay
DefaultAssay(object = obj) <- "tf.all.peaks"


###########################
# adding the atac PCA + UMAP 
###########################

# pca matrix
pca.mat <- reducedDims(peaks.cds)$LSI
# umap matrix
umap.mat <- reducedDims(peaks.cds)$UMAP

# adding the pca
obj[['TFallpca']] <- CreateDimReducObject(embeddings = pca.mat, key = 'TFallpca_', assay = DefaultAssay(obj))
# adding the umap
obj[['TFallumap']] <- CreateDimReducObject(embeddings = umap.mat, key = 'TFallumap_', assay = DefaultAssay(obj))

# running Seurat's L2Dim and FindNeighbors
print("running L2Dim and FindNeighbors")
kadj = round(0.5 * sqrt(ncol(obj)))
obj = obj %>%Seurat::L2Dim(reduction='TFallpca')%>%
Seurat::FindNeighbors(reduction='TFallpca.l2', dims=1:50, k=kadj)

saveRDS(obj, file = "./data/atac-rna.preprocessed-s3.signac.rds")

sessionInfo()
