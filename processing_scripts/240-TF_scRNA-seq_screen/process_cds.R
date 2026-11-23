# Loading in required packages
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(patchwork))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(scales))
suppressMessages(library(reshape2))
suppressMessages(library(biomaRt))
suppressMessages(library(monocle3))
suppressMessages(library(tidyr))
suppressMessages(library(optparse))

#-----------------------------------------------------
# Setting the options for the script (so that it can  
# easily be run with or without performing cell-cycle 
# regression)                                         
#-----------------------------------------------------

option_list = list(
    make_option("--cellcycle_regression", type="character", default=NULL,
              help="TRUE or FALSE of whether to regress out the effects of cell-cycle", metavar="character"),
    make_option("--out", type="character", default=NULL,
              help="name and path to write the output .rds file to", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Getting all of the different variables from the command line arguments
out.p <- opt$out
cc_regress <- opt$cellcycle_regression

#-----------------------------------------------------
# PREPROCESSING WITH SEURAT 
#-----------------------------------------------------
# - Adding the gRNA information
# - Adding predicted doublet scores from Scrublet,
#   and removing predicted doublets
# - Assigning cell-cycle phase so it can later be 
#   subtracted out with Monocle
#-----------------------------------------------------

# A function to load in each of the different 10X replicate datasets as a Seurat object
load_data <- function(replicate) {
    path <- paste('./data/RNA_counts/RNA_', replicate, '/outs/filtered_feature_bc_matrix', sep='')
    rep <- Read10X(data.dir = path)
    rep <- CreateSeuratObject(counts = rep, min.cells = 3, min.genes = 200, project = "10X")
    rep
}

# A function to add the gRNA assignments to the meta data
add_gRNA <- function(seurat, gRNA) {
    seurat@meta.data$CBC = rownames(seurat@meta.data)
    df1 <- left_join(seurat@meta.data, gRNA) 
    rownames(df1)=seurat$CBC
    seurat@meta.data = df1 #This Seurat object has the gRNA information
    seurat@meta.data
}

# Loading in all 12 replicates
seurat_R1 <- load_data('R1')
seurat_R2 <- load_data('R2')
seurat_R3 <- load_data('R3')
seurat_R4 <- load_data('R4')
seurat_R5 <- load_data('R5')
seurat_R6 <- load_data('R6')
seurat_R7 <- load_data('R7')
seurat_R8 <- load_data('R8')
seurat_R9 <- load_data('R9')
seurat_R10 <- load_data('R10')
seurat_R11 <- load_data('R11')
seurat_R12 <- load_data('R12')

# Merging all of the replicates into  a single Seurat object
replicates <- c(seurat_R2, seurat_R3, seurat_R4, seurat_R5, seurat_R6, seurat_R7, seurat_R8, seurat_R9, seurat_R10, seurat_R11, seurat_R12)
replicate_ids <- c('R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9', 'R10', 'R11', 'R12')
seurat.merged <- merge(seurat_R1, y = replicates, add.cell.ids = replicate_ids, project = "pioneer_screen")

# Adding a column to the meta data that has the replicate number of each cells
cbcs <- rownames(seurat.merged@meta.data)
replicates <- str_split_fixed(cbcs, "_", 2)
colnames(replicates) <- c("rep","barcode")
replicates <- as.data.frame(replicates)
replicate <- replicates$rep
seurat.merged@meta.data$replicate <- replicate

# Loading in the gRNA associations and adding it to the Seurat object meta data
gRNA <- read.table("./data/gRNA_counts/gRNA-associations.tsv", sep = "\t",header=T)
seurat.merged@meta.data = add_gRNA(seurat.merged, gRNA)

# A function to read in the doublet scores for each CBC that were predicted with Scrublet/scanpy
load_doublet <- function(replicate) {
    path <- paste("./data/RNA_counts/RNA_", replicate, "/outs/scrublet_metadata.csv", sep='')
    scrub.df <- read.csv(path, header = TRUE)
    prefix = paste0(replicate, '_')
    scrub.df$X = paste0(prefix, scrub.df$X)
    scrub.df
}

# Reading in all of the doublet score files for each replicate
R1.doublet = load_doublet('R1')
R2.doublet = load_doublet('R2')
R3.doublet = load_doublet('R3')
R4.doublet = load_doublet('R4')
R5.doublet = load_doublet('R5')
R6.doublet = load_doublet('R6')
R7.doublet = load_doublet('R7')
R8.doublet = load_doublet('R8')
R9.doublet = load_doublet('R9')
R10.doublet = load_doublet('R10')
R11.doublet = load_doublet('R11')
R12.doublet = load_doublet('R12')

# Combining all of the doublet scores for each replicate into a single object
doublet.merged <- rbind(R1.doublet, R2.doublet, R3.doublet, R4.doublet, R5.doublet, R6.doublet,
                       R7.doublet, R8.doublet, R9.doublet, R10.doublet, R11.doublet, R12.doublet)

# Adding the doublet scores to their respective cell barcodes in the Seurat objects meta data
seurat.merged@meta.data$doublet_score <- doublet.merged$doublet_score[match(seurat.merged@meta.data$CBC, doublet.merged$X)]

# Subsetting the Seurat object by removing cells/rows that have a doublet score greater than .2
seurat.merged <- subset(seurat.merged, subset = doublet_score < .2)

# The data need to be normalized and scaled for the cell-cycle scoring:
# First creating a new Seurat object to normalize and perform the cell-cycle scoring on and then adding the final phase 
# assignments to the the meta data of the original (unnormalized) Seurat object which we'll subsequently convert to a Monocle cds.
seurat.cell.cycle <- seurat.merged
seurat.cell.cycle <- NormalizeData(seurat.cell.cycle)
seurat.cell.cycle <- FindVariableFeatures(seurat.cell.cycle, selection.method = "vst")
seurat.cell.cycle <- ScaleData(seurat.cell.cycle, features = rownames(seurat.cell.cycle))

# A function to retrieve cell-cycle phase marker genes:
convertHumanGeneList <- function(x){
require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",host = "uswest.ensembl.org")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host = "uswest.ensembl.org")
genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
humanx <- unique(genesV2[, 2])
return(humanx)
}

# Getting S phase and G1 marker genes
m.s.genes <- convertHumanGeneList(cc.genes.updated.2019$s.genes)
m.g2m.genes <- convertHumanGeneList(cc.genes.updated.2019$g2m.genes)

# Performing the cell-cycle scoring with Seurat
seurat.cell.cycle <- CellCycleScoring(seurat.cell.cycle, s.features = m.s.genes, 
                                      g2m.features = m.g2m.genes, set.ident = TRUE)

# Creating a data frame that has cell barcodes and their corresponing cell-cycle phase assignment
cell.cycle.df <- seurat.cell.cycle@meta.data

# Adding the cell-cycle phase scores and assignments to the original Seurat object meta data
seurat.merged@meta.data$Phase <- cell.cycle.df$Phase[match(row.names(seurat.merged@meta.data), 
                                                           row.names(cell.cycle.df))]
seurat.merged@meta.data$S.Score <- cell.cycle.df$S.Score[match(row.names(seurat.merged@meta.data),
                                                           row.names(cell.cycle.df))]
seurat.merged@meta.data$G2M.Score <- cell.cycle.df$G2M.Score[match(row.names(seurat.merged@meta.data),
                                                           row.names(cell.cycle.df))]

# Lastly, adding a column to the meta data that has the percentage of mitochondrial reads for each cell
seurat.merged[["percent.mt"]] <- PercentageFeatureSet(seurat.merged, pattern = "^mt-")

seurat.merged@meta.data$log_nCount_RNA <- log(seurat.merged@meta.data$nCount_RNA)

#-----------------------------------------------------
# PREPROCESSING WITH MONOCLE3 
#-----------------------------------------------------
# - Creating a CDS object
# - Regressing out the impact of UMI and cell-cycle
#   phase with Monocle3
# - Performing dimensionality reduction and
#   clustering, and saving the Monocle object for
#   downstream analyses
#-----------------------------------------------------

# Extract count data, phenotype data, and feature data from the Seurat Object.
counts.data <- as(as.matrix(seurat.merged@assays$RNA@data), 'sparseMatrix')
pheno.data <- seurat.merged@meta.data
feature.data <- data.frame(gene_short_name = row.names(counts.data), row.names = row.names(counts.data))

# Creating a Monocle3 cds object from the count, phenotype, and feature data
cds <- new_cell_data_set(expression_data = counts.data, cell_metadata = pheno.data, gene_metadata = feature.data)

# Removing cells with greater than 10 percent mitochondrial reads
valid_cells <- row.names(subset(pData(cds),
            percent.mt < 10))
cds <- cds[,valid_cells]

# Removing cells with fewer than 1000 UMIs 
valid_cells <- row.names(subset(pData(cds),
            nCount_RNA >= 1000))
cds <- cds[,valid_cells]

# Preprocessing the data with Monocl3 
cds <- preprocess_cds(cds, num_dim = 100)

# Regressing out the impact of UMI count, cell-cycle phase, and replicate
if (cc_regress == "FALSE") {
            print("forgoing cell-cycle phase regression")
            cds <- align_cds(cds,
                   alignment_group = "replicate",
                   residual_model_formula_str = "~log_nCount_RNA")
} else {
            print("running cell-cycle phase regression")
            cds <- align_cds(cds,
                   alignment_group = "replicate",
                   residual_model_formula_str = "~G2M.Score + S.Score + log_nCount_RNA")
}

# Dimensionality reduction and clustering cells
cds <- reduce_dimension(cds,
                        umap.n_neighbors = 30,
                        umap.min_dist = 0.1,
                        reduction_method = "UMAP",
                        preprocess_method = "Aligned",
                        max_components = 3)

cds <- cluster_cells(cds)

# Saving the cds as an .rds file to easily be loaded in later
if (cc_regress == "FALSE") {
            print(paste0("writing the cds object WITHOUT cell-cycle regression to: ", out.p))
} else {
            print(paste0("writing the cds object WITH cell-cycle regression to: ", out.p)) 
}

saveRDS(cds, file = out.p)

