suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(patchwork))
suppressMessages(library(stringr))
suppressMessages(library(gridExtra))
suppressMessages(library(scales))
suppressMessages(library(reshape2))
suppressMessages(library(biomaRt))
suppressMessages(library(monocle3))
suppressMessages(library(tidyr))
suppressMessages(library(Signac))
suppressMessages(library(EnsDb.Mmusculus.v79))
suppressMessages(library(BSgenome.Mmusculus.UCSC.mm10))
suppressMessages(library(Matrix))
suppressMessages(library(magrittr))
suppressMessages(library(optparse))

option_list = list(
    make_option("--tf", type="character", default=NULL,
              help="a target transcription factor to call peaks for e.g. Gata6 gRNA cells", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
tf <- opt$tf

obj <- readRDS("./data/atac-rna.preprocessed-s1.signac.rds")

DefaultAssay(obj) <- "ATAC"

print(paste0("calling peaks for ",tf))
obj@meta.data$peak.group <- ifelse(grepl(paste0('(^|_)', tf, '(_|$)'), obj@meta.data$target) == TRUE, tf, 
                            ifelse(obj@meta.data$target == "NTC", "NTC", "other"))

Idents(object = obj) <- "peak.group"

obj <- subset(x = obj, peak.group != "other")

# Calling peaks on the seurat clusters
tf.peaks <- CallPeaks(
  object = obj,
  group.by = "peak.group",
  outdir = paste0("./data/peaks/",tf,"/"),
  name = tf,
  cleanup = FALSE
)

saveRDS(tf.peaks, file=paste0("./data/peaks/",tf,"/",tf,"_peaks.gr.rds"))

print("")

