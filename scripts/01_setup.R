
message(
  "===========================================================================",
  "\n",
  "Source dependendencies and load libraries ---------------------------------",
  "\n",
  "==========================================================================="
)

source("scripts/functions.R")
source("scripts/chem_functions.R")
source("scripts/met_disco_args.R")
suppressWarnings(
  suppressPackageStartupMessages({
    library(tidyverse)
    library(MSnbase)
    library(xcms)
    library(MsExperiment)
    library(RforMassSpectrometry)
    library(Spectra)
    library(Chromatograms)
    library(RColorBrewer)
    library(pheatmap)
    library(QFeatures)
    library(Rdisop)
    library(limma)
    library(BiocParallel)
    library(ggrepel)
    library(ComplexUpset)
    library(MetaboAnnotation)
    library(CompoundDb)
    library(MetaboCoreUtils)
    library(curl)
    library(gt)
    library(patchwork)
    library(AnnotationHub)
    library(optparse)

    # for examine_biotransformations_hits.R
    library(ComplexHeatmap)
    library(circlize)
    library(rcdk)
    library(rJava)
    library(fingerprint)
  })
)

message(
  "===========================================================================",
  "\n",
  "Create output folders -----------------------------------------------------",
  "\n",
  "==========================================================================="
)

folders <- c(
  "bpc",
  "internal_standard",
  "volcano",
  "upset",
  "feature_boxplot",
  "rtime",
  "filled_peaks",
  "pca",
  "feature_chromatogram_intensity",
  "per_sample_peaks",
  "features",
  "feature_pairs",
  "glycoside",
  "glycoside_feature_pairs"
)
dir.create(opt$output, FALSE, TRUE)
dir.create(file.path(opt$output, "objects"), FALSE, TRUE)
dir.create(file.path(opt$output, "tables"), FALSE, TRUE)
for (folder in folders) {
  dir.create(
    file.path(
      opt$output, "graphs", folder
    ),
    showWarnings = FALSE,
    recursive = TRUE
  )
}
dir.create(file.path("annotation_databases"), FALSE, TRUE)

# ==============================================================================
# Import metadata --------------------------------------------------------------
# ==============================================================================
message("Importing metadata...")
meta <- import_mzml(opt$data_path, opt$meta_file)
if (check_saved("ms_exp.rds")) {
  ms_exp <- readRDS(file = file.path(opt$output, "objects/ms_exp.rds"))
} else {
  ms_exp <- MsExperiment::readMsExperiment(
    spectraFiles = meta$path,
    sampleData = meta
  )
  saveRDS(object = ms_exp, file = file.path(opt$output, "objects/ms_exp.rds"))
}