# ==============================================================================
# Source dependendencies and load libraries ------------------------------------
# ==============================================================================
library(this.path)
library(cli)

cli::cli_h1(basename(this.path::this.path()))
cli::cli_alert_info("Sourcing dependendencies")

source("scripts/functions.R")
source("scripts/chem_functions.R")
source("scripts/met_disco_args.R")

cli::cli_alert_info("Loading libraries")
suppressWarnings(
  suppressPackageStartupMessages({
    library(cli)
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

# ==============================================================================
# Create output folders --------------------------------------------------------
# ==============================================================================
cli::cli_alert_info("Creating output folders")
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

meta <- import_mzml(opt$data_path, opt$meta_file)
if (check_saved("ms_exp.rds")) {
  ms_exp <- readRDS(file = file.path(opt$output, "objects/ms_exp.rds"))
  cli::cli_alert_success(
    paste0(
      "Read saved .mzml files and metadata from ",
      "{.path {file.path(opt$output, 'objects/ms_exp.rds')}}"
    )
  )
} else {
  cli::cli_alert_info("Importing .mzml files and metadata")
  ms_exp <- MsExperiment::readMsExperiment(
    spectraFiles = meta$path,
    sampleData = meta
  )
  saveRDS(object = ms_exp, file = file.path(opt$output, "objects/ms_exp.rds"))
  cli::cli_alert_success(
    paste0(
      "Saved ms experiment to ",
      "{.path {file.path(opt$output, 'objects/ms_exp.rds')}}"
    )
  )
}

# ==============================================================================
# Set colors for groups --------------------------------------------------------
# ==============================================================================
cli::cli_alert_info("Setting colors for groups")
groups_to_use <- unique(MsExperiment::sampleData(ms_exp)$group)
group_colors <- paste0(
  RColorBrewer::brewer.pal(
    n = length(groups_to_use), "Set1"
  )[seq_along(groups_to_use)]
)
group_colors <- setNames(group_colors, groups_to_use)