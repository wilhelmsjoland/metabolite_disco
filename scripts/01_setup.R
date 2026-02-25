# ==============================================================================
# Source functions and minimal startup parameters ------------------------------
# ==============================================================================
source("scripts/functions.R")
source("scripts/chem_functions.R")
source("scripts/met_disco_args.R")
start_log()
library(this.path)
library(cli)

# ==============================================================================
# Information on pipeline ------------------------------------------------------
# ==============================================================================
start_pipeline_msg()
cli::cli_h1(basename(this.path::this.path()))
cli::cli_text("Running pipeline with: ")
purrr::walk2(
  .x = opt,
  .y = names(opt),
  .f = ~ cli::cli_bullets(c("i" = paste0(.y, ": {.val {.x}}")))
)
cli::cli_rule()

# ==============================================================================
# Loading libraries ------------------------------------------------------------
# ==============================================================================
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
    library(future.apply)

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
cli::cli_progress_step("Creating output folders")
folders <- c(
  "bpc",
  "internal_standard",
  "volcano",
  "upset",
  "feature_boxplot",
  "quality_control",
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
cli::cli_progress_done()
# ==============================================================================
# Import metadata --------------------------------------------------------------
# ==============================================================================
meta <- import_mzml(opt$data_path, opt$meta_file)
ms_exp_path <- file.path(opt$output, "objects", "ms_exp.rds")
if (file.exists(ms_exp_path)) {
  ms_exp <- readRDS(file = ms_exp_path)
  cli::cli_alert_success(
    paste0(
      "Imported saved .mzml files and metadata object from ",
      "{.path {ms_exp_path}}"
    )
  )
} else {
  cli::cli_alert_info("Importing .mzml files and metadata")
  ms_exp <- MsExperiment::readMsExperiment(
    spectraFiles = meta$path,
    sampleData = meta
  )
  saveRDS(object = ms_exp, file = ms_exp_path)
  cli::cli_alert_success(
    paste0(
      "Saved ms experiment to ",
      "{.path {ms_exp_path}}"
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