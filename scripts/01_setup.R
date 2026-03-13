# ==============================================================================
# Source functions and minimal startup parameters ------------------------------
# ==============================================================================
source("scripts/functions.R")
dir.create(snakemake@params$output, FALSE, TRUE)
start_log(snakemake@params$output)

library(this.path)
library(cli)
set.seed(snakemake@params$seed)
register_parallel(snakemake@params$cores)
results_path <- snakemake@params$output
config <- snakemake@config

# ==============================================================================
# Information on pipeline ------------------------------------------------------
# ==============================================================================
start_pipeline_msg()
script_header()
cli::cli_text("Running pipeline with: ")
purrr::walk2(
  .x = config,
  .y = names(config),
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
    library(BiocParallel)
    library(dplyr)
    library(purrr)
    library(readr)
    library(MsExperiment)
    library(RColorBrewer)
    library(tibble)
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
  # "feature_boxplot",
  "quality_control",
  "rtime",
  # "filled_peaks",
  "pca",
  # "feature_chromatogram_intensity",
  # "per_sample_peaks",
  # "features",
  # "feature_pairs",
  # "glycoside",
  # "glycoside_feature_pairs"
)
dir.create(file.path(snakemake@params$output, "objects"), FALSE, TRUE)
dir.create(file.path(snakemake@params$output, "snakemake_objects"), FALSE, TRUE)
dir.create(file.path(snakemake@params$output, "tables"), FALSE, TRUE)
for (folder in folders) {
  dir.create(
    file.path(
      snakemake@params$output, "graphs", folder
    ),
    showWarnings = FALSE,
    recursive = TRUE
  )
}
# dir.create(file.path("annotation_databases"), FALSE, TRUE)
cli::cli_progress_done()
# ==============================================================================
# Import metadata --------------------------------------------------------------
# ==============================================================================
meta <- import_mzml(snakemake@params$data_path, snakemake@params$meta_file)
meta_path <- file.path(
  snakemake@params$output,
  "tables",
  "metadata.csv"
)
readr::write_csv(meta, meta_path)
cli::cli_alert_success("Saved metadata to {.path {meta_path}}")

ms_exp_path <- file.path(snakemake@params$output, "objects", "ms_exp.rds")
if (interactive() && file.exists(ms_exp_path)) {
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

saveRDS(
  object = list(
    ms_exp = ms_exp,
    group_colors = group_colors,
    meta = meta,
    results_path = results_path
  ),
  file = snakemake@output[[1]]
)
end_log()