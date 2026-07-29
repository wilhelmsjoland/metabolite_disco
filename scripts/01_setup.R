# ==============================================================================
# Source functions and minimal startup parameters ------------------------------
# ==============================================================================
source("scripts/functions.R")
dir.create(snakemake@params$output, FALSE, TRUE)
start_log(snakemake@params$output)

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
    library(BiocParallel)
    library(dplyr)
    library(purrr)
    library(readr)
    library(MsExperiment)
    library(MsBackendSql)
    library(RSQLite)
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
  "pca"
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
  # Materialize all spectra once here so downstream steps never re-read from
  # the mzML files. Repeated reads trigger an intermittent macOS-only
  # memory-corruption bug in mzR/proteowizard (sneumann/xcms#422).
  #
  # MsBackendMemory() (below, commented) works but embeds the full decoded
  # spectra data into the ms_exp object itself - since that object gets
  # saveRDS()'d again at nearly every later pipeline step, the payload was
  # getting duplicated into every checkpoint file (measured: ~0.86 GB per
  # checkpoint for this dataset, so 8-17+ GB per experiment across the
  # pipeline). MsBackendSql instead writes the data once to a SQLite file
  # on disk and keeps only a lightweight connection reference in the R
  # object, so every checkpoint after this one stays tiny (tens of KB)
  # regardless of how many steps re-save it (measured: ~1.17 GB total for
  # this dataset, one-time, vs. 1.72 GB after just 2 Memory-backed steps).
  # MsExperiment::spectra(ms_exp) <- retry_on_error(
  #   quote(
  #     Spectra::setBackend(
  #       MsExperiment::spectra(ms_exp),
  #       Spectra::MsBackendMemory()
  #     )
  #   )
  # )
  spectra_db_path <- file.path(
    snakemake@params$output,
    "objects",
    "spectra.sqlite"
  )
  MsExperiment::spectra(ms_exp) <- retry_on_error(
    quote({
      # a failed attempt can leave a partially-written db behind, which
      # collides with the next retry ("table already exists") - start clean
      if (file.exists(spectra_db_path)) file.remove(spectra_db_path)
      Spectra::setBackend(
        MsExperiment::spectra(ms_exp),
        MsBackendSql::MsBackendOfflineSql(),
        drv = RSQLite::SQLite(),
        dbname = spectra_db_path
      )
    })
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

script_footer()
end_log()