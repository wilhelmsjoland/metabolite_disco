# ==============================================================================
# Source functions and minimal startup parameters ------------------------------
# ==============================================================================
source("scripts/functions.R")
start_log(snakemake@params$output)
script_header()

set.seed(snakemake@params$seed)
register_parallel(snakemake@params$cores)
suppressWarnings(
  suppressPackageStartupMessages({
    library(cli)
    library(BiocParallel)
    library(xcms)
    library(MsExperiment)
    library(QFeatures)
    library(dplyr)
  })
)

gap_filling <- readRDS(snakemake@input[["gap_filling"]])
xchr8 <- gap_filling$xchr8


# ==============================================================================
# Filtering features and input to SummarizedExperiment -------------------------
# ==============================================================================
cli::cli_h3("Filtering features based on missingness")
group_factor <- as.factor(MsExperiment::sampleData(xchr8)$group)

xchr9_path <- file.path(
  snakemake@params$output,
  "objects",
  "xchr9.rds"
)
if (interactive() && file.exists(xchr9_path)) {
  xchr9 <- readRDS(file = xchr9_path)
  cli::cli_alert_success(
    paste0(
      "Imported missingness filtered object from: ",
      "{.path {xchr9_path}}"
    )
  )
} else {
  cli::cli_alert_info(
    "Filtering based on missingness of: {.val {snakemake@params$missingness}}"
  )
  xchr9 <- QFeatures::filterFeatures(
    xchr8,
    xcms::PercentMissingFilter(
      threshold = snakemake@params$missingness,
      f = group_factor
    )
  )
  saveRDS(object = xchr9, file = xchr9_path)
  cli::cli_alert_success(
    paste0(
      "Saved missingness filtered object to: ",
      "{.path {xchr9_path}}"
    )
  )
}

# ==============================================================================
# Filtering features with sn, beta_cor, beta_snr -------------------------------
# ==============================================================================
cli::cli_h3("Filtering features based on thresholds")

xchr9_filt_path <- file.path(
  snakemake@params$output,
  "objects",
  "xchr9_filt.rds"
)
if (interactive() && file.exists(xchr9_filt_path)) {
  xchr9_filt <- readRDS(xchr9_filt_path)
  cli::cli_alert_success(
    paste0(
      "Imported filtered feature object from: ",
      "{.path {xchr9_filt_path}}"
    )
  )
} else {
  cli::cli_bullets(
    c(
      "i" = "Filtering features with:",
      "*" = "sn: {.val {snakemake@params$sn_threshold}}",
      "*" = "beta_cor: {.val {snakemake@params$beta_cor_threshold}}",
      "*" = "beta_snr: {.val {snakemake@params$beta_snr_threshold}}"
    )
  )
  xchr9_filt <- filt_features(
    object = xchr9,
    sn_threshold = snakemake@params$sn_threshold,
    beta_cor_threshold = snakemake@params$beta_cor_threshold,
    beta_snr_threshold = snakemake@params$beta_snr_threshold
  )
  saveRDS(
    object = xchr9_filt,
    file = xchr9_filt_path
  )
  cli::cli_alert_success(
    paste0(
      "Saved filtered feature object to: ",
      "{.path {xchr9_filt_path}}"
    )
  )
}

# ==============================================================================
# Snakesave -------------------------------
# ==============================================================================
saveRDS(
  object = list(xchr9 = xchr9, xchr9_filt = xchr9_filt),
  file = snakemake@output[[1]]
)
end_log()