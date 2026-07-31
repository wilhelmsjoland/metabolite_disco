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
    library(RSQLite) # for sql backend
    library(MsBackendSql) # for sql backend
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
# Snakesave -------------------------------
# ==============================================================================
saveRDS(
  object = list(xchr9 = xchr9),
  file = snakemake@output[[1]]
)

script_footer()
end_log()