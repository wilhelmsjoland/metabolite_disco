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
    library(purrr)
    library(ComplexHeatmap)
    library(dplyr)
  })
)

upset_data <- readRDS(snakemake@input[["upset"]])
upset_distinct <- upset_data$upset_distinct
upset_intersect <- upset_data$upset_intersect

# ==============================================================================
# Finding intersecting features ------------------------------------------------
# ==============================================================================
cli::cli_progress_step("Extracting intersecting features")

# Generate all distinct upsets comparisons
upset_distinct_comps_path <- file.path(
  snakemake@params$output,
  "objects",
  "upset_distinct_comps.rds"
)
if (interactive() && file.exists(upset_distinct_comps_path)) {
  cli::cli_alert_info(
    "{.path {upset_distinct_comps_path}} already exists, skipping"
  )
  upset_distinct_comps <- readRDS(upset_distinct_comps_path)
} else {
  upset_distinct_comps <- extract_upset_comps(upset_distinct)
  saveRDS(upset_distinct_comps, upset_distinct_comps_path)
  cli::cli_alert_success(
    paste0(
      "Saved distinc upset comparisons to: ",
      "{.path {upset_distinct_comps_path}}"
    )
  )
}

# Generate all intersecting upsets comparisons
upset_intersect_comps_path <- file.path(
  snakemake@params$output,
  "objects",
  "upset_intersect_comps.rds"
)
if (interactive() && file.exists(upset_intersect_comps_path)) {
  cli::cli_alert_info(
    "{.path {upset_intersect_comps_path}} already exists, skipping"
  )
  upset_intersect_comps <- readRDS(upset_intersect_comps_path)
} else {
  upset_intersect_comps <- extract_upset_comps(upset_intersect)
  saveRDS(upset_intersect_comps, upset_intersect_comps_path)
  cli::cli_alert_success(
    paste0(
      "Saved intersecting upset comparisons to: ",
      "{.path {upset_intersect_comps_path}}"
    )
  )
}
cli::cli_progress_done()

# ==============================================================================
# Snakesave --------------------------------------------------------------------
# ==============================================================================
saveRDS(
  object = list(
    upset_distinct_comps = upset_distinct_comps,
    upset_intersect_comps = upset_intersect_comps
  ),
  file = snakemake@output[[1]]
)
end_log()