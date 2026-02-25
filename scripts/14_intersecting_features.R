cli::cli_h1(basename(this.path::this.path()))
# ==============================================================================
# Finding intersecting features ------------------------------------------------
# ==============================================================================
cli::cli_progress_step("Extracting intersecting features")

# Generate all distinct upsets comparisons
upset_distinct_comps_path <- file.path(
  opt$output,
  "objects",
  "upset_distinct_comps.rds"
)
if (file.exists(upset_distinct_comps_path)) {
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
  opt$output,
  "objects",
  "upset_intersect_comps.rds"
)
if (file.exists(upset_intersect_comps_path)) {
  cli::cli_alert_info(
    "{.path {upset_intersect_comps_path}} already exists, skipping"
  )
  upset_intersect_comps <- readRDS(upset_intersect_comps_path)
} else {
  upset_intersect_comps <- extract_upset_comps(upset_distinct)
  saveRDS(upset_intersect_comps, upset_intersect_comps_path)
  cli::cli_alert_success(
    paste0(
      "Saved intersecting upset comparisons to: ",
      "{.path {upset_intersect_comps_path}}"
    )
  )
}
cli::cli_progress_done()

# Only temporary
# Find the exact index of the groups, regardless of order
int_upset_comp <- c(
  "bu_mutant_apiin-bu_wt_apiin",
  "bu_mutant_control-bu_wt_apiin",
  "bu_wt_apiin-bu_wt_control"
)
upset_intersect_id <- extract_upset_id(upset_distinct, int_upset_comp)
upset_comp <- upset_distinct_comps[[upset_intersect_id]]
