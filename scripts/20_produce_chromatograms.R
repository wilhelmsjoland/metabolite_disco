cli::cli_h1(basename(this.path::this.path()))
# ==============================================================================
# Generating feature chromatograms ---------------------------------------------
# ==============================================================================
feature_chrs_path <- file.path(
  opt$output,
  "objects",
  "feature_chrs.rds"
)

if (file.exists(feature_chrs_path)) {
  feature_chrs <- readRDS(file = feature_chrs_path)
  cli::cli_alert_success(
    paste0(
      "Imported feature chromatograms from: ",
      "{.path {feature_chrs_path}}"
    )
  )
} else {
  cli_progress_step("Generating feature chromatograms")
  feature_chrs <- xcms::featureChromatograms(
    object = xchr9,
    expandRt = 0,
    expandMz = 0,
    aggregationFun = "sum",
    filled = TRUE,
    features = all_sig_diff,
    missing = 0,
    return.type = "XChromatograms"
  )
  cli::cli_progress_done()
  saveRDS(
    object = feature_chrs,
    file = paste0(opt$output, "/objects/feature_chrs.rds")
  )
  cli::cli_alert_success(
    paste0(
      "Saved feature chromatograms to :",
      "{.path {feature_chrs_path}}"
    )
  )
}

# Only temporary
# Find the exact index of the groups, regardless of order
int_upset_comp <- c(
  "bu_mutant_apiin-bu_wt_apiin",
  "bu_mutant_control-bu_wt_apiin",
  "bu_wt_apiin-bu_wt_control"
)
upset_intersect_id <- extract_upset_id(upset_distinct, int_upset_comp)
upset_comp <- upset_distinct_comps[[upset_intersect_id]]