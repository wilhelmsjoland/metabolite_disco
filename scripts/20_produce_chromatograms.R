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
    # TODO
    # Fix by filtering
    features = xchr9_filt$final.plotting.features,
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
