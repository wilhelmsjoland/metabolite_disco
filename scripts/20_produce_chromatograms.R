# ==============================================================================
# Produce chromatograms --------------------------------------------------------
# ==============================================================================

message("Producing feature chromatograms...")
if (check_saved("feature_chrs.rds")) {
  feature_chrs <- readRDS(
    file = paste0(opt$output, "/objects/feature_chrs.rds")
  )
} else {
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
  saveRDS(
    object = feature_chrs,
    file = paste0(opt$output, "/objects/feature_chrs.rds")
  )
}

