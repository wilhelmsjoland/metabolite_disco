# ==============================================================================
# Filtering features and input to SummarizedExperiment -------------------------
# ==============================================================================
message("Filtering features based on missingness...")
group_factor <- MsExperiment::sampleData(xchr8)$grRoup
group_factor <- as.factor(group_factor)
if (check_saved("xchr9.rds")) {
  xchr9 <- readRDS(file = paste0(opt$output, "/objects/xchr9.rds"))
} else {
  xchr9 <- QFeatures::filterFeatures(
    xchr8,
    xcms::PercentMissingFilter(
      threshold = opt$missingness,
      f = group_factor
    )
  )
  saveRDS(object = xchr9, file = paste0(opt$output, "/objects/xchr9.rds"))
}

# TODO Fix - this was moved from later in the pipeline
# ==============================================================================
# Filtering features -----------------------------------------------------------
# ==============================================================================
message(sprintf(
  "Filtering features with sn: %s, beta_cor: %s, beta_snr: %s",
  opt$sn_threshold,
  opt$beta_cor_threshold,
  opt$beta_snr_threshold
))

if (file.exists(file.path(opt$output, "objects", "xchr9_filt.rds"))) {
  xchr9_filt <- readRDS(file.path(opt$output, "objects", "xchr9_filt.rds"))
} else {
  xchr9_filt <- filt_features(
    object = xchr9,
    sn_threshold = opt$sn_threshold,
    beta_cor_threshold = opt$beta_cor_threshold,
    beta_snr_threshold = opt$beta_snr_threshold
  )
  saveRDS(
    object = xchr9_filt,
    file = paste0(opt$output, "/objects/xchr9_filt.rds")
  )
}
