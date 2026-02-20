# ==============================================================================
# Gap filling
# ==============================================================================
# Checking features
feat_def <- tibble::as_tibble(
  xcms::featureDefinitions(xchr7),
  rownames = "feature"
)
feat_val <- tibble::as_tibble(
  xcms::featureValues(xchr7, method = "sum"),
  rownames = "feature"
)

# Extract features with nas for peak filling
feat_with_na <- feat_val %>%
  tidyr::pivot_longer(cols = 2:ncol(.)) %>%
  dplyr::filter(is.na(value)) %>%
  dplyr::pull(feature) %>%
  unique(.)

# Filter feature defintions to features with nas
feat_def_nas <- feat_def %>%
  dplyr::filter(feature %in% feat_with_na)

# Create a list for checking the features chromatograms
feat_def_nas_vals <- feat_def_nas %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    feat_extract = rbind(
      c(
        mzmed - 0.0015,
        mzmed + 0.0015,
        rtmin - 2,
        rtmax + 2
      )
    )
  ) %>%
  dplyr::pull(feat_extract)

message(
  "Features with NAs prior to gap filling: ",
  sum(is.na(featureValues(xchr7)))
)

# Perform gap filling
if (check_saved("xchr8.rds")) {
  xchr8 <- readRDS(file = paste0(opt$output, "/objects/xchr8.rds"))
} else {
  xchr8 <- xcms::fillChromPeaks(
    object = xchr7,
    param = xcms::ChromPeakAreaParam()
  )
  saveRDS(object = xchr8, file = paste0(opt$output, "/objects/xchr8.rds"))
}

# Number of missing values after gap filling
message(
  "Features with NAs after gap filling: ",
  sum(is.na(xcms::featureValues(xchr8)))
)

# Number of filled peaks
message(
  "Number of filled peaks: ",
  sum(is.na(featureValues(xchr7))) - sum(is.na(xcms::featureValues(xchr8)))
)

# Plot all non-filled peaks
# Extract the m/z - rt regions for these features
# Extract features with nas for peak filling
feat_with_na_after <- tibble::as_tibble(
  xcms::featureValues(xchr8, method = "sum"),
  rownames = "feature"
) %>%
  tidyr::pivot_longer(cols = 2:ncol(.)) %>%
  dplyr::filter(is.na(value)) %>%
  dplyr::pull(feature) %>%
  unique(.)

chrs_na_feat <- xcms::featureArea(xchr8, features = feat_with_na_after)

# Expand the retention time by 1 second on both sides
chrs_na_feat[, "rtmin"] <- chrs_na_feat[, "rtmin"] - 1
chrs_na_feat[, "rtmax"] <- chrs_na_feat[, "rtmax"] + 1

if (check_saved("chrs_na.rds")) {
  chrs_na <- readRDS(file = paste0(opt$output, "/objects/chrs_na.rds"))
} else {
  chrs_na <- xcms::chromatogram(
    xchr8,
    mz = chrs_na_feat[, c("mzmin", "mzmax")],
    rt = chrs_na_feat[, c("rtmin", "rtmax")] # probably increase this a little
  )
  saveRDS(object = chrs_na, file = paste0(opt$output, "/objects/chrs_na.rds"))
}