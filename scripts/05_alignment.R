cli::cli_h1(basename(this.path::this.path()))

# ==============================================================================
# Grouping of peak groups
# ==============================================================================
cli::cli_h3("Performing first peak grouping")

xchr5_path <- file.path(opt$output, "objects", "xchr5.rds")
if (file.exists(xchr5_path)) {
  xchr5 <- readRDS(xchr5_path)
  cli::cli_alert_success(
    paste0(
      "Imported saved peak grouping object from: ",
      "{.path {xchr5_path}}"
    )
  )
} else {
  cli::cli_alert_info("First peak grouping")
  xchr5 <- xcms::groupChromPeaks(
    object = xchr, # xchr4
    param = xcms::PeakDensityParam(
      sampleGroups = MsExperiment::sampleData(xchr)$group,
      bw = opt$bw_first_grouping,
      minFraction = 0.5, # 0.7
      binSize = 0.01,
      maxFeatures = 1000, # 200
      ppm = opt$ppm_global,
      minSamples = 2 # 1
    )
  )
  saveRDS(object = xchr5, file = xchr5_path)
  cli::cli_alert_success(
    paste0(
      "Saved peak grouping object to: ",
      "{.path {xchr5_path}}"
    )
  )
}

# ==============================================================================
# Print first peak grouping params used to console -----------------------------
# ==============================================================================
xchr5_params <- xchr5@processHistory[[2]]@param

cli::cli_alert_success("First peak grouping performed with:")
purrr::walk(
  .x = slotNames(xchr5_params),
  .f = ~ cli::cli_bullets(
    c(
      "i" = paste0(
        .x, ": ",
        paste( # needed to not repeat the names twice for > 1 vectors
          slot(xchr5_params, .x),
          collapse = ", "
        )
      )
    )
  )
)
# TODO
# FROM HERE ====================================================================

# ==============================================================================
# Alignment of retention times
# ==============================================================================
# TODO Use this part to check if anchor peaks cover the whole RT
# Get the anchor peaks that would be selected
pgm <- xcms::adjustRtimePeakGroups(
  xchr5,
  xcms::PeakGroupsParam(minFraction = 0.9)
)

# Evaluate distribution of anchor peaks' rt in the first sample
# TODO compare this to the full range of retention times from the runs
chrom_peaks5 <- xcms::chromPeaks(xchr5)
max_rt_time <- max(chrom_peaks5[, "rtmax"], na.rm = TRUE)
min_rt_time <- min(chrom_peaks5[, "rtmin"], na.rm = TRUE)

# TODO Add a check here as well!
# quantile(pgm[, 1], na.rm = TRUE) # Remove na.rm = TRUE?
# End of checking for anchor peak rt distribution

cli::cli_h3("Aligning retention times")
# Alignment
if (check_saved("xchr6.rds")) {
  xchr6 <- readRDS(paste0(opt$output, "/objects/xchr6.rds"))
} else {
  xchr6 <- xcms::adjustRtime(
    object = xchr5,
    param = xcms::PeakGroupsParam(
      minFraction = 0.8, # 0.9
      extraPeaks = 0, # 1
      smooth = "loess",
      span = 0.6, # 0.6
      family = "gaussian",
      peakGroupsMatrix = matrix(nrow = 0, ncol = 0),
      subset = integer(),
      subsetAdjust = c("average", "previous")
    )
  )
  saveRDS(object = xchr6, file = paste0(opt$output, "/objects/xchr6.rds"))
}

# Checking for retention drift
# Extract base peak chromatograms
if (check_saved("bpc_after.rds")) {
  bpc_after <- readRDS(file = paste0(opt$output, "/objects/bpc_after.rds"))
} else {
  bpc_after <- xcms::chromatogram(
    xchr6,
    aggregationFun = "max",
    chromPeaks = "none"
  )
  saveRDS(
    object = bpc_after,
    file = paste0(opt$output, "/objects/bpc_after.rds")
  )
}

pdf(paste0(opt$output, "/graphs/rtime/before_after_alignment.pdf"))
par(mfrow = c(2, 1))
# Before retention time alignment
plot(
  bpcs,
  col = group_colors[MsExperiment::sampleData(xchr6)$group],
  main = "Before retention time alignment"
)

# After retention time alignment
plot(
  bpc_after,
  col = group_colors[MsExperiment::sampleData(xchr6)$group],
  main = "After retention time alignment"
)
invisible(dev.off())

# Checking for retention drift in IS
if (check_saved("is_drift_check_before.rds")) {
  is_drift_check_before <- readRDS(
    file = paste0(opt$output, "/objects/is_drift_check_before.rds")
  )
} else {
  is_drift_check_before <- xchr %>%
    Spectra::filterRt(ranges$rt_range) %>%
    Spectra::filterMzRange(ranges$mz_range) %>%
    xcms::chromatogram(
      aggregationFun = "max",
      chromPeaks = "none"
    )
  saveRDS(
    object = is_drift_check_before,
    file = paste0(opt$output, "/objects/is_drift_check_before.rds")
  )
}

if (check_saved("is_drift_check_after.rds")) {
  is_drift_check_after <- readRDS(
    file = paste0(opt$output, "/objects/is_drift_check_after.rds")
  )
} else {
  is_drift_check_after <- xchr6 %>%
    Spectra::filterRt(ranges$rt_range) %>%
    Spectra::filterMzRange(ranges$mz_range) %>%
    xcms::chromatogram(
      aggregationFun = "max",
      chromPeaks = "none"
    )
  saveRDS(
    object = is_drift_check_after,
    file = paste0(opt$output, "/objects/is_drift_check_after.rds")
  )
}

# Checking the adjustment in the IS peak
pdf(paste0(opt$output, "/graphs/rtime/is_before_after_alignment.pdf"))
par(mfrow = c(1, 2))
plot(
  is_drift_check_before,
  col = group_colors[MsExperiment::sampleData(xchr6)$group],
  main = "Before:\nRT: 130 - 175 (s)\nM/z range: 226.9 - 228",
  lwd = 3
)

plot(
  is_drift_check_after,
  col = group_colors[MsExperiment::sampleData(xchr6)$group],
  main = "After:\nRT: 130 - 175 (s)\nM/z range: 226.9 - 228",
  lwd = 3
)
invisible(dev.off())

# From Sattely paper
# Retention time correction was performed using the obiwarp method, with a
# step size of m/z 0.5. Peak alignment was performed with bandwidth
# of 3 seconds and minimum fraction (minfrac) of samples
# necessary for a valid group of 0.5.