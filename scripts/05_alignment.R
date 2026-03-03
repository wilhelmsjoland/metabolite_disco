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
param_msg(xchr5_params)

# ==============================================================================
# Evaluating distribution of anchor peaks --------------------------------------
# ==============================================================================
# Conundrum: unreliable peaks at start and end
# Should I filter these?
# Or should I just take the less variable ones? -> pgm_filt
# Use this part to check if anchor peaks cover the whole RT
pgm <- xcms::adjustRtimePeakGroups(
  xchr5,
  xcms::PeakGroupsParam(minFraction = 0.9)
)
anc_peak_dist_before_path <- file.path(
  opt$output,
  "graphs",
  "rtime",
  "anchor_peak_dist_before.pdf"
)
pdf(anc_peak_dist_before_path)
hist(as.vector(pgm), breaks = 30, main = "Anchor peak RT distribution")
invisible(dev.off())
cli::cli_alert_success(
  "Saved anchor peak distribution to: {.path {anc_peak_dist_before_path}}"
)

# Remove rows where the anchor RT variance across samples with variance cutoff


pgm_apply <- apply(
  X = pgm,
  MARGIN = 1,
  FUN = function(x) {
    sd(x, na.rm = TRUE) < opt$peak_anchor_sd
  }
)
pgm_filt <- pgm[pgm_apply, ]

cli::cli_alert_info("Plotting anchor peak distribution after filtering")
anc_peak_dist_after_path <- file.path(
  opt$output,
  "graphs",
  "rtime",
  "anchor_peak_dist_after.pdf"
)
pdf(anc_peak_dist_after_path)
hist(
  as.vector(pgm_filt),
  breaks = 30,
  main = "Anchor peak RT distribution after filtering"
)
invisible(dev.off())
cli::cli_alert_success(
  paste0(
    "Saved anchor peak distribution after filtering to: ",
    "{.path {anc_peak_dist_after_path}}"
  )
)
# ==============================================================================
# Alignment of retention times
# ==============================================================================
cli::cli_h3("Aligning retention times")

xchr6_path <- file.path(opt$output, "objects", "xchr6.rds")
if (file.exists(xchr6_path)) {
  xchr6 <- readRDS(file = xchr6_path)
  cli::cli_alert_success(
    paste0(
      "Imported grouped peaks object from: ",
      "{.path {xchr6_path}}"
    )
  )
} else {
  cli::cli_alert_info("Grouping peaks")
  xchr6 <- xcms::adjustRtime(
    object = xchr5,
    param = xcms::PeakGroupsParam(
      minFraction = opt$min_fraction_align, # 0.8
      extraPeaks = opt$extra_peaks, # 0 
      smooth = "loess",
      peakGroupsMatrix = pgm_filt,
      span = opt$span, # 0.6 # 0.8
      family = "gaussian",
      # peakGroupsMatrix = matrix(nrow = 0, ncol = 0),
      subset = integer(),
      subsetAdjust = c("average", "previous")
    )
  )
  saveRDS(object = xchr6, file = paste0(opt$output, "/objects/xchr6.rds"))
  cli::cli_alert_success(
    paste0(
      "Saved grouped peaks object to: ",
      "{.path {xchr6_path}}"
    )
  )
}

# ==============================================================================
# Generating alignment base peak chromatograms ---------------------------------
# ==============================================================================
bpc_after_path <- file.path(opt$output, "objects", "bpc_after.rds")
# Extract base peak chromatograms
if (file.exists(bpc_after_path)) {
  bpc_after <- readRDS(file = bpc_after_path)
  cli::cli_alert_success(
    paste0(
      "Imported base peak chromatogram from after ",
      "retention time adjustment from: ",
      "{.path {bpc_after_path}}"
    )
  )
} else {
  cli::cli_alert_info(
    "Generating retention time adjusted base peak chromatogram"
  )
  bpc_after <- xcms::chromatogram(
    xchr6,
    aggregationFun = "max",
    chromPeaks = "none"
  )
  saveRDS(
    object = bpc_after,
    file = bpc_after_path
  )
  cli::cli_alert_success(
    paste0(
      "Saved retention time adjusted base peak chromatogram to: ",
      "{.path {bpc_after_path}}"
    )
  )
}

# ==============================================================================
# Plotting retention time drift ------------------------------------------------
# ==============================================================================
before_after_alignment_path <- file.path(
  opt$output,
  "graphs",
  "rtime",
  "before_after_alignment.pdf"
)
pdf(before_after_alignment_path)
par(mfrow = c(2, 1))
plot(
  bpcs,
  col = group_colors[MsExperiment::sampleData(xchr6)$group],
  main = "Before retention time alignment"
)
plot(
  bpc_after,
  col = group_colors[MsExperiment::sampleData(xchr6)$group],
  main = "After retention time alignment"
)
invisible(dev.off())
cli::cli_alert_success(
  paste0(
    "Saved alignment for base peak chromatogram to ",
    "{.path {before_after_alignment_path}}"
  )
)

# ==============================================================================
# Generating alignment base peak chromatograms for internal standards ----------
# ==============================================================================
is_drift_check_before_path <- file.path(
  opt$output,
  "objects",
  "is_drift_check_before.rds"
)

if (file.exists(is_drift_check_before_path)) {
  is_drift_check_before <- readRDS(file = is_drift_check_before_path)
  cli::cli_alert_success(
    paste0(
      "Imported base peak chromatogram from before ",
      "retention time adjustment for internal standards from: ",
      "{.path {is_drift_check_before_path}}"
    )
  )
} else {
  cli::cli_alert_info(
    paste0(
      "Generating non-retention time adjusted base peak chromatogram ",
      "for internal standards"
    )
  )
  is_drift_check_before <- xchr %>%
    Spectra::filterRt(ranges$rt_range) %>%
    Spectra::filterMzRange(ranges$mz_range) %>%
    xcms::chromatogram(
      aggregationFun = "max",
      chromPeaks = "none"
    )
  saveRDS(
    object = is_drift_check_before,
    file = is_drift_check_before_path
  )
  cli::cli_alert_success(
    paste0(
      "Saved non-retention time adjusted base peak chromatogram ",
      "for internal standards to: ",
      "{.path {is_drift_check_before_path}}"
    )
  )
}

is_drift_check_after_path <- file.path(
  opt$output,
  "objects",
  "is_drift_check_after.rds"
)

if (file.exists(is_drift_check_after_path)) {
  is_drift_check_after <- readRDS(file = is_drift_check_after_path)
  cli::cli_alert_success(
    paste0(
      "Imported base peak chromatogram from after ",
      "retention time adjustment for internal standards from: ",
      "{.path {is_drift_check_after_path}}"
    )
  )
} else {
  cli::cli_alert_info(
    paste0(
      "Generating retention time adjusted base peak chromatogram ",
      "for internal standards"
    )
  )
  is_drift_check_after <- xchr6 %>%
    Spectra::filterRt(ranges$rt_range) %>%
    Spectra::filterMzRange(ranges$mz_range) %>%
    xcms::chromatogram(
      aggregationFun = "max",
      chromPeaks = "none"
    )
  saveRDS(
    object = is_drift_check_after,
    file = is_drift_check_after_path
  )
  cli::cli_alert_success(
    paste0(
      "Saved retention time adjusted base peak chromatogram ",
      "for internal standards to: ",
      "{.path {is_drift_check_after_path}}"
    )
  )
}

# ==============================================================================
# Plotting retention time drift in internal standards --------------------------
# ==============================================================================
is_before_after_alignment_path <- file.path(
  opt$output,
  "graphs",
  "rtime",
  "is_before_after_alignment.pdf"
)
pdf(is_before_after_alignment_path)
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
cli::cli_alert_success(
  paste0(
    "Saved alignment for base peak chromatogram ",
    "of internal standards to: ",
    "{.path {is_before_after_alignment_path}}"
  )
)

# ==============================================================================
# Plotting anchor peaks before and after adjustment ----------------------------
# ==============================================================================
# Before
anchor_peaks_before_path <- file.path(
  opt$output,
  "graphs",
  "rtime",
  "anchor_peaks_before.pdf"
)
pdf(anchor_peaks_before_path)
plotAdjustedRtime(xchr6, adjusted = FALSE)
invisible(dev.off())

# After
anchor_peaks_after_path <- file.path(
  opt$output,
  "graphs",
  "rtime",
  "anchor_peaks_after.pdf"
)
pdf(anchor_peaks_after_path)
plotAdjustedRtime(xchr6, adjusted = TRUE)
invisible(dev.off())