# ==============================================================================
# Correspondence
# TODO
# Extract a chromatogram for a m/z range containing internal standard
# Test these settings on the extracted slice
# Do this for several
# ==============================================================================
message("Producing simulated bandwidth plots...")
# Check bandwidth
if (check_saved("chr_1.rds")) {
  chr_1 <- readRDS(file = paste0(opt$output, "/objects/chr_1.rds"))
} else {
  chr_1 <- xcms::chromatogram(
    object = xchr6,
    mz = ranges$mz_range,
    rt = ranges$rt_range + c(-16, 16)
  )
  saveRDS(object = chr_1, file = paste0(opt$output, "/objects/chr_1.rds"))
}

# Test these settings on the extracted slice
pdf(
  paste0(opt$output, "/graphs/internal_standard/is_simul_first_grouping.pdf")
)
density_simul_p <- xcms::plotChromPeakDensity(
  object = chr_1,
  param = xcms::PeakDensityParam(
    sampleGroups = MsExperiment::sampleData(xchr6)$group,
    bw = opt$bw_first_grouping
  )
)
invisible(dev.off())

pdf(
  paste0(opt$output, "/graphs/internal_standard/is_simul_second_grouping.pdf")
)
density_simul_p <- xcms::plotChromPeakDensity(
  object = chr_1,
  param = xcms::PeakDensityParam(
    sampleGroups = MsExperiment::sampleData(xchr6)$group,
    bw = opt$bw_second_grouping
  )
)
invisible(dev.off())
# End of checking for bandwidth to use

message("Performing correspondence...")
# Correspondence
if (check_saved("xchr7.rds")) {
  xchr7 <- readRDS(file = paste0(opt$output, "/objects/xchr7.rds"))
} else {
  xchr7 <- xcms::groupChromPeaks(
    object = xchr6,
    param = xcms::PeakDensityParam(
      sampleGroups = MsExperiment::sampleData(xchr6)$group,
      bw = opt$bw_second_grouping, # 0.5
      minFraction = 0.5, # T0.7
      binSize = 0.01,
      maxFeatures = 1000, # 200
      ppm = opt$ppm_global,
      minSamples = 2 # 1
    )
  )
  saveRDS(object = xchr7, file = paste0(opt$output, "/objects/xchr7.rds"))
}

# TODO Check this for several ions as well
# Extract chromatogram including signal for is
if (check_saved("chr_2.rds")) {
  chr_2 <- readRDS(file = paste0(opt$output, "/objects/chr_2.rds"))
} else {
  chr_2 <- xcms::chromatogram(
    object = xchr7,
    mz = ranges$mz_range,
    rt = ranges$rt_range + c(-16, 16),
    aggregationFun = "max"
  )
  saveRDS(object = chr_2, file = paste0(opt$output, "/objects/chr_2.rds"))
}

message("Producing second grouping bandwidth plots...")
# Setting simulate = FALSE to show the actual correspondence results
pdf(
  paste0(
    opt$output,
    "/graphs/internal_standard/is_non_simul_second_grouping.pdf"
  )
)
density_non_simul_p <- xcms::plotChromPeakDensity(
  object = chr_2,
  simulate = FALSE
)
invisible(dev.off())