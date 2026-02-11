# ==============================================================================
# Source dependendencies and load libraries -----------------------------------
# ==============================================================================
source("scripts/functions.R")
source("scripts/chem_functions.R")
source("scripts/met_disco_args.R")
suppressWarnings(
  suppressPackageStartupMessages({
    library(tidyverse)
    library(MSnbase)
    library(xcms)
    library(MsExperiment)
    library(RforMassSpectrometry)
    library(Spectra)
    library(Chromatograms)
    library(RColorBrewer)
    library(pheatmap)
    library(QFeatures)
    library(Rdisop)
    library(limma)
    library(BiocParallel)
    library(ggrepel)
    library(ggfortify)
    library(ComplexUpset)
    library(MetaboAnnotation)
    library(CompoundDb)
    library(MetaboCoreUtils)
    library(curl)
    library(gt)
    library(patchwork)
    library(AnnotationHub)
    library(optparse)

    # for examine_biotransformations_hits.R
    library(ComplexHeatmap)
    library(circlize)
    library(rcdk)
    library(rJava)
    library(fingerprint)

  })
)

# ==============================================================================
# Create output folders --------------------------------------------------------
# ==============================================================================
folders <- c(
  "bpc",
  "internal_standard",
  "volcano",
  "upset",
  "feature_boxplot",
  "rtime",
  "filled_peaks",
  "pca",
  "feature_chromatogram_intensity",
  "per_sample_peaks",
  "features",
  "feature_pairs",
  "glycoside",
  "glycoside_feature_pairs"
)
dir.create(res_folder, FALSE, TRUE)
dir.create(file.path(res_folder, "objects"), FALSE, TRUE)
dir.create(file.path(res_folder, "tables"), FALSE, TRUE)
for (folder in folders) {
  dir.create(
    file.path(
      res_folder, "graphs", folder
    ),
    showWarnings = FALSE,
    recursive = TRUE
  )
}
dir.create(file.path("annotation_databases"), FALSE, TRUE)

# ==============================================================================
# Import metadata --------------------------------------------------------------
# ==============================================================================
message("Importing metadata...")
meta <- import_mzml(data_path, meta_file)
if (check_saved("ms_exp.rds")) {
  ms_exp <- readRDS(file = file.path(res_folder, "objects/ms_exp.rds"))
} else {
  ms_exp <- MsExperiment::readMsExperiment(
    spectraFiles = meta$path,
    sampleData = meta
  )
  saveRDS(object = ms_exp, file = file.path(res_folder, "objects/ms_exp.rds"))
}

# ==============================================================================
# Set colors for groups
# ==============================================================================
message("Setting colors for groups...")
groups_to_use <- unique(MsExperiment::sampleData(ms_exp)$group)
group_colors <- paste0(
  RColorBrewer::brewer.pal(
    n = length(groups_to_use),
    "Set1")[seq_along(groups_to_use)]
)
group_colors <- setNames(group_colors, groups_to_use)

# ==============================================================================
# Create and plot base peak chromatograms
# ==============================================================================
message("Creating base peak chromatograms...")
if (check_saved("bpcs.rds")) {
  bpcs <- readRDS(file = file.path(res_folder, "objects/bpcs.rds"))
} else {
  bpcs <- xcms::chromatogram(ms_exp, aggregationFun = "max")
  saveRDS(object = bpcs, file = file.path(res_folder, "objects/bpcs.rds"))
}

message("Plotting base peak chromatograms...")
pdf(file.path(res_folder, "graphs/bpc/raw_bpc.pdf"))
par(mar = c(4, 4, 3, 2))
plot(
  x = bpcs,
  col = group_colors[MsExperiment::sampleData(ms_exp)$group],
  main = "Base peak chromatogram"
)
invisible(dev.off())

# ==============================================================================
# Create a heatmap of base peak intensities
# ==============================================================================
# Calculate correlation on the log2 transformed base peak intensities
bpcs_bin <- bin(bpcs, binSize = 1)
cormat <- cor(
  log2(
    do.call(
      cbind,
      lapply(bpcs_bin, intensity)
    )
  ),
  use = "complete.obs" # Because NAs
)

cormat_rownames <- basename(Biobase::pData(xcms::phenoData(bpcs))$path)
colnames(cormat) <- rownames(cormat) <- cormat_rownames

# Define which phenodata columns should be highlighted in the plot
ann <- data.frame(group = bpcs_bin$group)
rownames(ann) <- cormat_rownames

# Perform the cluster analysis

cormat_p <- pheatmap::pheatmap(
  cormat,
  annotation = ann,
  annotation_color = list(group = group_colors),
  silent = TRUE
)

ggplot2::ggsave(
  filename = paste0(res_folder, "/graphs/bpc/raw_bpc_hmp.pdf"),
  plot = cormat_p,
  device = "pdf",
  height = 10,
  width = 10,
  units = "in"
)

# ==============================================================================
# Inspect internal standard prior to peak-calling
# Define the rt and m/z range of the peak area
# ==============================================================================
message("Inspecting internal standard peaks prior to peak-calling...")
mz_theory <- get_theory_mz(chem_form = internal_standard, adduct = adduct)
mz_range <- get_short_mz_range(mz_theory, mz_window = 0.02)
if (check_saved("is_chr.rds")) {
  is_chr <- readRDS(file = paste0(res_folder, "/objects/is_chr.rds"))
} else {
  is_chr <- xcms::chromatogram(
    object = ms_exp,
    mz = mz_range,
    aggregationFun = "sum"
  )
  saveRDS(object = is_chr, file = paste0(res_folder, "/objects/is_chr.rds"))
}
ranges <- get_rt_mz_range(chromatogram = is_chr, rt_window = 0.02)

# Wide IS chromatogram
pdf(paste0(res_folder, "/graphs/internal_standard/all_is_wide.pdf"))
plot(x = is_chr, col = group_colors[is_chr$group], lwd = 3)
legend("topright", legend = names(group_colors), col = group_colors, pch = 16)
invisible(dev.off())

# Get the IS XIC
if (check_saved("is_eic.rds")) {
  is_eic <- readRDS(file = paste0(res_folder, "/objects/is_eic.rds"))
} else {
  is_eic <- xcms::chromatogram(
    object = ms_exp,
    mz = ranges$mz_range,
    rt = ranges$rt_range,
    aggregationFun = "sum"
  )
  saveRDS(object = is_eic, file = paste0(res_folder, "/objects/is_eic.rds"))
}

# All IS EICs together
pdf(paste0(res_folder, "/graphs/internal_standard/all_is.pdf"))
plot(x = is_eic, col = group_colors[is_eic$group], lwd = 3)
legend("topleft", legend = names(group_colors), col = group_colors, pch = 16)
invisible(dev.off())

# Individual IS XICs
for (i in seq_along(is_eic)) {
  plot(
    x = is_eic[, i],
    col = group_colors[
      names(group_colors) %in% dplyr::filter(
        meta,
        rownames(meta) == colnames(bpcs)[i]
      )$group
    ],
    lwd = 3,
    main = colnames(bpcs)[i]
  )
  # stop_loop <- readline("Enter for next, 'break' for stop: ")
  stop_loop <- "break"
  if (stop_loop == "break") {
    break
  }
}
invisible(dev.off())

# ==============================================================================
# - Ensure IS peak is chosen
# - Determine max and min peak width for CentWaveParam() from the IS.
# TODO
# do this for several peaks and not only the IS
# ==============================================================================
message(
  "Determining minimal and maximal peakwidth ",
  "based on the internal standard..."
)

if (check_saved("is_eic_wide.rds")) {
  is_eic_wide <- readRDS(file = paste0(res_folder, "/objects/is_eic_wide.rds"))
} else {
  is_eic_wide <- xcms::chromatogram(
    ms_exp,
    mz = mz_range + c(-0.05, 0.05),
    rt = ranges$rt_range + c(-16, 16),
    aggregationFun = "sum"
  )
  saveRDS(
    object = is_eic_wide,
    file = paste0(res_folder, "/objects/is_eic_wide.rds")
  )
}

# Run peak detection on the EIC
if (check_saved("is_chr2.rds")) {
  is_chr2 <- readRDS(file = paste0(res_folder, "/objects/is_chr2.rds"))
} else {
  is_chr2 <- xcms::findChromPeaks(
    object = is_eic_wide,
    param = xcms::CentWaveParam(
      ppm = ppm_global,
      peakwidth = c(2, 20),
      prefilter = c(1, 1),
      snthresh = sn_threshold, # 10
      mzCenterFun = "wMean",
      mzdiff = 0.001,
      integrate = 2,
      noise = 1000,
      verboseBetaColumns = TRUE
    ),
    ms_level = 1
  )
  saveRDS(object = is_chr2, file = paste0(res_folder, "/objects/is_chr2.rds"))
}

# Calculate peakwidth
is_peaks <- tibble::as_tibble(
  x = xcms::chromPeaks(is_chr2),
  rownames = "rownames"
) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(delta_rt = rtmax - rtmin) %>%
  dplyr::ungroup()


# Min: to half of some peaks in the datasets
# Max: 2-4x times the average size
is_min_peak_width <- min(is_peaks$delta_rt, na.rm = TRUE)
is_max_peak_width <- max(is_peaks$delta_rt, na.rm = TRUE)

min_peak_width <- quantile(is_peaks$delta_rt, 0.05, na.rm = TRUE) * 0.3 # 0.3
max_peak_width <- quantile(is_peaks$delta_rt, 0.95, na.rm = TRUE) * 4 # 4

message(
  "===========================================================================",
  "\nInternal standard: ", "C7H8O2", ", Theoretical m/z: ", round(mz_theory, 3),
  "\n\tMinimal peak width of IS: ", round(is_min_peak_width, 3),
  "\n\tMaximal peak width of IS: ", round(is_max_peak_width, 3),
  "\nPeak widths used for peak picking: ",
  "\n\t Minimal width: ", round(min_peak_width, 3),
  "\n\t Maximal width: ", round(max_peak_width, 3), "\n",
  "===========================================================================",
  sep = ""
)

# ==============================================================================
# Call peaks on whole dataset with parameters
# ==============================================================================
message("Calling peaks...")

if (check_saved("xchr.rds")) {
  xchr <- readRDS(file = paste0(res_folder, "/objects/xchr.rds"))
} else {
  xchr <- xcms::findChromPeaks(
    object = ms_exp,
    BPPARAM = BiocParallel::bpparam(),
    return.type = "XCMSnExp",
    ms_level = 1L,
    param = xcms::CentWaveParam(
      ppm = ppm_global,
      peakwidth = c(min_peak_width, max_peak_width),
      snthresh = sn_threshold, # 10
      prefilter = c(4, 1000), # k pks (left) over intens (right) # c(3, 100)
      mzCenterFun = "wMean",
      integrate = 2,
      mzdiff = 0.001,
      fitgauss = FALSE,
      noise = 1000,
      verboseColumns = TRUE,
      roiList = list(),
      firstBaselineCheck = TRUE,
      roiScales = numeric(),
      extendLengthMSW = FALSE,
      verboseBetaColumns = TRUE
    )
  )
  saveRDS(object = xchr, file = paste0(res_folder, "/objects/xchr.rds"))
}

plot(
  bpcs,
  col = group_colors[MsExperiment::sampleData(xchr)$group],
  main = "Base peak chromatogram after peak picking"
)
legend(
  "topright",
  col = unique(group_colors[MsExperiment::sampleData(xchr)$group]),
  legend = unique(names(group_colors[MsExperiment::sampleData(ms_exp)$group])),
  pch = 16
)
invisible(dev.off())

# Use this for the beta-distribution parameters
xchr_data <- tibble::as_tibble(xcms::chromPeaks(xchr), rownames = "peak")

# ==============================================================================
# Inspect peaks
# TODO
# Extract some peaks here and check quality of peak picking
# ==============================================================================
message("Inspecting called peaks...")

# peaks_to_inspect <- as_tibble(chromPeaks(xchr[1:2]), rownames = "peak") %>%
#     tidyr::drop_na(beta_snr) %>%
#     dplyr::arrange(desc(into)) %>%
#     dplyr::slice(1:10) %>%
#     dplyr::pull(peak) %>%
#     stringr::str_remove_all(., "[A-Za-z]") %>%
#     as.numeric(.)
#
# for (peak in peaks_to_inspect) {
#     inspect_peaks <- inspect_peak(
#         chromatogram = xchr[1:2],
#         peak_data = as_tibble(chromPeaks(xchr[1:2]), rownames = "peak"),
#         peak_idx = peak,
#         sample_no = FALSE,
#         save_graph = TRUE
#     )
# }


# Poor peaks: beta_cor < 0.5 (or even < 0.2
# Good peaks: beta_snr > 7
# Keep signal to noise at 10, and filter by beta_cor < 0.5, and beta_snr > 7?

# These are the per-sample peak counts
inspect_peak_intensity(chr_data = xchr, value = into, save_graph = TRUE)

# ==============================================================================
# Refine peaks
# TODO
# Let these be part of the pipeline as well and as a choice to do or not
# ==============================================================================
if (check_saved("xchr2.rds")) {
  xchr2 <- readRDS(paste0(res_folder, "/objects/xchr2.rds"))
} else {
  xchr2 <- xcms::refineChromPeaks(
    object = xchr,
    param = xcms::CleanPeaksParam(maxPeakwidth = max_peak_width)
  )
  saveRDS(object = xchr2, file = paste0(res_folder, "/objects/xchr2.rds"))
}

if (check_saved("xchr3.rds")) {
  xchr3 <- readRDS(paste0(res_folder, "/objects/xchr3.rds"))
} else {
  xchr3 <- xcms::refineChromPeaks(
    object = xchr2,
    param = xcms::FilterIntensityParam(
      threshold = 0,
      nValues = 1L,
      value = "maxo"
    )
  )
  saveRDS(object = xchr3, file = paste0(res_folder, "/objects/xchr3.rds"))
}

if (check_saved("xchr4.rds")) {
  xchr4 <- readRDS(paste0(res_folder, "/objects/xchr4.rds"))
} else {
  xchr4 <- xcms::refineChromPeaks(
    object = xchr, # xchr3
    param = xcms::MergeNeighboringPeaksParam(
      expandRt = 0.25,
      expandMz = 0,
      ppm =  ppm_global,
      minProp = 0.95 # between 0 & 1
    )
  )
  saveRDS(object = xchr4, file = paste0(res_folder, "/objects/xchr4.rds"))
}

# TODO
# Make a check for if these exist
clean_removed_peaks <- dplyr::anti_join(
  x = tibble::as_tibble(xcms::chromPeaks(xchr), rownames = "peaks"),
  y = tibble::as_tibble(xcms::chromPeaks(xchr2), rownames = "peaks"),
  by = "peaks"
)
clean_removed_peaks %>% nrow()

intensity_removed_peaks <- dplyr::anti_join(
  x = tibble::as_tibble(xcms::chromPeaks(xchr2), rownames = "peaks"),
  y = tibble::as_tibble(xcms::chromPeaks(xchr3), rownames = "peaks"),
  by = "peaks"
)
intensity_removed_peaks %>% nrow()

merged_peaks <- dplyr::anti_join(
  x = tibble::as_tibble(xcms::chromPeaks(xchr3), rownames = "peaks"),
  y = tibble::as_tibble(xcms::chromPeaks(xchr4), rownames = "peaks"),
  by = "peaks"
)

# ==============================================================================
# Alignment of retention times
# ==============================================================================
message("Aligning retention times across samples...")

# TODO Change these to more broad so I don't get a million anchor peaks?
if (check_saved("xchr5.rds")) {
  xchr5 <- readRDS(paste0(res_folder, "/objects/xchr5.rds"))
} else {
  xchr5 <- xcms::groupChromPeaks(
    object = xchr, # xchr4
    param = xcms::PeakDensityParam(
      sampleGroups = MsExperiment::sampleData(xchr)$group,
      bw = bw_first_grouping,
      minFraction = 0.5, # 0.7
      binSize = 0.01,
      maxFeatures = 1000, # 200
      ppm = ppm_global,
      minSamples = 2 # 1
    )
  )
  saveRDS(object = xchr5, file = paste0(res_folder, "/objects/xchr5.rds"))
}

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

# Alignment
if (check_saved("xchr6.rds")) {
  xchr6 <- readRDS(paste0(res_folder, "/objects/xchr6.rds"))
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
  saveRDS(object = xchr6, file = paste0(res_folder, "/objects/xchr6.rds"))
}

# Checking for retention drift
# Extract base peak chromatograms
if (check_saved("bpc_after.rds")) {
  bpc_after <- readRDS(file = paste0(res_folder, "/objects/bpc_after.rds"))
} else {
  bpc_after <- xcms::chromatogram(
    xchr6,
    aggregationFun = "max",
    chromPeaks = "none"
  )
  saveRDS(
    object = bpc_after,
    file = paste0(res_folder, "/objects/bpc_after.rds")
  )
}

pdf(paste0(res_folder, "/graphs/rtime/before_after_alignment.pdf"))
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
    file = paste0(res_folder, "/objects/is_drift_check_before.rds")
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
    file = paste0(res_folder, "/objects/is_drift_check_before.rds")
  )
}

if (check_saved("is_drift_check_after.rds")) {
  is_drift_check_after <- readRDS(
    file = paste0(res_folder, "/objects/is_drift_check_after.rds")
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
    file = paste0(res_folder, "/objects/is_drift_check_after.rds")
  )
}

# Checking the adjustment in the IS peak
pdf(paste0(res_folder, "/graphs/rtime/is_before_after_alignment.pdf"))
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
  chr_1 <- readRDS(file = paste0(res_folder, "/objects/chr_1.rds"))
} else {
  chr_1 <- xcms::chromatogram(
    object = xchr6,
    mz = ranges$mz_range,
    rt = ranges$rt_range + c(-16, 16)
  )
  saveRDS(object = chr_1, file = paste0(res_folder, "/objects/chr_1.rds"))
}

# Test these settings on the extracted slice
pdf(
  paste0(res_folder, "/graphs/internal_standard/is_simul_first_grouping.pdf")
)
density_simul_p <- xcms::plotChromPeakDensity(
  object = chr_1,
  param = xcms::PeakDensityParam(
    sampleGroups = MsExperiment::sampleData(xchr6)$group,
    bw = bw_first_grouping
  )
)
invisible(dev.off())

pdf(
  paste0(res_folder, "/graphs/internal_standard/is_simul_second_grouping.pdf")
)
density_simul_p <- xcms::plotChromPeakDensity(
  object = chr_1,
  param = xcms::PeakDensityParam(
    sampleGroups = MsExperiment::sampleData(xchr6)$group,
    bw = bw_second_grouping
  )
)
invisible(dev.off())
# End of checking for bandwidth to use

message("Performing correspondence...")
# Correspondence
if (check_saved("xchr7.rds")) {
  xchr7 <- readRDS(file = paste0(res_folder, "/objects/xchr7.rds"))
} else {
  xchr7 <- xcms::groupChromPeaks(
    object = xchr6,
    param = xcms::PeakDensityParam(
      sampleGroups = MsExperiment::sampleData(xchr6)$group,
      bw = bw_second_grouping, # 0.5
      minFraction = 0.5, # T0.7
      binSize = 0.01,
      maxFeatures = 1000, # 200
      ppm = ppm_global,
      minSamples = 2 # 1
    )
  )
  saveRDS(object = xchr7, file = paste0(res_folder, "/objects/xchr7.rds"))
}

# TODO Check this for several ions as well
# Extract chromatogram including signal for is
if (check_saved("chr_2.rds")) {
  chr_2 <- readRDS(file = paste0(res_folder, "/objects/chr_2.rds"))
} else {
  chr_2 <- xcms::chromatogram(
    object = xchr7,
    mz = ranges$mz_range,
    rt = ranges$rt_range + c(-16, 16),
    aggregationFun = "max"
  )
  saveRDS(object = chr_2, file = paste0(res_folder, "/objects/chr_2.rds"))
}

message("Producing second grouping bandwidth plots...")
# Setting simulate = FALSE to show the actual correspondence results
pdf(
  paste0(
    res_folder,
    "/graphs/internal_standard/is_non_simul_second_grouping.pdf"
  )
)
density_non_simul_p <- xcms::plotChromPeakDensity(
  object = chr_2,
  simulate = FALSE
)
invisible(dev.off())

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
  dplyr::mutate(feat_extract = rbind(
    c(
      mzmed - 0.0015,
      mzmed + 0.0015,
      rtmin - 2,
      rtmax + 2
    ))
  ) %>%
  dplyr::pull(feat_extract)

message(
  "Features with NAs prior to gap filling: ",
  sum(is.na(featureValues(xchr7)))
)

# Perform gap filling
if (check_saved("xchr8.rds")) {
  xchr8 <- readRDS(file = paste0(res_folder, "/objects/xchr8.rds"))
} else {
  xchr8 <- xcms::fillChromPeaks(
    object = xchr7,
    param = xcms::ChromPeakAreaParam()
  )
  saveRDS(object = xchr8, file = paste0(res_folder, "/objects/xchr8.rds"))
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
  chrs_na <- readRDS(file = paste0(res_folder, "/objects/chrs_na.rds"))
} else {
  chrs_na <- xcms::chromatogram(
    xchr8,
    mz = chrs_na_feat[, c("mzmin", "mzmax")],
    rt = chrs_na_feat[, c("rtmin", "rtmax")] # probably increase this a little
  )
  saveRDS(object = chrs_na, file = paste0(res_folder, "/objects/chrs_na.rds"))
}

# ==============================================================================
# Filtering features and input to SummarizedExperiment -------------------------
# ==============================================================================
message("Filtering features based on missingness...")
group_factor <- MsExperiment::sampleData(xchr8)$group
group_factor <- as.factor(group_factor)
if (check_saved("xchr9.rds")) {
  xchr9 <- readRDS(file = paste0(res_folder, "/objects/xchr9.rds"))
} else {
  xchr9 <- QFeatures::filterFeatures(
    xchr8,
    xcms::PercentMissingFilter(
      threshold = missing_threshold,
      f = group_factor
    )
  )
  saveRDS(object = xchr9, file = paste0(res_folder, "/objects/xchr9.rds"))
}

# ==============================================================================
# Median scaling & PCA ------------------------------------------------------
# ==============================================================================
message("Producing PCA plots...")

res <- xcms::quantify(
  xchr9,
  method = "sum",
  value = "into",
  filled = FALSE,
  missing = "rowmin_half" # 0 ?
)

SummarizedExperiment::assays(res)$raw_filled <- xcms::featureValues(
  xchr9,
  method = "sum",
  value = "into",
  filled = TRUE,
  missing = "rowmin_half" # 0 ?
)

# Compute median and generate normalization factor
mdns <- apply(
  SummarizedExperiment::assay(res, "raw"),
  MARGIN = 2,
  median,
  na.rm = TRUE
)
nf_mdn <- mdns / median(mdns)

# Dividing dataset by median of median and creating a new assay
SummarizedExperiment::assays(res)$norm <- sweep(
  SummarizedExperiment::assay(res, "raw"),
  MARGIN = 2,
  nf_mdn,
  "/"
)

# Compute median and generate normalization factor
mdns <- apply(
  SummarizedExperiment::assay(res, "raw_filled"),
  MARGIN = 2,
  median,
  na.rm = TRUE
)
nf_mdn <- mdns / median(mdns)

# Dividing dataset by median of median and creating a new assay
SummarizedExperiment::assays(res)$norm_filled <- sweep(
  SummarizedExperiment::assay(res, "raw_filled"),
  MARGIN = 2,
  nf_mdn,
  "/"
)

# Log2 transform and scale data
vals <- SummarizedExperiment::assay(res, "raw_filled") %>%
  log2() %>%
  t() %>%
  scale(center = TRUE, scale = TRUE) %>%
  as.matrix(.)

pca_res <- prcomp(vals, scale = FALSE, center = FALSE)

# Data before normalization
vals_st <- cbind(vals, phenotype = res$group)
pca_raw <- autoplot(
  pca_res,
  data = vals_st,
  colour = "phenotype",
  scale = 0,
  size = 3) +
  ggplot2::scale_color_manual(values = group_colors) +
  ggplot2::theme_bw() +
  ggplot2::labs(title = "Before normalization")

# Data after normalization
vals_norm <- SummarizedExperiment::assay(res, "norm_filled") %>%
  log2() %>%
  t() %>%
  scale(center = TRUE, scale = TRUE) %>%
  as.matrix(.)

pca_res_norm <- prcomp(vals_norm, scale = FALSE, center = FALSE)
vals_st_norm <- cbind(vals_norm, phenotype = res$group)
pca_adj <- autoplot(
  pca_res_norm,
  data = vals_st_norm,
  colour = "phenotype",
  scale = 0,
  size = 3) +
  ggplot2::scale_color_manual(values = group_colors) +
  ggplot2::theme_bw() +
  ggplot2::labs(title = "After normalization")

norm_filled_12_pca_p <- pca_raw / pca_adj +
  patchwork::plot_layout(guides = "collect")

ggplot2::ggsave(
  filename = paste0(res_folder, "/graphs/pca/norm_filled_pca_1_2.pdf"),
  plot = norm_filled_12_pca_p,
  device = "pdf",
  height = 6,
  width = 6,
  units = "in"
)

# PC 3-4 before & after normalization
pca_raw <- autoplot(
  pca_res,
  data = vals_st,
  colour = "phenotype",
  x = 3,
  y = 4,
  scale = 0,
  size = 3) +
  ggplot2::scale_color_manual(values = group_colors) +
  ggplot2::theme_bw() +
  ggplot2::labs(title = "Before normalization")

pca_adj <- autoplot(
  pca_res_norm,
  data = vals_st_norm,
  colour = "phenotype",
  x = 3,
  y = 4,
  scale = 0,
  size = 3) +
  ggplot2::scale_color_manual(values = group_colors) +
  ggplot2::theme_bw() +
  ggplot2::labs(title = "After normalization")

norm_filled_34_pca_p  <- pca_raw / pca_adj +
  patchwork::plot_layout(guides = "collect")

ggplot2::ggsave(
  filename = paste0(res_folder, "/graphs/pca/norm_filled_pca_3_4.pdf"),
  plot = norm_filled_34_pca_p,
  device = "pdf",
  height = 6,
  width = 6,
  units = "in"
)

# ==============================================================================
# Limnear models using limma
# ==============================================================================
message("Running linear models...")

group_used <- factor(meta$group)
design <- model.matrix(~ 0 + group_used)
colnames(design) <- levels(group_used)

comparisons <- combn(
  x = levels(group_used),
  m = 2,
  simplify = TRUE) %>%
  t(.) %>%
  tibble::as_tibble(.) %>%
  dplyr::mutate(comp = paste0(V1, "-", V2)) %>%
  dplyr::pull(comp)

contrasts_mat <- limma::makeContrasts(
  contrasts = comparisons,
  levels = design
)

fit <- limma::lmFit(
  log2(SummarizedExperiment::assay(res, "norm")), # don't use imputed
  design = design
)
fit <- limma::contrasts.fit(fit, contrasts_mat)
fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)

limma_res <- list()
for (i in comparisons) {
  tmp <- limma::topTable(
    fit = fit,
    coef = i,
    number = Inf,
    adjust.method = "BH",
    sort.by = "none") %>%
    tibble::as_tibble(., rownames = "feature") %>%
    dplyr::mutate(contrast = i) %>%
    dplyr::left_join(
      x = .,
      y = SummarizedExperiment::rowData(res) %>%
        tibble::as_tibble(., rownames = "feature"),
      by = "feature"
    )
  limma_res[[i]] <- tmp
}

full_limma <- tibble::tibble()
for (i in names(limma_res)) {
  tmp_tib <- limma_res[[i]]
  full_limma <- dplyr::bind_rows(full_limma, tmp_tib)
}

message("Producing volcano plots...")
volc_plot_list <- list()
for (i in names(limma_res)) {

  tmp <- limma_res[[i]]

  tmp_tib <- tmp %>%
    dplyr::mutate(
      label.p = dplyr::if_else(
        adj.P.Val < p_value_global & abs(logFC) > quantile(abs(logFC), 0.99),
        mzmed,
        NA
      ),
      direction = dplyr::case_when(
        logFC >= 0 & adj.P.Val < p_value_global ~ "Up",
        logFC < 0 & adj.P.Val < p_value_global ~ "Down",
        adj.P.Val >= p_value_global ~ "ns",
        TRUE ~ as.character("check")
      )
    ) %>%
    dplyr::mutate(direction = forcats::fct_relevel(
      direction,
      c(
        "Up",
        "ns",
        "Down"
      )
    )) %>%
    tidyr::drop_na(logFC)

  tmp_p <- tmp_tib %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = logFC,
        y = -log10(adj.P.Val),
        color = direction
      )
    ) +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = c(
      "Up" = "firebrick",
      "Down" = "cornflowerblue",
      "ns" = "grey"
    )) +
    ggrepel::geom_label_repel(
      data = tidyr::drop_na(tmp_tib, label.p) %>%
        dplyr::arrange(dplyr::desc(abs(logFC))) %>%
        dplyr::slice(1:50),
      ggplot2::aes(label = round(label.p, 2)),
      size = 3,
      max.overlaps = 100,
      box.padding = 0.5,
      color = "black"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    ) +
    ggplot2::labs(
      title = i,
      subtitle = paste0(
        "Rounded mass-to-charge ratios with adj.p < ",
        p_value_global,
        " are labelled"
      ),
      x = "Log2 fold change",
      y = "-Log10 adjusted p-value"
    )

  volc_plot_list[[i]] <- tmp_p

  ggplot2::ggsave(
    filename = paste0(res_folder, "/graphs/volcano/", i, ".pdf"),
    plot = tmp_p,
    device = "pdf",
    height = 10,
    width = 10,
    units = "in"
  )
}

message("Saving intensity information to tables...")
# Creating tables of all output data and saving to tables
# TODO Fix this dumb logic here or use as witch statement?
assay_names <- names(SummarizedExperiment::assays(res))
for (i in assay_names) {
  if (i %in% c("norm", "norm_filled")) {
    full_data <- dplyr::left_join(
      x = SummarizedExperiment::rowData(res) %>%
        tibble::as_tibble(., rownames = "feature"),
      y = SummarizedExperiment::assay(res, i) %>%
        log2() %>%
        t() %>%
        scale(., center = TRUE, scale = TRUE) %>%
        t() %>%
        tibble::as_tibble(., rownames = "feature"),
      by = "feature"
    )
  } else if (i %in% c("raw", "raw_filled")) {
    full_data <- dplyr::left_join(
      x = SummarizedExperiment::rowData(res) %>%
        tibble::as_tibble(., rownames = "feature"),
      y = SummarizedExperiment::assay(res, i) %>%
        tibble::as_tibble(., rownames = "feature"),
      by = "feature"
    )
  }

  assign(
    x = paste0("full_", i),
    value = full_data,
    envir = .GlobalEnv
  )

  readr::write_csv(
    x = full_data,
    file = paste0(res_folder, "/tables/full_", i, ".csv"),
    na = "NA",
    col_names = TRUE,
    append = FALSE
  )
}

# ==============================================================================
# Produce upset plots ----------------------------------------------------------
# ==============================================================================
message("Producing upset plots...")
upset_tib <- full_limma %>%
  dplyr::select(feature, contrast, adj.P.Val) %>%
  tidyr::pivot_wider(
    names_from = "contrast",
    values_from = "adj.P.Val"
  ) %>%
  dplyr::mutate(
    dplyr::across(
      .cols = 2:ncol(.),
      .fns = ~ dplyr::if_else(
        . < p_value_global,
        TRUE,
        FALSE
      )
    )
  )

upset_p <- upset_tib %>%
  ComplexUpset::upset(
    intersect = comparisons,
    name = paste0("Features with p adjusted < ", p_value_global),
    width_ratio = 0.15,
    base_annotations = list(
      "Intersecting features" = ComplexUpset::intersection_size()
    )
  )

# upset_tib %>%
#     dplyr::select(c("feature", tidyselect::all_of(upset_comps))) %>%
#     ComplexUpset::upset(
#         intersect = upset_comps,
#         name = paste0("Features with p adjusted < ", p_value_global),
#         width_ratio = 0.15,
#         base_annotations = list(
#             "Intersecting features" = ComplexUpset::intersection_size()
#         ),
#         min_degree = 3
#     )

ggplot2::ggsave(
  filename = paste0(res_folder, "/graphs/upset/upset_pdf"),
  plot = upset_p,
  device = "pdf",
  height = 7,
  width = 22,
  units = "in"
)

# ==============================================================================
# Find intersecting features ---------------------------------------------------
# ==============================================================================
message("Finding intersecting features...")
# TODO Create all interesting comparisons
upset_comp <- c(
  "bu_wt_apiin-bu_wt_control", # 1
  "bu_mutant_control-bu_wt_apiin", # 1
  "bu_mutant_apiin-bu_wt_apiin" # , # 1
  # "bu_mutant_apiin-bu_wt_control" # 2
)

# Add to list
upset_comps <- list()
intersecting_feats <- find_intersect_feat(
  data = upset_tib,
  set = upset_comp,
  full_set = comparisons
)
upset_comps[[
  stringr::str_flatten(upset_comp, collapse = "*")
]] <- intersecting_feats$feature
# TODO END of part that needs fixing

all.int.comps <- upset_tib %>%
  tidyr::pivot_longer(cols = 2:ncol(.)) %>%
  dplyr::filter(value == TRUE) %>%
  dplyr::group_by(feature) %>%
  dplyr::filter(dplyr::n() >= 2) %>%
  dplyr::pull(feature) %>%
  unique(.)

intersect_data <- full_norm_filled %>%
  tidyr::pivot_longer(cols = tidyselect::all_of(rownames(meta))) %>%
  dplyr::filter(feature %in% intersecting_feats$feature) %>%
  dplyr::left_join(
    x = .,
    y = meta %>%
      tibble::as_tibble(., rownames = "sample"),
    by = c("name" = "sample")
  ) %>%
  dplyr::mutate(group = forcats::fct_relevel(
    group,
    c(
      "bu_wt_control",
      "bu_wt_apiin",
      "bu_mutant_control",
      "bu_mutant_apiin"
    )
  ))

limma_p_res <- full_limma %>%
  dplyr::select(feature, adj.P.Val, contrast) %>%
  dplyr::mutate(
    group1 = stringr::str_split_i(contrast, "-", 1),
    group2 = stringr::str_split_i(contrast, "-", 2),
  ) %>%
  dplyr::select(-contrast) %>%
  dplyr::filter(feature %in% intersecting_feats$feature) %>%
  rstatix::add_significance(p.col = "adj.P.Val") %>%
  find_y_position(
    test_df = .,
    df = intersect_data,
    formula = "value ~ group",
    fun_data = "max",
    grouping = "feature"
  )

message("Producing significant intersecting feature boxplots...")
# Intersecting feats individually
for (i in intersecting_feats$feature) {
  tmp_inter_data <- full_norm_filled %>%
    tidyr::pivot_longer(cols = tidyselect::all_of(rownames(meta))) %>%
    dplyr::filter(feature %in% i) %>%
    dplyr::left_join(
      x = .,
      y = meta %>%
        tibble::as_tibble(., rownames = "sample"),
      by = c("name" = "sample")
    ) %>%
    dplyr::mutate(group = forcats::fct_relevel(
      group,
      c(
        "bu_wt_control",
        "bu_wt_apiin",
        "bu_mutant_control",
        "bu_mutant_apiin"
      ))
    )

  tmp_inter_limma <- full_limma %>%
    dplyr::select(feature, adj.P.Val, contrast) %>%
    dplyr::mutate(
      group1 = stringr::str_split_i(contrast, "-", 1),
      group2 = stringr::str_split_i(contrast, "-", 2),
    ) %>%
    dplyr::select(-contrast) %>%
    dplyr::filter(feature %in% i) %>%
    rstatix::add_significance(p.col = "adj.P.Val") %>%
    find_y_position(
      test_df = .,
      df = tmp_inter_data,
      formula = "value ~ group",
      fun_data = "max"
    )

  tmp_inter_p <- tmp_inter_data %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = group,
        y = value
      )) +
    ggplot2::geom_boxplot(
      ggplot2::aes(fill = group),
      outliers = FALSE
    ) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(
      facets = ~ feature,
      scales = "free_y"
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(c(0.1, 0.15))) +
    ggplot2::guides(x = ggplot2::guide_axis(angle = -45)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      y = "Log2 median-scaled intensity"
    ) +
    ggpubr::geom_bracket(
      data = tmp_inter_limma %>%
        dplyr::filter(!adj.P.Val.signif %in% c("ns")),
      ggplot2::aes(
        xmin = group1,
        xmax = group2,
        label = adj.P.Val.signif,
        y.position = y.pos * 1.005
      ),
      step.increase = 0.1
    )

  ggplot2::ggsave(
    filename = paste0(res_folder, "/graphs/feature_boxplot/", i, ".pdf"),
    plot = tmp_inter_p,
    device = "pdf",
    height = 5,
    width = 5,
    units = "in"
  )
}

# ==============================================================================
# m/z predictions all ----------------------------------------------------------
# ==============================================================================
message("Expanding possible adducts...")
xchr9_defs <- xcms::featureDefinitions(xchr9) %>%
  tibble::as_tibble(., rownames = "feature")

xchr9_mzs <- xchr9_defs$mzmed
names(xchr9_mzs) <- xchr9_defs$feature

possible_adducts <- MetaboCoreUtils::mz2mass(
  xchr9_mzs,
  adduct = adducts(polarity = polarity)
) %>%
  tibble::as_tibble(., rownames = "feature") %>%
  tidyr::pivot_longer(
    cols = 2:ncol(.),
    names_to = "adduct",
    values_to = "mass"
  ) %>%
  dplyr::left_join(
    x = .,
    y = xchr9_defs,
    by = "feature"
  )

if (nrow(xchr9_defs) * 17 == nrow(possible_adducts)) {
  message("Adducts expanded")
} else {
  warning("Adducts did not correctly expanded")
}

message("Importing biotransformation file...")
bio_transf <- import_biotransform_meta(
  file = paste0(data_path, "/", biotransf_file)
)

if (!exists("biotransf_append")) {
  source("R/rpairs_parse.R")
} else {
  "'biotransf_append' already exists."
}

bio_transf2 <- dplyr::bind_rows(
  bio_transf,
  biotransf_append
)

message(
  "Predicting potential biotransformations based on:\n\t",
  "Biotransformation database: ", biotransf_file,
  "\n\tppm: ", ppm_global,
  sep = ""
)

# TODO
# CHECK IF THIS MAKES SENSE NOW!!!!!
# Fix so the observed ppm is added
# Also filter noisy features with the filt_features() function I made
matched_diffs <- pred_biot(
  data = possible_adducts,
  biotransf_data = bio_transf, # bio_transf2
  tolerance_ppm = 1 # ppm_global
)

message("Writing predictions to table...")
writexl::write_xlsx(
  x = matched_diffs,
  path = paste0(res_folder, "/tables/matched_diffs.xlsx"),
  col_names = TRUE,
  format_headers = TRUE,
  use_zip64 = FALSE
)

message("
 based on significant features in contrasts:\n\t",
  paste0(upset_comps, "\n\t"),
  sep = ""
)

filt_match_diffs <- matched_diffs %>%
  dplyr::filter(
    dplyr::if_any(
      tidyselect::all_of(c("feat1", "feat2")),
      # ~ .x %in% intersecting_feats$feature
      ~ .x %in% all.int.comps
    )
  ) %>%
  dplyr:::mutate(pair = purrr::map2(feat1, feat2, ~ c(.x, .y))) %>%
  dplyr::filter(grepl("1 x", name))

message("Writing filtered predictions to table...")
writexl::write_xlsx(
  x = filt_match_diffs %>%
    dplyr::select(-c("pair")),
  path = paste0(res_folder, "/tables/filt_matched_diffs.xlsx"),
  col_names = TRUE,
  format_headers = TRUE,
  use_zip64 = FALSE
)

# ==============================================================================
# m/z predictions subset -------------------------------------------------------
# ==============================================================================
# Checking specifically for the glycoside anad aglycone m/zs
glycoside2 <- MetaboCoreUtils::mass2mz(
  # Rdisop::getMass(Rdisop::getMolecule(glycoside_form)),
  MetaboCoreUtils::calculateMass(glycoside_form)[[1]],
  adduct = MetaboCoreUtils::adducts(polarity = polarity)) %>%
  t() %>%
  tibble::as_tibble(., rownames = "adduct") %>%
  dplyr::rename("glycoside" = V1)

aglycone <- MetaboCoreUtils::mass2mz(
  # Rdisop::getMass(Rdisop::getMolecule(aglycone_form)),
  MetaboCoreUtils::calculateMass(aglycone_form)[[1]],
  adduct = MetaboCoreUtils::adducts(polarity = polarity)) %>%
  t() %>%
  tibble::as_tibble(., rownames = "adduct") %>%
  dplyr::rename("aglycone" = V1)

range_tol <- ppm_to_num(glycoside_ppm)

gly_agly_adducts <- glycoside %>%
  dplyr::left_join(
    x = .,
    y = aglycone,
    by = "adduct"
  ) %>%
  dplyr::mutate(
    glycoside_min = glycoside - range_tol,
    glycoside_max = glycoside + range_tol,
    aglycone_min = aglycone - range_tol,
    aglycone_max = aglycone + range_tol
  )

gly_agly <- tibble::tibble()
for (i in seq_along(gly_agly_adducts)) {
  tmp <- full_raw_filled %>%
    dplyr::filter(
      dplyr::between(
        mzmed,
        gly_agly_adducts[i, ]$glycoside_min,
        gly_agly_adducts[i, ]$glycoside_max
      ) |
        dplyr::between(
          mzmed,
          gly_agly_adducts[i, ]$aglycone_min,
          gly_agly_adducts[i, ]$aglycone_max
      )
    ) %>%
    dplyr::mutate(adduct = gly_agly_adducts[i, ]$adduct) %>%
    dplyr::relocate(adduct, .after = "feature")

  gly_agly <- bind_rows(gly_agly, tmp)
}

pot_glycosides <- unique(gly_agly$feature)
# Only for testing right now
pot_glycosides <- c("FT02089", "FT08181", "FT08191")
matched_diffs2 <- pred_biot_subset(
  data = possible_adducts,
  biotransf_data = bio_transf2, # bio_transf
  tolerance_ppm = ppm_global, # glycoside_ppm
  feat_filt = pot_glycosides
)

filt_match_diffs2 <- matched_diffs2 %>%
  dplyr::rowwise() %>%
  dplyr::mutate(pair = list(c(feat1, feat2))) %>%
  dplyr::ungroup()

glycoside_pairs <- unique(c(filt_match_diffs2$feat1, filt_match_diffs2$feat2))

# ==============================================================================
# Matching m/z's against databases -------------------------------------------
# ==============================================================================
int_mets <- c(filt_match_diffs$feat1, filt_match_diffs$feat2)

peaks_used <- full_norm_filled %>%
  dplyr::select(feature, mzmed, rtmed) %>%
  dplyr::rename("mz" = "mzmed", "rtime" = "rtmed") %>%
  dplyr::filter(feature %in% int_mets) %>%
  tibble::column_to_rownames(var = "feature")

peaks_used$peak_id <- rownames(peaks_used) # keep XCMS peak IDs

annotation_hub <- AnnotationHub()
# query(annotation_hub, "CompDb")
cdb <- annotation_hub[["AH119519"]]

# dbname  <- "CompDb.Hsapiens.HMDB.5.0.sqlite"
# db_file <- file.path("annotation_databases", dbname) # temp.dir()
# if (!file.exists(db_file)) {
#     curl::curl_download(
#         url = paste0(
#             "https://github.com/jorainer/MetaboAnnotationTutorials/",
#             "releases/download/2021-11-02/", dbname
#         ),
#         destfile = db_file
#     )
# }
# cdb <- CompoundDb::CompDb(db_file)

target_df <- ProtGenerics::compounds(
  cdb,
  columns = c(
    "compound_id",
    "name",
    "formula",
    "exactmass",
    "smiles",
    "inchi",
    "inchikey",
    "cas",
    "pubchem"
  )
)

# parameters to match by
mz_match_param <- MetaboAnnotation::Mass2MzParam(
  adducts = c(MetaboCoreUtils::adductNames(polarity = polarity)),
  ppm = ppm_global
)

matches <- MetaboAnnotation::matchValues(
  query = peaks_used,
  target = target_df,
  param = mz_match_param
)

anno <- MetaboAnnotation::matchedData(matches) %>%
  tibble::as_tibble(., rownames = "feature") %>%
  dplyr::mutate(abs_score = abs(score)) %>%
  dplyr::arrange(abs_score)

# ==============================================================================
# Filtering chromatograms ------------------------------------------------------
# ==============================================================================
message(sprintf(
  "Filtering features with sn: %s, beta_cor: %s, beta_snr: %s",
  sn_threshold,
  beta_cor_threshold,
  beta_snr_threshold
))

if (check_saved("xchr9_filt.rds")) {
  xchr9_filt <- readRDS(file = paste0(res_folder, "/objects/xchr9_filt.rds"))
} else {
  xchr9_filt <- filt_features(
    object = xchr9,
    sn_threshold = sn_threshold,
    beta_cor_threshold = beta_cor_threshold,
    beta_snr_threshold = beta_snr_threshold,
    filt_vector = all.int.comps
  )
  saveRDS(
    object = xchr9_filt,
    file = paste0(res_folder, "/objects/xchr9_filt.rds")
  )
}

xchr9_filt$final.plotting.features <- unique(c(
  glycoside_pairs,
  xchr9_filt$final.plotting.features
))

# ==============================================================================
# Plotting features ------------------------------------------------------------
# ==============================================================================
message("Producing feature chromatograms...")
if (check_saved("feature_chrs.rds")) {
  feature_chrs <- readRDS(
    file = paste0(res_folder, "/objects/feature_chrs.rds")
  )
} else {
  feature_chrs <- xcms::featureChromatograms(
    object = xchr9,
    expandRt = 0,
    expandMz = 0,
    aggregationFun = "sum",
    filled = TRUE,
    features = xchr9_filt$final.plotting.features,
    # features = rownames(xcms::featureDefinitions(xchr9)),
    missing = 0,
    return.type = "XChromatograms"
  )
  saveRDS(
    object = feature_chrs,
    file = paste0(res_folder, "/objects/feature_chrs.rds")
  )
}

# # Note
# # The EIC data of a feature is extracted from every sample using the
# # same m/z - rt area. The EIC in a sample does thus not exactly represent the
# # signal of the actually identified chromatographic peak in that sample.
# # The chromPeakChromatograms() function would allow to extract the actual EIC
# # of the chromatographic peak in a specific sample. See also examples below.

# message("Writing feature chromatograms to plots...")
# # Features
# plot_chrom_intensity(
#   chromatogram = feature_chrs,
#   chrom_object = xchr9,
#   save_loc = "/graphs/features/",
#   # amount = 1,
#   peaks_or_feats = "features"
# )

# message("Writing all gap filled peaks to plots...")
# # Gap filled peaks only
# plot_chrom_intensity(
#   chromatogram = chrs_na,
#   chrom_object = xchr8,
#   save_loc = "/graphs/filled_peaks/",
#   # amount = 1,
#   peaks_or_feats = "peaks"
# )

# message("Writing feature chromatograms and intensity boxplots...")
# feats_to_plot <- rownames(xcms::featureDefinitions(feature_chrs))
# feats_to_plot <- xchr9_filt$filt.sig.features
# for (i in feats_to_plot) {
#   ft_p <- plot_feat_chrom_int(
#     feature_chrom = feature_chrs,
#     feature = i,
#     method = "sum",
#     value = "into",
#     filled = TRUE,
#     missing = 0,
#     ms_level = 1,
#     save_loc = "/graphs/feature_chromatogram_intensity/",
#     device = "pdf",
#     feat_pairs = FALSE
#   )
# }

# message("Plotting feature pairs in filtered biotransformations...")
# for (i in seq_len(nrow(xchr9_filt$biot_filt_sig_features_tib))) {
#   ft_pair_p <- plot_feature_pairs(
#     feature_chrom = feature_chrs,
#     filt_match_row = xchr9_filt$biot_filt_sig_features_tib[i, ],
#     method = "sum",
#     value = "into",
#     filled = TRUE,
#     missing = 0,
#     ms_level = 1,
#     save_pairs_loc = "/graphs/feature_pairs/",
#     device = "pdf"
#   )
# }

# message(
#   "Writing feature chromatograms and intensity boxplots "
#   "for glycosides/aglycones..."
# )
# for (i in pot_glycosides) {
#   ft_p <- plot_feat_chrom_int(
#     feature_chrom = feature_chrs,
#     feature = i,
#     method = "sum",
#     value = "into",
#     filled = TRUE,
#     missing = 0,
#     ms_level = 1,
#     save_loc = "/graphs/glycoside/",
#     device = "pdf",
#     feat_pairs = FALSE
#   )
# }

# message("Plotting glycoside/aglycone feature pairs biotransformations...")
# for (i in seq_len(nrow(filt_match_diffs2))) {
#   ft_pair_p <- plot_feature_pairs(
#     feature_chrom = feature_chrs,
#     filt_match_row = filt_match_diffs2[i, ],
#     method = "sum",
#     value = "into",
#     filled = TRUE,
#     missing = 0,
#     ms_level = 1,
#     save_pairs_loc = "/graphs/glycoside_feature_pairs/",
#     device = "pdf"
#   )
# }
