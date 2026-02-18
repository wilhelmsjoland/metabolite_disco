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
# Inspect internal standard prior to peak-calling ------------------------------
# Define the rt and m/z range of the peak area ---------------------------------
# ==============================================================================

message(
  "===========================================================================",
  "\n",
  "Inspecting internal standard peaks prior to peak-calling ------------------",
  "\n",
  "==========================================================================="
)
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

# All IS XICs together
pdf(paste0(res_folder, "/graphs/internal_standard/all_is.pdf"))
plot(x = is_eic, col = group_colors[is_eic$group], lwd = 3)
legend("topleft", legend = names(group_colors), col = group_colors, pch = 16)
invisible(dev.off())

# Individual IS XICs
for (i in seq_along(is_eic)) {
  pdf(
    paste0(
      res_folder, 
      "/graphs/internal_standard/", 
      colnames(bpcs)[i], 
      ".pdf"
    )
  )

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
  legend("topleft", legend = names(group_colors), col = group_colors, pch = 16)
  invisible(dev.off())
}

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
message(
  "===========================================================================",
  "\n",
  "Calling peaks -------------------------------------------------------------",
  "\n",
  "==========================================================================="
)

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

# Use this for the beta-distribution parameters
xchr_data <- tibble::as_tibble(xcms::chromPeaks(xchr), rownames = "peak")

# ==============================================================================
# Inspect peaks
# TODO
# Extract some peaks here and check quality of peak picking
# ==============================================================================
message(
  "===========================================================================",
  "\n",
  "Inspecting called peaks ---------------------------------------------------",
  "\n",
  "==========================================================================="
)

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
# Alignment of retention times
# ==============================================================================
message(
  "===========================================================================",
  "\n",
  "Aligning retention times---------------------------------------------------",
  "\n",
  "==========================================================================="
)

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
group_factor <- MsExperiment::sampleData(xchr8)$grRoup
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
message(
  "===========================================================================",
  "\n",
  "Generating PCAs -----------------------------------------------------------",
  "\n",
  "==========================================================================="
)

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

# Data before normalization
vals <- SummarizedExperiment::assay(res, "raw_filled") %>%
  log2() %>%
  t() %>%
  scale(center = TRUE, scale = TRUE) %>%
  as.matrix(.)

pca_res <- prcomp(vals, scale = FALSE, center = FALSE)
pca_raw <- plot_pca(
  prcomp_res = pca_res,
  metad = meta,
  x = PC1,
  y = PC2
) +
  ggplot2::labs(title = "Before normalization")

# Data after normalization
vals_norm <- SummarizedExperiment::assay(res, "norm_filled") %>%
  log2() %>%
  t() %>%
  scale(center = TRUE, scale = TRUE) %>%
  as.matrix(.)

pca_res_norm <- prcomp(vals_norm, scale = FALSE, center = FALSE)
pca_adj <- plot_pca(
  prcomp_res = pca_res_norm,
  metad = meta,
  x = PC1,
  y = PC2
) +
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
pca_raw <- plot_pca(
  prcomp_res = pca_res,
  metad = meta,
  x = PC3,
  y = PC4
) +
  ggplot2::labs(title = "Before normalization")

pca_adj <- plot_pca(
  prcomp_res = pca_res_norm,
  metad = meta,
  x = PC3,
  y = PC4
) +
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
message(
  "===========================================================================",
  "\n",
  "Running linear models -----------------------------------------------------",
  "\n",
  "==========================================================================="
)

group_used <- factor(meta$group)
design <- model.matrix(~ 0 + group_used)
colnames(design) <- levels(group_used)

comparisons <- combn(
  x = levels(group_used),
  m = 2,
  simplify = TRUE) %>%
  t(.) %>%
  tibble::as_tibble(., rownames = "rownumber") %>%
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
message(
  "===========================================================================",
  "\n",
  "Generating upset plots ----------------------------------------------------",
  "\n",
  "==========================================================================="
)

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
message(
  "===========================================================================",
  "\n",
  "Finding intersecting features ---------------------------------------------",
  "\n",
  "==========================================================================="
)

# Create all comparisons
combinations <- unlist(
  lapply(
    X = seq_along(comparisons),
    FUN = function(x) combn(comparisons, x, simplify = FALSE)
  ),
  recursive = FALSE
)

upset_comps <- list()
for (i in seq_along(combinations)) {
  tmp_intersect_feats <- find_intersect_feat(
    data = upset_tib,
    set = combinations[[i]],
    full_set = comparisons
  )

  upset_comps[[
    stringr::str_flatten(combinations[[i]], collapse = "*")
  ]] <- tmp_intersect_feats$feature
}

# Only temporary
upset_comp <- upset_comps[[
  paste0(
    "bu_mutant_apiin-bu_wt_apiin*bu_mutant_control",
    "-",
    "bu_wt_apiin*bu_wt_apiin-bu_wt_control"
  )
]]

# ==============================================================================
# Prepare biotransformation checking -------------------------------------------
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
  data <- readr::read_tsv(
    "scripts/search_compounds/output/rpairs.tsv",
    show_col_types = FALSE
  ) %>%
  dplyr::mutate(rpair_num = paste0("RP", dplyr::row_number()))

  df <- data %>%
    tidyr::drop_na(formula1, formula2) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      dplyr::across(
        .cols = dplyr::all_of(c("formula1", "formula2")),
        .fns = ~ MetaboCoreUtils::standardizeFormula(.)
      )
    ) %>%
    dplyr::mutate(
      delta_formula = {
        # Try subtracting both directions
        diff1 <- MetaboCoreUtils::subtractElements(formula1, formula2)
        diff2 <- MetaboCoreUtils::subtractElements(formula2, formula1)
        # Use whichever worked
        if (!is.na(diff1)) {
          diff1
        } else if (!is.na(diff2)) {
          diff2
        } else {
          paste0("\u00B1 ", formula1, " <=> ", formula2)
        }
      }
    ) %>%
    dplyr::mutate(delta_mass = Rdisop::getMonoisotopic(delta_formula)) %>%
    dplyr::ungroup() %>%
    dplyr::relocate(rpair_num, .before = "entry") %>%
    dplyr::mutate(
      across(
        .cols = tidyselect::all_of(c("name1", "name2")),
        .fns = ~ stringr::str_trim(stringr::str_extract(., "^[^;]+"))
      )
    ) %>%
    dplyr::mutate(allowed_n = 1) %>%
    tidyr::uncount(
      data = .,
      weights = allowed_n,
      .id = "multiplier",
      .remove = FALSE
    ) %>%
    dplyr::select(rpair_num, delta_formula, allowed_n, multiplier, delta_mass)

  biotransf_append <- df %>%
    dplyr::group_by(delta_mass) %>%
    dplyr::summarize(
      delta_formula = paste(sort(unique(delta_formula)), collapse = ", "),
      rpair_nums = paste(sort(unique(rpair_num)), collapse = ", "),
      n_rpairs   = dplyr::n_distinct(rpair_num),
      allowed_n  = dplyr::first(allowed_n),
      multiplier = dplyr::first(multiplier),
      .groups = "drop"
    ) %>%
    dplyr::select(
      name = rpair_nums,
      delta_formula,
      allowed_n,
      multiplier,
      delta_mass
    ) %>%
    dplyr::filter(!grepl("<=>", delta_formula)) %>%
    # nothing changed? or nothing has changed at least according
    # since the delta formulas probably only had 1n or 1c or similar and were
    # subtracted
    tidyr::drop_na(delta_mass)
} else {
  message("'biotransf_append' already exists.")
}

bio_transf2 <- dplyr::bind_rows(
  bio_transf,
  biotransf_append
)

# Potentially filter noisy features with the filt_features() function I made
all_sig_diff <- sort(unique(unlist(upset_comps, use.names = FALSE)))
possible_adducts_signif <- possible_adducts %>%
  dplyr::filter(feature %in% all_sig_diff)

# ==============================================================================
# Filtering features -----------------------------------------------------------
# ==============================================================================
message(sprintf(
  "Filtering features with sn: %s, beta_cor: %s, beta_snr: %s",
  sn_threshold,
  beta_cor_threshold,
  beta_snr_threshold
))

if (file.exists(file.path(res_folder, "objects", "xchr9_filt.rds"))) {
  xchr9_filt <- readRDS(file.path(res_folder, "objects", "xchr9_filt.rds"))
} else {
  xchr9_filt <- filt_features(
    object = xchr9,
    sn_threshold = sn_threshold,
    beta_cor_threshold = beta_cor_threshold,
    beta_snr_threshold = beta_snr_threshold
  )
  saveRDS(
    object = xchr9_filt,
    file = paste0(res_folder, "/objects/xchr9_filt.rds")
  )
}

# ==============================================================================
# m/z predictions subset -------------------------------------------------------
# ==============================================================================
message(
  "Predicting potential biotransformations based on:\n\t",
  "Biotransformation database: ", biotransf_file, # + the other kegg stuff
  "\n\tppm: ", ppm_match,
  sep = ""
)

# Checking specifically for the glycoside anad aglycone m/zs
glycoside <- MetaboCoreUtils::mass2mz(
  MetaboCoreUtils::calculateMass(glycoside_form)[[1]],
  adduct = MetaboCoreUtils::adducts(polarity = polarity)) %>%
  t() %>%
  tibble::as_tibble(., rownames = "adduct") %>%
  dplyr::rename("glycoside" = V1)

aglycone <- MetaboCoreUtils::mass2mz(
  MetaboCoreUtils::calculateMass(aglycone_form)[[1]],
  adduct = MetaboCoreUtils::adducts(polarity = polarity)) %>%
  t() %>%
  tibble::as_tibble(., rownames = "adduct") %>%
  dplyr::rename("aglycone" = V1)

# This is okay for now since it's only looking for the glycone and aglycone
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
subset_matched_diffs <- pred_biot(
  data = possible_adducts_signif,
  biotransf_data = bio_transf2, # bio_transf
  tolerance_ppm = ppm_match, # glycoside_ppm
  features_of_interest = pot_glycosides,
  parallel = TRUE
) %>%
  dplyr::mutate(
    pair = purrr::map2(
      .x = feat1,
      .y = feat2,
      .f = c
    ),
    obs_diff = abs(obs_delta_mass - delta_mass)
  )

glycoside_pairs <- unique(
  c(
    subset_matched_diffs$feat1,
    subset_matched_diffs$feat2
  )
)

# ==============================================================================
# m/z predictions all ----------------------------------------------------------
# ==============================================================================

# Turn this around
if (file.exists(file.path(res_folder, "objects", "matched_diffs.rds"))) {
  matched_diffs <- readRDS(
    file = file.path(res_folder, "objects", "matched_diffs.rds")
  )
} else {
  matched_diffs <- pred_biot(
    data = possible_adducts_signif,
    biotransf_data = bio_transf2,
    tolerance_ppm = ppm_match, # try 5 and 10, # 15 too much
    parallel = TRUE
  )
  saveRDS(
    object = matched_diffs,
    file = paste0(res_folder, "/objects/matched_diffs.rds")
  )
}

matched_diffs2 <- matched_diffs %>%
  dplyr::filter(
    dplyr::if_all(
      .cols = dplyr::all_of(c("mass1", "mass2")),
      .fns = ~ . > 0
    )
  ) %>%
  dplyr::filter(feat1 != feat2)
  # Too slow just do for a few when filtered
  # dplyr::mutate(
  #   pair = purrr::map2(
  #     .x = feat1,
  #     .y = feat2,
  #     .f = c
  #   ),
  #   obs_diff = abs(obs_delta_mass - delta_mass)
  # )

message("Writing predictions to table...")
# TODO This needs filtering first
# readr::write_csv(
#   x = matched_diffs2,
#   file = file.path(res_folder, "tables", "matched_diffs.csv")
# )

# ==============================================================================
# Matching m/z's against databases -------------------------------------------
# ==============================================================================
peaks_used <- full_norm_filled %>%
  dplyr::select(feature, mzmed, rtmed) %>%
  dplyr::rename("mz" = "mzmed", "rtime" = "rtmed") %>%
  # All comps with at least one significant difference
  dplyr::filter(feature %in% all_sig_diff) %>%
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
  ppm = ppm_match
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
# Biotransformer.jar -----------------------------------------------------------
# ==============================================================================
if (!file.exists(file.path(getwd(), res_folder, "tables", "prediction.csv"))) {
  run_biotransformer(
    bt_dir = biotransformer_path, # "biotransformer3.0jar",
    # Needs optparse smiles argument
    smiles = opt$smiles,
    b_type = "superbio",
    k_task = "pred",
    output_file = "prediction"
  )
} else {
  biot_pred <- readr::read_csv(
    file = file.path(res_folder, "tables", "prediction.csv"),
    show_col_types = FALSE
  )
}

biot_dedup <- biot_pred %>%
  dplyr::group_by(InChIKey) %>%
  dplyr::summarize(
    dplyr::across(
      .cols = setdiff(colnames(.), "InChIKey"),
      .fns  = ~ paste(unique(.x), collapse = ", ")
    ),
    .groups = "keep"
  )

biot_mass <- biot_dedup %>%
  dplyr::mutate(
    mass = MetaboCoreUtils::calculateMass(`Molecular formula`)
  ) %>%
  dplyr::relocate(mass, .before = "InChI")

biot_mets <- biot_mass$mass
names(biot_mets) <- biot_mass$InChIKey

biot_final <- MetaboCoreUtils::mass2mz(
  x = biot_mets,
  adduct = MetaboCoreUtils::adducts(polarity = polarity)
) %>%
  tibble::as_tibble(., rownames = "InChIKey") %>%
  tidyr::pivot_longer(
    cols = 2:ncol(.),
    names_to = "adduct",
    values_to = "mz"
  ) %>%
  dplyr::arrange(InChIKey) %>%
  dplyr::left_join(
    x = .,
    y = biot_mass,
    by = "InChIKey",
    relationship = "many-to-one"
  )

biot_mass_len <- length(biot_mass$InChIKey) *
  nrow(adducts(polarity = polarity))

if (biot_mass_len != nrow(biot_final)) {
  warning("The transformation prediction dataframes are not the same length.")
} else if (biot_mass_len == nrow(biot_final)) {
  message("The transformation prediction dataframes are the same length.")
}

def_tib <- xchr9_filt

# biot_final = mass to mzs - > match the m/zs to the m/zs in the data
predicted_feats <- biot_final %>%
  dplyr::inner_join(
    x = .,
    y = def_tib %>%
      dplyr::mutate(
        tol = MsCoreUtils::ppm(mzmed, ppm_match),
        mz_lo = mzmed - tol,
        mz_hi = mzmed + tol
      ) %>%
      dplyr::relocate(c("mz_lo", "mz_hi"), .after = "feature"),
    by = dplyr::join_by(dplyr::between(mz, mz_lo, mz_hi))
  )

pred_peak_ids <- sort(unique(predicted_feats$feature))

if (file.exists(file.path(res_folder, "objects", "pred_chrs.rds"))) {
  pred_chrs <- readRDS(file.path(res_folder, "objects", "pred_chrs.rds"))
} else {
  pred_chrs <- xcms::featureChromatograms(
    object = xchr9,
    expandRt = 0,
    expandMz = 0,
    aggregationFun = "sum",
    filled = TRUE,
    features = pred_peak_ids,
    missing = 0,
    return.type = "XChromatograms"
  )
  saveRDS(
    object = pred_chrs,
    file = file.path(res_folder, "objects", "pred_chrs.rds")
  )
}

# ==============================================================================
# Molecular similarity m/z matching --------------------------------------------
# ==============================================================================
anno_filt <- anno %>%
  # TODO
  # DO A more intelligent filtering than this based on
  # how much information is available in all the rows
  dplyr::group_by(adduct) %>% # Think this might have fixed the loss of same
  dplyr::distinct(target_inchikey, .keep_all = TRUE)

# TODO
# ==============================================================================
# # For checking
# # These are all the base feature without "." added
# anno %>%
#   dplyr::filter(is.na(adduct)) %>%
#   print(n = 300)

# # looks good
# anno %>%
#   dplyr::group_by(peak_id, adduct) %>%
#   dplyr::distinct(target_inchikey, .keep_all = TRUE) %>%
#   dplyr::filter(!is.na(adduct)) %>%
#   dplyr::group_by(feature, adduct, target_inchikey) %>%
#   dplyr::summarize(n = n()) %>%
#   dplyr::filter(n != 1)
# ==============================================================================

anno_smiles <- anno_filt$target_smiles
names(anno_smiles) <- anno_filt$feature

anno_sims <- mol_similarity(
  query_smiles = opt$smiles,
  target_smiles = anno_smiles,
  kekulise = TRUE, # parsing incorrect smiles with electrons
  omit_nulls = TRUE,
  fingerprint = "circular",
  circular_type = "ECFP6",
  method = "tanimoto"
) %>%
  dplyr::left_join(
    x = .,
    y = anno,
    by = "feature"
  ) %>%
  # TODO
  # Arbitrary for now
  dplyr::filter(sim > 0.1)

# Get names for individual inchikey
inchi_ks_split <- split(
  anno_sims$target_inchikey,
  ceiling(seq_along(anno_sims$target_inchikey) / 500)
)
inchi_ks <- lapply(
  X = inchi_ks_split,
  FUN = function(x) {
    paste0(x, collapse = ",")
  }
)

full_inchikey_map <- tibble::tibble()
for (i in seq_along(inchi_ks)) {
  cmd <- paste0(
    "curl -s -d \"inchikey=", inchi_ks[[i]], "\" ",
    "\"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/property/",
    "Title,MolecularFormula,InChIKey,ExactMass,XLogP,TPSA/",
    "CSV\""
  )
  results <- system(cmd, intern = TRUE)
  # TODO
  # vroom warning here
  inchikey_map <- readr::read_csv(file = paste0(results, collapse = "\n"))
  full_inchikey_map <- dplyr::bind_rows(full_inchikey_map, inchikey_map)
}

anno_sims2 <- anno_sims %>%
  dplyr::left_join(
    x = .,
    y = full_inchikey_map,
    by = c("target_inchikey" = "InChIKey")
  )

# ==============================================================================
# Molecular similarity biotransformer ------------------------------------------
# ==============================================================================
chem_pred_feats <- predicted_feats %>%
  # This is needed since one inchikey can be annotated to
  # several features because they can have different adducts
  # ---------------------------------------------- #
  dplyr::mutate(met_id = as.character(dplyr::row_number())) %>%
  dplyr::relocate(met_id, .before = "InChIKey") %>%
  # ---------------------------------------------- #

  # Here don't need unique since it's not matched against a db
  # so shouldnt have duplicates for the same feature and adduct
  # ---------------------------------------------- #
  dplyr::group_by(feature, adduct) %>%
  dplyr::distinct(InChIKey, .keep_all = TRUE)
  # ---------------------------------------------- #

pred_smiles <- chem_pred_feats$SMILES
names(pred_smiles) <- chem_pred_feats$met_id

pred_sims <- mol_similarity(
  query_smiles = opt$smiles,
  target_smiles = pred_smiles,
  kekulise = TRUE, # parsing incorrect smiles with electrons
  omit_nulls = TRUE,
  fingerprint = "circular",
  circular_type = "ECFP6",
  method = "tanimoto"
) %>%
  # rename to not overlap in the dplyr::left_join()
  dplyr::rename("met_id" = "feature") %>%
  dplyr::left_join(
    x = .,
    y = chem_pred_feats,
    by = c("met_id")
  )

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
    # TODO Fix
    features = xchr9_filt$final.plotting.features,
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

message("Writing feature chromatograms to plots...")
# # Features
# plot_chrom_intensity(
#   chromatogram = feature_chrs,
#   chrom_object = xchr9,
#   save_loc = "/graphs/features/",
#   # amount = 1,
#   peaks_or_feats = "features"
# )

message("Writing all gap filled peaks to plots...")
# # Gap filled peaks only
# plot_chrom_intensity(
#   chromatogram = chrs_na,
#   chrom_object = xchr8,
#   save_loc = "/graphs/filled_peaks/",
#   # amount = 1,
#   peaks_or_feats = "peaks"
# )

message("Writing feature chromatograms and intensity boxplots...")
# feats_to_plot <- sort(unique(all_sig_diff))
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

message("Plotting feature pairs in filtered biotransformations...")
# TODO
# FIX this
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

message(
  "Writing feature chromatograms and intensity boxplots ",
  "for glycosides/aglycones..."
)
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

message("Plotting glycoside/aglycone feature pairs biotransformations...")
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

message("Producing significant intersecting feature boxplots...")
# for (i in upset_comp) {
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

message("Plotting annotated features after molecular similarity")
# anno_sims2 %>%
#   # TODO
#   # arbitrary for now
#   dplyr::filter(sim > 0.3) %>%
#   dplyr::mutate(
#     Title = forcats::fct_reorder( # target_name
#       .f = Title, # target_name
#       .x = sim,
#       .fun = "mean",
#       .desc = TRUE
#     )
#   ) %>%
#   ggplot2::ggplot(
#     ggplot2::aes(
#       x = Title, # target_name
#       y = sim
#     )
#   ) +
#   # because there are duplicates - just choose the best one
#   ggplot2::geom_col(
#     # aes(fill = peak_id),
#     stat = "summary",
#     fun = "max",
#     color = "black",
#     position = ggplot2::position_dodge()
#   ) +
#   ggplot2::guides(x = ggplot2::guide_axis(angle = -45)) +
#   ggplot2::scale_y_continuous(
#     expand = ggplot2::expansion(c(0, 0)),
#     limits = c(0, 1)
#   ) +
#   ggplot2::theme_bw() +
#   ggplot2::theme(
#     axis.title.x = ggplot2::element_blank(),
#     axis.title.y = ggplot2::element_text(angle = -90),
#     legend.title = ggplot2::element_blank()
#   ) +
#   ggplot2::labs(y = "Tanimoto similarity")

message("Plotting predicted features after molecular similarity")
# pred_sims %>%
#   dplyr::mutate(
#     met_id = forcats::fct_reorder( # target_name
#       .f = met_id, # target_name
#       .x = sim,
#       .fun = "mean",
#       .desc = TRUE
#     )
#   ) %>%
#   ggplot2::ggplot(
#     ggplot2::aes(
#       x = met_id, # target_name
#       y = sim
#     )
#   ) +
#   # because there are duplicates - just choose the best one
#   ggplot2::geom_col(
#     aes(fill = feature),
#     stat = "summary",
#     fun = "max",
#     color = "black",
#     position = ggplot2::position_dodge()
#   ) +
#   ggplot2::guides(x = ggplot2::guide_axis(angle = -45)) +
#   ggplot2::scale_y_continuous(
#     expand = ggplot2::expansion(c(0, 0)),
#     limits = c(0, 1)
#   ) +
#   ggplot2::theme_bw() +
#   ggplot2::theme(
#     axis.title.x = ggplot2::element_blank(),
#     axis.title.y = ggplot2::element_text(angle = -90),
#     legend.title = ggplot2::element_blank()
#   ) +
#   ggplot2::labs(y = "Tanimoto similarity")
