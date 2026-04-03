# ==============================================================================
# Source functions and minimal startup parameters ------------------------------
# ==============================================================================
source("scripts/functions.R")
source("analyses/standard_disco_args.R")
suppressWarnings(
  suppressPackageStartupMessages({
    library(tidyverse)
    library(xcms)
    library(MsExperiment)
    library(optparse)
    library(writexl)
    library(BiocParallel)
    library(MetaboCoreUtils)
    library(MsCoreUtils)
    library(Rdisop)
  })
)

dir.create(file.path(opt$output, "peaks"), FALSE, TRUE)

stds <- import_mzml(data_path = opt$data_path, meta_file = opt$meta_file)

if (file.exists(file.path(opt$output, "stds_exp.rds"))) {
  stds_exp <- readRDS(file.path(opt$output, "stds_exp.rds"))
} else {
  stds_exp <- MsExperiment::readMsExperiment(
    spectraFiles = stds$path,
    sampleData = stds
  )
  saveRDS(stds_exp, file.path(opt$output, "stds_exp.rds"))
}

# ==============================================================================
# Source functions and minimal startup parameters ------------------------------
# ==============================================================================
std_mz_theory <- get_theory_mz(
  chem_form = opt$is_formula,
  adduct = opt$is_adduct
)
std_mz_range <- get_short_mz_range(std_mz_theory, mz_window = 0.02)

if (file.exists(file.path(opt$output, "std_chr_wide.rds"))) {
  std_chr <- readRDS(file.path(opt$output, "std_chr_wide.rds"))
} else {
  std_chr <- xcms::chromatogram(
    object = stds_exp,
    BPPARAM = BiocParallel::MulticoreParam(opt$cores),
    chunkSize = opt$cores,
    mz = std_mz_range,
    aggregationFun = "sum"
  )
  saveRDS(std_chr, file.path(opt$output, "std_chr_wide.rds"))
}

# Individual IS XICs
# SAVE THEM HERE ONLY
for (i in seq_along(std_chr)) {
  samp_idx <- which(colnames(std_chr)[i] == rownames(stds))
  plot(
    x = std_chr[, i],
    lwd = 3,
    main = paste0(
      stds$group[samp_idx], "\n",
      colnames(std_chr)[i]
    )
  )
  # SAVE THEM HERE ONLY
}
invisible(dev.off())

colnames(std_chr)[24]

# GGPLOT WAY OF PLOTTING ALL CHROMATOGRAMS
tibble(
  sample = colnames(std_chr)[24],
  rt = std_chr[, 24]@rtime,
  intensity = std_chr[, 24]@intensity,
  mzmin = std_chr[, 24]@mz[1],
  mzmax = std_chr[, 24]@mz[2]
) %>%
  dplyr::mutate(intensity = tidyr::replace_na(intensity, replace = 0)) %>%
  ggplot(
    aes(
      x = rt,
      y = intensity
    )
  ) +
  geom_line() +
  theme_classic() +
  scale_y_continuous(expand = expansion(c(0, 0.05)))
# ENd of plotting of all graphds

# ==============================================================================
# Source functions and minimal startup parameters ------------------------------
# ==============================================================================
std_ranges <- get_rt_mz_range(chromatogram = std_chr, rt_window = 0.05)
# Get the IS XIC
if (file.exists(file.path(opt$output, "std_eic_wide.rds"))) {
  std_eic_wide <- readRDS(file.path(opt$output, "std_eic_wide.rds"))
} else {
  std_eic_wide <- xcms::chromatogram(
    object = stds_exp,
    BPPARAM = BiocParallel::MulticoreParam(opt$cores),
    chunkSize = opt$cores,
    mz = std_ranges$mz_range,
    rt = std_ranges$rt_range,
    aggregationFun = "sum"
  )
  saveRDS(std_eic_wide, file.path(opt$output, "std_eic_wide.rds"))
}

# plotting chromatograms
tibble(
  sample = colnames(std_eic_wide)[24],
  rt = std_eic_wide[, 24]@rtime,
  intensity = std_eic_wide[, 24]@intensity,
  mzmin = std_eic_wide[, 24]@mz[1],
  mzmax = std_eic_wide[, 24]@mz[2]
) %>%
  dplyr::mutate(intensity = tidyr::replace_na(intensity, replace = 0)) %>%
  ggplot(
    aes(
      x = rt,
      y = intensity
    )
  ) +
  geom_line() +
  theme_classic() +
  scale_y_continuous(expand = expansion(c(0, 0.05)))
# end of plotting ghromatograms

# ==============================================================================
# Finding peaks ------------------------------
# ==============================================================================
if (file.exists(file.path(opt$output, "std_peaks.rds"))) {
  std_peaks <- readRDS(file.path(opt$output, "std_peaks.rds"))
} else {
  # Run peak detection on the EIC
  std_peaks <- xcms::findChromPeaks(
    object = stds_exp,
    BPPARAM = BiocParallel::MulticoreParam(opt$cores),
    chunkSize = opt$cores,
    msLevel = 1L,
    param = xcms::CentWaveParam(
      ppm = opt$ppm_global,
      peakwidth = c(1, 50),
      prefilter = c(1, 1),
      snthresh = 100, # set
      mzCenterFun = "wMean",
      mzdiff = 0.01, # set
      integrate = 2,
      noise = 5000, # set
      verboseBetaColumns = FALSE
    )
  )
  saveRDS(std_peaks, file.path(opt$output, "std_peaks.rds"))
}

grouped_peaks <- xcms::groupChromPeaks(
  object = std_peaks,
  param = xcms::PeakDensityParam(
    sampleGroups = MsExperiment::sampleData(std_peaks)$group,
    bw = 0.2, # 0.5
    minFraction = 0.5,
    binSize = 0.01,
    maxFeatures = 1000, # 200
    ppm = opt$ppm_global,
    minSamples = 2 # 1
  )
)

map_samples <- tibble::tibble(
  sample = seq_along(MSnbase::fileNames(grouped_peaks)),
  filename = basename(MSnbase::fileNames(grouped_peaks))
) %>%
  dplyr::left_join(
    x = .,
    y = tibble::as_tibble(stds, rownames = "filename"),
    by = "filename"
  ) %>%
  tidyr::separate_rows(formula, group, sep = ", ") %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    mass = Rdisop::getMonoisotopic(Rdisop::getMolecule(formula))
  ) %>%
  dplyr::ungroup()

# Expected m/z for all negative adducts, one row per formula × adduct
expected_mz <- MetaboCoreUtils::mass2mz(
  x = setNames(
    dplyr::distinct(map_samples, formula, mass)$mass,
    dplyr::distinct(map_samples, formula, mass)$formula
  ),
  adduct = MetaboCoreUtils::adducts(polarity = "negative")
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("formula") %>%
  tidyr::pivot_longer(
    cols = -formula,
    names_to = "adduct",
    values_to = "mz_theory"
  ) %>%
  dplyr::filter(!is.na(mz_theory))

formula_to_group <- map_samples %>%
  dplyr::distinct(formula, group)

# Map sample files to their original (combined) group names — one row per file
map_samples_orig <- tibble::tibble(
  filename = basename(MSnbase::fileNames(grouped_peaks))
) %>%
  dplyr::left_join(
    y = tibble::as_tibble(stds, rownames = "filename") %>%
      dplyr::select(filename, group),
    by = "filename"
  )

# Map chromPeak index → feature using peakidx from featureDefinitions
feat_peak_map <- xcms::featureDefinitions(grouped_peaks) %>%
  tibble::as_tibble(rownames = "feature") %>%
  dplyr::select(feature, peakidx) %>%
  tidyr::unnest(cols = peakidx)

# sample index → original group
sample_group_map <- tibble::tibble(
  sample = seq_along(MSnbase::fileNames(grouped_peaks)),
  filename = basename(MSnbase::fileNames(grouped_peaks))
) %>%
  dplyr::left_join(y = map_samples_orig, by = "filename")

# Assign each feature to the original sample group with highest mean peak
# intensity, then split combined group names to individual compounds
feat_in_compound_samples <- xcms::chromPeaks(grouped_peaks) %>%
  tibble::as_tibble(rownames = "peak") %>%
  dplyr::mutate(peak_idx = dplyr::row_number()) %>%
  dplyr::left_join(y = feat_peak_map, by = c("peak_idx" = "peakidx")) %>%
  dplyr::left_join(y = sample_group_map, by = "sample") %>%
  dplyr::group_by(feature, group) %>%
  dplyr::summarise(mean_into = mean(into, na.rm = TRUE), .groups = "drop") %>%
  dplyr::group_by(feature) %>%
  dplyr::slice_max(mean_into, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  tidyr::separate_rows(group, sep = ", ") %>%
  dplyr::select(feature, group)

# Match features to standards: mz within ppm AND detected in correct samples
std_feature_matches <- xcms::featureDefinitions(grouped_peaks) %>%
  tibble::as_tibble(rownames = "feature") %>%
  dplyr::select(feature, mzmed, mzmin, mzmax, rtmed, rtmin, rtmax, npeaks) %>%
  dplyr::cross_join(y = expected_mz) %>%
  dplyr::filter(abs(mzmed - mz_theory) / mz_theory * 1e6 <= opt$ppm_global) %>%
  dplyr::left_join(y = formula_to_group, by = "formula") %>%
  dplyr::inner_join(
    x = .,
    y = feat_in_compound_samples,
    by = c("feature", "group")
  ) %>%
  dplyr::arrange(group, adduct)

feat_values_by_group <- xcms::featureValues (
  object = grouped_peaks,
  method = "sum",
  missing = 0,
  value = "into"
) %>%
  tibble::as_tibble(rownames = "feature") %>%
  tidyr::pivot_longer(
    cols = -feature,
    names_to = "filename",
    values_to = "into"
  ) %>%
  dplyr::inner_join(
    y = dplyr::distinct(map_samples, filename, group),
    by = "filename"
  ) %>%
  dplyr::group_by(feature, group) %>%
  dplyr::summarise(into = sum(into, na.rm = TRUE), .groups = "drop")

final_standards <- std_feature_matches %>%
  dplyr::left_join(y = feat_values_by_group, by = c("feature", "group")) %>%
  dplyr::arrange(group, dplyr::desc(into))

best_stds <- final_standards %>%
  dplyr::slice_max(into, n = 1, by = group)

feat_group_map <- final_standards %>%
  dplyr::slice_max(into, n = 1, by = group) %>%
  dplyr::select(feature, group, adduct, mz_theory, rtmed)

feats <- feat_group_map$feature

test <- xcms::featureChromatograms(
  object = grouped_peaks,
  features = feats
)

tester <- purrr::map_dfr(seq_len(nrow(test)), function(i) {
  purrr::map_dfr(seq_len(ncol(test)), function(j) {
    tibble::tibble(
      rt        = rtime(test[i, j]),
      intensity = intensity(test[i, j]),
      filename  = colnames(test)[j],
      feature   = feats[i]
    )
  })
}) %>%
  dplyr::mutate(intensity = tidyr::replace_na(intensity, 0)) %>%
  dplyr::left_join(y = map_samples_orig, by = "filename") %>%
  dplyr::left_join(
    y = dplyr::distinct(feat_group_map, feature, .keep_all = TRUE),
    by = "feature"
  ) %>%
  dplyr::filter(stringr::str_detect(group.x, group.y)) %>%
  dplyr::select(-group.x) %>%
  dplyr::rename(group = group.y)


tester %>%
  dplyr::filter(adduct == "[M-H]-") %>%
  ggplot(
    aes(
      x = rt,
      y = intensity,
      color = group,
      group = filename
    )
  ) +
  geom_line() +
  facet_wrap(~ group, scales = "free") +
  theme_classic()

## no duplicates anymore
final_standards %>%
  dplyr::group_by(feature) %>%
  dplyr::filter(dplyr::n_distinct(group) > 1) %>%
  dplyr::distinct(feature, group, mzmed, rtmed) %>%
  dplyr::arrange(feature)


final_peaks <- std_peak_matches %>%
  dplyr::group_by(group, adduct) %>%
  dplyr::slice_max(into, n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(group, dplyr::desc(into))



# ==============================================================================
# Scratch / exploratory --------------------------------------------------------
# ==============================================================================

# Filter detected peaks to only those matching a standard in that sample
std_peak_matches <- xcms::chromPeaks(std_peaks) %>%
  tibble::as_tibble(rownames = "peak") %>%
  dplyr::left_join(
    dplyr::left_join(map_samples, expected_mz, by = "formula"),
    by = "sample",
    relationship = "many-to-many"
  ) %>%
  dplyr::filter(abs(mz - mz_theory) / mz_theory * 1e6 <= opt$ppm_global) %>%
  dplyr::arrange(group, dplyr::desc(into)) %>%
  dplyr::select(-path)

std_peak_matches %>%
  dplyr::arrange(dplyr::desc(into))

std_peak_matches %>%
  dplyr::filter(group == "Apigenin") %>%
  dplyr::pull(peak)

# Plotting peaks from peak list
p <- xcms::chromPeaks(std_peaks)["CP05290", ]

test <- xcms::chromatogram(
  object = xcms::filterFile(std_peaks, p[["sample"]]),
  mz = c(p[["mzmin"]], p[["mzmax"]]),
  rt = c(p[["rtmin"]] - 16, p[["rtmax"]] + 16)
)

tibble::tibble(
  rt = test[1, 1]@rtime,
  int = test[1, 1]@intensity,
  filename = colnames(test)
) %>%
  dplyr::mutate(int = tidyr::replace_na(int, replace = 0)) %>%
  dplyr::left_join(
    x = .,
    y = map_samples,
    by = "filename"
  ) %>%
  ggplot(
    aes(
      x = rt,
      y = int,
      color = group
    )
  ) +
  geom_line() +
  theme_classic() +
  scale_y_continuous()
