# ==============================================================================
# Source functions and minimal startup parameters ------------------------------
# ==============================================================================
source("scripts/functions.R")
start_log(snakemake@params$output)
script_header()

set.seed(snakemake@params$seed)
register_parallel(snakemake@params$cores)
suppressWarnings(
  suppressPackageStartupMessages(
    {
      library(cli)
      library(BiocParallel)
      library(xcms)
      library(MsExperiment)
      library(dplyr)
      library(tibble)
      library(RSQLite)
      library(MsBackendSql)
    }
  )
)

setup <- readRDS(snakemake@input[["setup"]])
bpc_data <- readRDS(snakemake@input[["bpc"]])
ms_exp <- setup$ms_exp
group_colors <- setup$group_colors
meta <- setup$meta
bpcs <- bpc_data$bpcs

# ==============================================================================
# Inspect internal standard prior to peak-calling ------------------------------
# Define the rt and m/z range of the peak area ---------------------------------
# ==============================================================================
cli::cli_alert_info("Inspecting internal standard peaks prior to peak-calling")
mz_theory <- get_theory_mz(
  chem_form = snakemake@params$internal_standard,
  adduct = snakemake@params$is_adduct
)
mz_range <- get_short_mz_range(mz_theory, mz_window = 0.02)

is_chr_path <- file.path(snakemake@params$output, "objects", "is_chr.rds")
if (interactive() && file.exists(is_chr_path)) {
  is_chr <- readRDS(file = is_chr_path)
  cli::cli_alert_success(
    paste0(
      "Imported saved full internal standard chromatogram from: ",
      " {.path {is_chr_path}}"
    )
  )
} else {
  cli::cli_alert_info("Creating full internal standard chromatogram")
  is_chr <- xcms::chromatogram(
    BPPARAM = BiocParallel::SerialParam(),
    chunkSize = 1L,
    object = ms_exp,
    mz = mz_range,
    aggregationFun = "sum"
  )
  saveRDS(object = is_chr, file = is_chr_path)
  cli::cli_alert_success(
    paste0(
      "Saved full internal standard chromatograms to: ",
      "{.path {is_chr_path}}"
    )
  )
}

ranges <- get_rt_mz_range(chromatogram = is_chr, rt_window = 0.02)

# ==============================================================================
# Plotting full IS chromatogram ------------------------------------------------
# ==============================================================================
all_is_full_p_path <- file.path(
  snakemake@params$output,
  "graphs",
  "internal_standard",
  "all_is_full.pdf"
)
pdf(all_is_full_p_path)
plot(x = is_chr, col = group_colors[is_chr$group], lwd = 3)
legend("topright", legend = names(group_colors), col = group_colors, pch = 16)
invisible(dev.off())
cli::cli_alert_success(
  paste0(
    "Saved full internal standard chromatogram plot to: ",
    "{.path {all_is_full_p_path}}"
  )
)

# ==============================================================================
# Getting a narrow internal standard XIC ---------------------------------------
# ==============================================================================
is_eic_path <- file.path(snakemake@params$output, "objects", "is_eic.rds")
if (interactive() && file.exists(is_eic_path)) {
  is_eic <- readRDS(file = is_eic_path)
  cli::cli_alert_success(
    paste0(
      "Imported the saved, narrow internal standard XIC from: ",
      "{.path {is_eic_path}}"
    )
  )
} else {
  cli::cli_alert_info("Creating a narrow internal standard XIC")
  is_eic <- xcms::chromatogram(
    BPPARAM = BiocParallel::SerialParam(),
    chunkSize = 1L,
    object = ms_exp,
    mz = ranges$mz_range,
    rt = ranges$rt_range,
    aggregationFun = "sum"
  )
  saveRDS(object = is_eic, file = is_eic_path)
  cli::cli_alert_success(
    paste0(
      "Saved a narrow internal standard XIC to: ",
      "{.path {is_eic_path}}"
    )
  )
}

# ==============================================================================
# Plot all the internal standard XICs together ---------------------------------
# ==============================================================================
all_is_p_path <- file.path(
  snakemake@params$output,
  "graphs",
  "internal_standard",
  "all_is.pdf"
)
pdf(all_is_p_path)
plot(x = is_eic, col = group_colors[is_eic$group], lwd = 3)
legend("topleft", legend = names(group_colors), col = group_colors, pch = 16)
invisible(dev.off())
cli::cli_alert_success(
  paste0(
    "Saved internal standard chromatogram plot to: ",
    "{.path {all_is_p_path}}"
  )
)


# ==============================================================================
# Plot all the internal standard XICs individually -----------------------------
# ==============================================================================
for (i in seq_along(is_eic)) {
  tmp_is_eic_p_path <- file.path(
    snakemake@params$output,
    "graphs",
    "internal_standard",
    paste0(
      colnames(bpcs)[i],
      ".pdf"
    )
  )

  pdf(tmp_is_eic_p_path)
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
cli::cli_alert_success(
  paste0(
    "Saved all individual internal standard chromatogram plots to: ",
    "{.path {dirname(tmp_is_eic_p_path)}}"
  )
)

# ==============================================================================
# Determining minimal and maximal peakwidth based on the internal standard
# - Ensure IS peak is chosen
# - Determine max and min peak width for CentWaveParam() from the IS.
# TODO
# do this for several peaks and not only the IS
# ==============================================================================

cli::cli_alert_info(
  paste0(
    "Determining minimal and maximal peak width based on the internal standard"
  )
)

is_chr_wide_path <- file.path(
  snakemake@params$output,
  "objects",
  "is_chr_wide.rds"
)
if (interactive() && file.exists(is_chr_wide_path)) {
  is_chr_wide <- readRDS(file = is_chr_wide_path)
  cli::cli_alert_success(
    paste0(
      "Imported saved wide internal standard chromatograms from: ",
      "{.path {is_chr_wide_path}}"
    )
  )
} else {
  cli::cli_alert_info("Creating wide internal standard chromatogram")
  is_chr_wide <- xcms::chromatogram(
    BPPARAM = BiocParallel::SerialParam(),
    chunkSize = 1L,
    object = ms_exp,
    mz = mz_range + c(-0.05, 0.05),
    rt = ranges$rt_range + c(-16, 16),
    aggregationFun = "sum"
  )
  saveRDS(object = is_chr_wide, file = is_chr_wide_path)
  cli::cli_alert_success(
    paste0(
      "Saved wide internal standard chromatogram to: ",
      "{.path {is_chr_wide_path}}"
    )
  )
}

# ==============================================================================
# Run peak detection on the IS EIC ---------------------------------------------
# ==============================================================================
is_chr_wide_peaks_path <- file.path(
  snakemake@params$output,
  "objects",
  "is_chr_wide_peaks.rds"
)
if (interactive() && file.exists(is_chr_wide_peaks_path)) {
  is_chr_wide_peaks <- readRDS(file = is_chr_wide_peaks_path)
  cli::cli_alert_success(
    paste0(
      "Imported saved wide internal standard chromatogram peaks from: ",
      "{.path {is_chr_wide_peaks_path}}"
    )
  )
} else {
  cli::cli_alert_info("Calling peaks on wide internal standard chromatograms")
  is_chr_wide_peaks <- xcms::findChromPeaks(
    BPPARAM = BiocParallel::SerialParam(),
    chunkSize = 1,
    object = is_chr_wide,
    param = xcms::CentWaveParam(
      ppm = snakemake@params$ppm_global,
      peakwidth = c(2, 20),
      prefilter = c(1, 1),
      snthresh = snakemake@params$sn_threshold, # 10
      mzCenterFun = "wMean",
      mzdiff = 0.001,
      integrate = 2,
      noise = 1000,
      verboseBetaColumns = TRUE
    ),
    ms_level = 1
  )
  saveRDS(object = is_chr_wide_peaks, file = is_chr_wide_peaks_path)
  cli::cli_alert_success(
    paste0(
      "Saved wide internal standard chromatogram peaks to: ",
      "{.path {is_chr_wide_peaks_path}}"
    )
  )
}

# ==============================================================================
# Calculate peakwidths to use for full call ------------------------------------
# ==============================================================================
is_peaks <- tibble::as_tibble(
  x = xcms::chromPeaks(is_chr_wide_peaks),
  rownames = "rownames"
) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(delta_rt = rtmax - rtmin) %>%
  dplyr::ungroup()


# Min: to half of some peaks in the datasets
# Max: 2-4x times the average size
is_min_peak_width <- min(is_peaks$delta_rt, na.rm = TRUE)
is_max_peak_width <- max(is_peaks$delta_rt, na.rm = TRUE)
min_peak_width <- unname(quantile(is_peaks$delta_rt, 0.05, na.rm = TRUE) * 0.1)
max_peak_width <- unname(quantile(is_peaks$delta_rt, 0.95, na.rm = TRUE) * 8)

cli::cli_ul("Internal standard: {snakemake@params$internal_standard}")
cli::cli_ul("Internal standard adduct: {snakemake@params$is_adduct}")
cli::cli_ul("Theoretical m/z: {round(mz_theory, 2)}")
cli::cli_ul("Min IS peak width: {round(is_min_peak_width, 2)}")
cli::cli_ul("Max IS peak width: {round(is_max_peak_width, 2)}")

saveRDS(
  object = list(
    ranges = ranges,
    min_peak_width = min_peak_width,
    max_peak_width = max_peak_width
  ),
  file = snakemake@output[[1]]
)

script_footer()
end_log()