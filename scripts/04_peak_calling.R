# ==============================================================================
# Source functions and minimal startup parameters ------------------------------
# ==============================================================================
source("scripts/functions.R")
start_log(snakemake@params$output)
script_header()

set.seed(snakemake@params$seed)
register_parallel(snakemake@params$cores)
suppressWarnings(
  suppressPackageStartupMessages({
    library(cli)
    library(BiocParallel)
    library(xcms)
    library(tibble)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
  })
)

setup <- readRDS(snakemake@input[["setup"]])
is_data <- readRDS(snakemake@input[["internal_std"]])
ms_exp <- setup$ms_exp
group_colors <- setup$group_colors
meta <- setup$meta
min_peak_width <- is_data$min_peak_width
max_peak_width <- is_data$max_peak_width

# ==============================================================================
# Call peaks on whole dataset with parameters
# ==============================================================================
cli::cli_h3("Calling peaks: ")

xchr_path <- file.path(snakemake@params$output, "objects", "xchr.rds")
if (interactive() && file.exists(xchr_path)) {
  xchr <- readRDS(file = xchr_path)
  cli::cli_alert_success(
    paste0(
      "Imported saved peak calling object from: ",
      "{.path {xchr_path}}"
    )
  )
} else {
  cli::cli_alert_info("Calling peaks")
  xchr <- xcms::findChromPeaks(
    object = ms_exp,
    BPPARAM = BiocParallel::bpparam(),
    chunkSize = snakemake@params$cores,
    mslevel = 1L,
    param = xcms::CentWaveParam(
      ppm = snakemake@params$ppm_global,
      peakwidth = c(min_peak_width, max_peak_width),
      snthresh = snakemake@params$sn_threshold,
      prefilter = c(3, 1000), # k pks (left) over intens (right) # c(4, 1000)
      mzCenterFun = "wMean",
      integrate = 2,
      mzdiff = snakemake@params$mzdiff, # 0.001
      fitgauss = FALSE, # TRUE
      noise = 1000,
      verboseColumns = FALSE,
      roiList = list(),
      firstBaselineCheck = TRUE,
      roiScales = numeric(),
      extendLengthMSW = FALSE,
      verboseBetaColumns = TRUE
    )
  )
  saveRDS(object = xchr, file = xchr_path)
  cli::cli_alert_success(
    paste0(
      "Saved peak calling object to: ",
      "{.path {xchr_path}}"
    )
  )
}

# ==============================================================================
# Print peak calling params used to console ------------------------------------
# ==============================================================================
xchr_params <- xchr@processHistory[[1]]@param
cli::cli_alert_success("Called peaks with:")
param_msg(xchr_params)

# ==============================================================================
# Inspect peaks ----------------------------------------------------------------
# ==============================================================================
cli::cli_h3("Inspecting peaks with the largest area for each sample")
peaks_to_inspect <- tibble::as_tibble(chromPeaks(xchr), rownames = "peak") %>%
  tidyr::drop_na(beta_cor, beta_snr) %>%
  dplyr::filter(beta_cor >= snakemake@params$beta_cor_threshold) %>%
  dplyr::filter(beta_snr >= snakemake@params$beta_snr_threshold) %>%
  dplyr::arrange(desc(into)) %>%
  # Keeps the top peak for every sample
  dplyr::distinct(sample, .keep_all = TRUE) %>%
  dplyr::pull(peak)

inspect_peak_p_path <- file.path(
  snakemake@params$output,
  "graphs",
  "quality_control"
)
for (peak in peaks_to_inspect) {
  if (!paste0(peak, ".pdf") %in% list.files(inspect_peak_p_path)) {
    inspect_peaks <- inspect_peak(
      chrom_obj = xchr,
      peak = peak,
      save_loc = inspect_peak_p_path
    )
    cli::cli_alert_success(
      "Saved peak plot: {.val {peak}} to {inspect_peak_p_path}"
    )
  } else {
    cli::cli_alert_info(
      paste0(
        "{.val {peak}} already saved to",
        " {.path {paste0(inspect_peak_p_path, '/', {peak}, '.pdf')}}"
      )
    )
  }
}

# ==============================================================================
# Plotting per-sample peak counts ----------------------------------------------
# ==============================================================================
# These are the per-sample peak counts
exp_inten_p_path <- file.path(
  snakemake@params$output,
  "graphs",
  "quality_control"
)
exp_inten_p <- plot_experiment_intensities(
  chrom_obj = xchr,
  value = into,
  save_loc = exp_inten_p_path
)
cli::cli_alert_success(
  "Saved per-sample peak counts to {.path {exp_inten_p_path}}"
)

saveRDS(
  object = list(xchr = xchr),
  file = snakemake@output[[1]]
)
end_log()