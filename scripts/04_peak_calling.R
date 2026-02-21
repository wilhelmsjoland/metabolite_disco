cli::cli_h1(basename(this.path::this.path()))

# ==============================================================================
# Call peaks on whole dataset with parameters
# ==============================================================================
cli::cli_h3("Calling peaks: ")

if (check_saved("xchr.rds")) {
  xchr <- readRDS(file = paste0(opt$output, "/objects/xchr.rds"))
} else {
  xchr <- xcms::findChromPeaks(
    object = ms_exp,
    BPPARAM = BiocParallel::bpparam(),
    return.type = "XCMSnExp",
    ms_level = 1L,
    param = xcms::CentWaveParam(
      ppm = opt$ppm_global,
      peakwidth = c(min_peak_width, max_peak_width),
      snthresh = opt$sn_threshold,
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
  saveRDS(object = xchr, file = paste0(opt$output, "/objects/xchr.rds"))
}

# ==============================================================================
# Print peak calling params used to console ------------------------------------
# ==============================================================================
xchr_params <- xchr@processHistory[[1]]@param

# instead of peak_calling_param_msg()
# peak_calling_param_msg()
# but kind of confusing I'll admit
purrr::walk(
  .x = slotNames(xchr_params),
  .f = ~ cli::cli_bullets(
    c(
      "i" = paste0(
        .x, ": ",
        paste( # needed to not repeat the names twice for > 1 vectors
          slot(xchr_params, .x),
          collapse = ", "
        )
      )
    )
  )
)

# ==============================================================================
# Inspect peaks ----------------------------------------------------------------
# ==============================================================================
cli::cli_h3("Inspecting peaks with the largest area for each sample")
peaks_to_inspect <- tibble::as_tibble(chromPeaks(xchr), rownames = "peak") %>%
  tidyr::drop_na(beta_cor, beta_snr) %>%
  dplyr::filter(beta_cor >= opt$beta_cor_threshold) %>%
  dplyr::filter(beta_cor >= opt$beta_cor_threshold) %>%
  dplyr::arrange(desc(into)) %>%
  # Keeps the top peak for every sample
  dplyr::distinct(sample, .keep_all = TRUE) %>%
  dplyr::pull(peak)

inspect_peak_p_path <- file.path(opt$output, "graphs", "quality_control")
for (peak in peaks_to_inspect) {
  if (!paste0(peak, ".pdf") %in% list.files(inspect_peak_p_path)) {
    inspect_peaks <- inspect_peak(
      chrom_obj = xchr,
      peak = peak,
      save_loc = file.path(opt$output, "graphs", "quality_control")
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
exp_inten_p_path <- file.path(opt$output, "graphs", "quality_control")
exp_inten_p <- plot_experiment_intensities(
  chrom_obj = xchr,
  value = into,
  save_loc = exp_inten_p_path
)
cli::cli_alert_success(
  "Saved per-sample peak counts to {.path {exp_inten_p_path}}"
)