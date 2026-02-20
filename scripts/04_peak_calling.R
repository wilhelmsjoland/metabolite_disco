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

message(
  "===========================================================================",
  "\nPeak widths used for peak picking: ",
  "\n\t Minimal width: ", round(min_peak_width, 3),
  "\n\t Maximal width: ", round(max_peak_width, 3), "\n",
  "===========================================================================",
  sep = ""
)

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
      snthresh = opt$sn_threshold, # 10
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