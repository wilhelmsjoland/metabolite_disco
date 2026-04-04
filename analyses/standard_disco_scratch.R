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

final_peaks <- std_peak_matches %>%
  dplyr::group_by(group, adduct) %>%
  dplyr::slice_max(into, n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(group, dplyr::desc(into))

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
