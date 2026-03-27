################################################################################
# Only features that are interesting mabye????
################################################################################
cp <- xcms::chromPeaks(anno_chrs) %>%
  tibble::as_tibble(rownames = "peak")
fd <- xcms::featureDefinitions(anno_chrs) %>%
  tibble::as_tibble(rownames = "feature")

fd_filtered <- fd %>%
  dplyr::filter(feature %in% shared_features) %>%
  tidyr::unnest(peakidx) %>%
  dplyr::left_join(
    x = .,
    y = cp %>%
      dplyr::mutate(peakidx = row_number()) %>%
      dplyr::select(peakidx, beta_cor),
    by = "peakidx"
  ) %>%
  dplyr::group_by(feature) %>%
  dplyr::summarise(
    median_beta_cor = median(beta_cor, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::filter(median_beta_cor > 0.3)

fd_final <- fd %>%
  dplyr::filter(feature %in% shared_features) %>%
  dplyr::filter(feature %in% fd_filtered$feature)