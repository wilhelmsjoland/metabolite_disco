filt_features_old <- function(
  object = NULL,
  beta_cor_threshold = 0.3,
  beta_snr_threshold = 6,
  sn_threshold = 10,
  filt_vector = NULL
) {
  filt_chrompeaks_tib <- tibble::as_tibble(
    xcms::chromPeaks(object),
    rownames = "feature"
  ) %>%
    dplyr::filter(!is.na(beta_cor) & !is.na(beta_snr)) %>%
    dplyr::group_by(feature) %>%
    dplyr::filter(
      beta_cor >= beta_cor_threshold &
        beta_snr >= beta_snr_threshold
    ) %>%
    dplyr::filter(sn >= sn_threshold) %>%
    dplyr::mutate(feature2 = as.numeric(gsub("[A-Za-z]", "", feature))) %>%
    dplyr::relocate(feature2, .after = feature)

  filt_features_tib <- tibble::as_tibble(
    xcms::featureDefinitions(object),
    rownames = "feature"
  ) %>%
    tidyr::unnest(peakidx) %>%
    dplyr::filter(peakidx %in% filt_chrompeaks_tib$feature2) %>%
    tidyr::nest(data = peakidx)

  filt_sig_features_tib <- filt_features_tib %>%
    # Added so we don't remove interesting ones
    # that don't have a predicted biotransformation
    dplyr::filter(feature %in% filt_vector)

  filt_sig_features <- filt_sig_features_tib$feature

  # Now only look at features that have at least
  # one significantly different feature
  biot_filt_sig_feature_tib <- matched_diffs %>%
    dplyr::filter(
      dplyr::if_any(
        tidyselect::all_of(c("feat1", "feat2")),
        ~ .x %in% filt_sig_features
      )
    ) %>%
    dplyr:::mutate(
      pair = purrr::map2(feat1, feat2, ~ c(.x, .y)),
      # or ppm global
      mz1_forms = purrr::map(
        mz1, ~ Rdisop::getFormula(Rdisop::decomposeMass(.x, ppm = 0))
      ),
      mz2_forms = purrr::map(
        mz2, ~ Rdisop::getFormula(Rdisop::decomposeMass(.x, ppm = 0))
      )
    )

  biot_filt_sig_features <- unique(
    c(
      biot_filt_sig_feature_tib$feat1,
      biot_filt_sig_feature_tib$feat2
    )
  )

  final_plotting_features <- unique(
    c(
      filt_sig_features,
      biot_filt_sig_features
    )
  )

  filt_list <- list(
    # filtering features below:

    # quality filtered peak tib
    "filt_chrompeaks_tib" = filt_chrompeaks_tib,
    # quality filtered feature tib
    "filt_features_tib" = filt_features_tib,
    # quality + sig filtered feature tib
    "filt_sig_features_tib" = filt_sig_features_tib,
    # quality + sig filtered features
    "filt_sig_features" = filt_sig_features,
    # biotransformation features below:

    # quality + sig filtered biotransf tib
    "biot_filt_sig_features_tib" = biot_filt_sig_feature_tib,
    # quality + sig filtered feature
    "biot_filt_sig_features" = biot_filt_sig_features,
    # final features for plotting: below

    # all feats in biot and in filt_sig_features
    "final_plotting_features" = final_plotting_features
  )

  return(filt_list)

}

filt_features <- function(
  object = NULL,
  beta_cor_threshold = 0.8,
  beta_snr_threshold = 3,
  sn_threshold = 10
) {
  filt_chrompeaks_tib <- tibble::as_tibble(
    xcms::chromPeaks(object),
    rownames = "feature"
  ) %>%
    dplyr::filter(!is.na(beta_cor) & !is.na(beta_snr)) %>%
    dplyr::group_by(feature) %>%
    dplyr::filter(
      beta_cor >= beta_cor_threshold &
        beta_snr >= beta_snr_threshold
    ) %>%
    dplyr::filter(sn >= sn_threshold) %>%
    dplyr::mutate(feature2 = as.numeric(gsub("[A-Za-z]", "", feature))) %>%
    dplyr::relocate(feature2, .after = feature)

  filt_features_tib <- tibble::as_tibble(
    xcms::featureDefinitions(object),
    rownames = "feature"
  ) %>%
    tidyr::unnest(peakidx) %>%
    dplyr::filter(peakidx %in% filt_chrompeaks_tib$feature2) %>%
    tidyr::nest(data = peakidx)

  return(filt_features_tib)

}