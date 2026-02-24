cli::cli_h1(basename(this.path::this.path()))
# ==============================================================================
# Biotransformer.jar -----------------------------------------------------------
# ==============================================================================
cli::cli_h3("Predicting biotransformations from SMILES")
prediction_path <- file.path(opt$output, "tables", "prediction.csv")
if (file.exists(prediction_path)) {
  biot_pred <- readr::read_csv(
    file = prediction_path,
    show_col_types = FALSE,
    progress = FALSE
  )
  cli::cli_alert_success(
    paste0(
      "Imported predictions for {.val {opt$smiles}} from: ",
      "{.path {prediction_path}}"
    )
  )
} else {
  cli::cli_alert_info(
    paste0(
      "Predicting biotransformations from {.val {opt$smiles}} ",
      "with biotransformer.jar"
    )
  )
  run_biotransformer(
    bt_dir = opt$biot_dir,
    smiles = opt$smiles,
    b_type = "superbio",
    k_task = "pred",
    output_file = "prediction"
  )

  cli::cli_alert_success(
    paste0(
      "Saved biotransformer predictions to: ",
      "{.path {prediction_path}}"
    )
  )

  biot_pred <- readr::read_csv(
    file = prediction_path,
    show_col_types = FALSE,
    progress = FALSE
  )
}

cli::cli_alert_success("Merging predicted features with feature definitions")
predicted_feats_path <- file.path(
  opt$output,
  "tables",
  "predicted_annotated_feats.csv"
)
if (file.exists(predicted_feats_path)) {
  predicted_feats <- readr::read_csv(
    file = predicted_feats_path,
    show_col_types = FALSE,
    progress = FALSE
  )
  cli::cli_alert_success(
    paste0(
      "Imported biotransformer predicted features from: ",
      "{.path {predicted_feats_path}}"
    )
  )
} else {
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
    adduct = MetaboCoreUtils::adducts(polarity = opt$polarity)
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
    nrow(MetaboCoreUtils::adducts(polarity = opt$polarity))

  if (biot_mass_len != nrow(biot_final)) {
    cli::cli_alert_danger(
      "The transformation prediction dataframes are not the same length."
    )
  }

  def_tib <- xchr9_filt
  # biot_final = mass to mzs - > match the m/zs to the m/zs in the data
  predicted_feats <- biot_final %>%
    dplyr::inner_join(
      x = .,
      y = def_tib %>%
        dplyr::mutate(
          tol = MsCoreUtils::ppm(mzmed, opt$ppm_match),
          mz_lo = mzmed - tol,
          mz_hi = mzmed + tol
        ) %>%
        dplyr::relocate(c("mz_lo", "mz_hi"), .after = "feature"),
      by = dplyr::join_by(dplyr::between(mz, mz_lo, mz_hi))
    )
  cli::cli_alert_success("Merged predicted features with feature definitions")

  readr::write_csv(
    x = predicted_feats,
    file = predicted_feats_path
  )
  cli::cli_alert_success(
    paste0(
      "Saved biotransformer predicted features to: ",
      "{.path {predicted_feats_path}}"
    )
  )
}

pred_peak_ids <- sort(unique(predicted_feats$feature))

pred_chrs_path <- file.path(
  opt$output,
  "objects",
  "pred_chrs.rds"
)
if (file.exists(pred_chrs_path)) {
  pred_chrs <- readRDS(pred_chrs_path)
  cli::cli_alert_success(
    paste0(
      "Imported biotransformer predicted feature chromatograms from: ",
      "{.path {pred_chrs_path}}"
    )
  )
} else {
  cli::cli_alert_info(
    paste0(
      "Generating chromatograms for annotated ",
      "biotransformer predicted features"
    )
  )
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
    file = pred_chrs_path
  )
  cli::cli_alert_success(
    paste0(
      "Saved biotransformer predicted feature chromatograms: ",
      "{.path {pred_chrs_path}}"
    )
  )
}