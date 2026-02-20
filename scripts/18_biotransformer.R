# ==============================================================================
# Biotransformer.jar -----------------------------------------------------------
# ==============================================================================
if (!file.exists(file.path(getwd(), opt$output, "tables", "prediction.csv"))) {
  run_biotransformer(
    bt_dir = opt$biot_dir,
    smiles = opt$smiles,
    b_type = "superbio",
    k_task = "pred",
    output_file = "prediction"
  )
} else {
  biot_pred <- readr::read_csv(
    file = file.path(opt$output, "tables", "prediction.csv"),
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
        tol = MsCoreUtils::ppm(mzmed, opt$ppm_match),
        mz_lo = mzmed - tol,
        mz_hi = mzmed + tol
      ) %>%
      dplyr::relocate(c("mz_lo", "mz_hi"), .after = "feature"),
    by = dplyr::join_by(dplyr::between(mz, mz_lo, mz_hi))
  )

pred_peak_ids <- sort(unique(predicted_feats$feature))

if (file.exists(file.path(opt$output, "objects", "pred_chrs.rds"))) {
  pred_chrs <- readRDS(file.path(opt$output, "objects", "pred_chrs.rds"))
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
    file = file.path(opt$output, "objects", "pred_chrs.rds")
  )
}