if (!file.exists(file.path(getwd(), res_folder, "tables", "apigenin.csv"))) {
  run_biotransformer(
    bt_dir = "biotransformer3.0jar",
    smiles = "C1=CC(=CC=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)O",
    b_type = "superbio",
    k_task = "pred",
    output_file = "apigenin"
  )
}

biot_pred <- read_csv(
  file = file.path(res_folder, "tables", "apigenin.csv"),
  show_col_types = FALSE
)

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
  dplyr::mutate(mass = MetaboCoreUtils::calculateMass(`Molecular formula`)) %>%
  dplyr::relocate(mass, .before = "InChI")

biot_mets <- biot_mass$mass
names(biot_mets) <- biot_mass$InChIKey

biot_final <- MetaboCoreUtils::mass2mz(
  x = biot_mets,
  adduct = MetaboCoreUtils::adducts(polarity = "negative")
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
  nrow(adducts(polarity = "negative"))

if (biot_mass_len != nrow(biot_final)) {
  warning("The transformation prediction dataframes are not the same length.")
} else {
  message("The transformation prediction dataframes are the same length.")
}

biot_final

def_tib <- xchr9.filt$filt.features.tib

# biot_final = mass to mzs - > match the m/zs to the m/zs in the data
predicted_feats <- biot_final %>%
  dplyr::inner_join(
    x = .,
    y = def_tib %>% # xchr9.defs
      dplyr::mutate(
        tol = MsCoreUtils::ppm(mzmed, 15),
        mz_lo = mzmed - tol,
        mz_hi = mzmed + tol
      ) %>%
      dplyr::relocate(c("mz_lo", "mz_hi"), .after = "feature"),
    by = dplyr::join_by(dplyr::between(mz, mz_lo, mz_hi))
  )

pred_peak_ids <- unique(predicted_feats$feature)

if (!exists("pred_chrs")) {
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
}

for (i in pred_peak_ids) {
  tmp_anno_chr <- plotFeatChrInt(
    feature_chrom = pred_chrs,
    feature = i,
    method = "sum",
    value = "into",
    filled = TRUE,
    missing = 0, # "rowmin_half"
    msLevel = 1,
    save_loc = NULL,
    device = NULL,
    feat_pairs = FALSE
  )
  print(tmp_anno_chr$combined)
  stop_loop <- readline("Enter for next, break for end: ")
  if (stop_loop == "break") {
    break
  }
}

# FT06025
# FT08376

# REALLY INTERESTING!
to.look.at <- c(
  "FT03581", # 1
  "FT02413", # 2
  "FT08376",
  "FT09870",
  "FT06149",
  "FT04835"
)

# FIX the filtering first
test <- predicted_feats %>%
  dplyr::filter(feature %in% to.look.at) %>%
  dplyr::arrange(InChIKey) %>%
  dplyr::mutate(name = paste(feature, InChIKey)) %>%
  # only for now
  dplyr::distinct(InChIKey, .keep_all = TRUE) %>%
  # only for now
  dplyr::filter(grepl("FT03581", feature))

for (i in test$SMILES) {
  mol.2d <- parse.smiles(i)[[1]]
  img <- view.image.2d(mol.2d)
  img.grob <- grid::rasterGrob(img, interpolate = TRUE)
  img.p <- wrap_elements(img.grob, clip = TRUE)
  print(img.p) 
    stop_loop <- readline("Enter for next, break for end: ")
    if (stop_loop == "break") {
      break
    }
}

biot.sims <- test$SMILES
names(biot.sims) <- test$name
apigenin.smiles <- "C1=CC(=CC=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)O"

biot.sim.final <- mol_similarity(
  query_smiles = apigenin.smiles,
  target_smiles = biot.sims,
  kekulise = TRUE,
  omit_nulls = TRUE,
  fingerprint = "circular",
  circular_type = "ECFP6",
  method = "tanimoto"
) 

biot.sim.final %>%
  dplyr::left_join(
    x = .,
    y = test,
    by = c("feature" = "name")
  ) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = feature, # target_name
      y = sim
    )
  ) +
  # because there are duplicates - just choose the best one
  ggplot2::geom_col(
      stat = "summary",
      fun = "max",
      color = "black",
      position = ggplot2::position_dodge()
  ) +
  ggplot2::guides(x = ggplot2::guide_axis(angle = -45)) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(c(0, 0)),
    limits = c(0, 1)
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_text(angle = -90),
    legend.title = ggplot2::element_blank()
  ) +
  ggplot2::labs(y = "Tanimoto similarity")
  
