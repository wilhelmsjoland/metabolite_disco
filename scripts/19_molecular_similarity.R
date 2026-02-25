cli::cli_h1(basename(this.path::this.path()))
# ==============================================================================
# Molecular similarity m/z matching --------------------------------------------
# ==============================================================================

anno_sims_final_path <- file.path(
  opt$output,
  "tables",
  "anno_similarities.csv"
)
if (file.exists(anno_sims_final_path)) {
  anno_sims_final <- readr::read_csv(
    file = anno_sims_final_path,
    show_col_types = FALSE,
    progress = FALSE
  )
  cli::cli_alert_success(
    paste0(
      "{.path {anno_sims_final_path}} already exists, ",
      "imported molecular similarities of database annotated features"
    )
  )
} else {
  cli::cli_progress_step(
    paste0(
      "Calculate molecular similarity of databased annotated features to: ",
      "{.val {opt$smiles}} using similarity cutoff 0.1"
    )
  )
  anno_filt <- anno %>%
    dplyr::group_by(adduct) %>%
    dplyr::distinct(target_inchikey, .keep_all = TRUE)

  anno_smiles <- anno_filt$target_smiles
  names(anno_smiles) <- anno_filt$feature

  anno_sims <- mol_similarity(
    query_smiles = opt$smiles,
    target_smiles = anno_smiles,
    kekulise = TRUE, # parsing incorrect smiles with electrons
    omit_nulls = TRUE,
    fingerprint = "circular",
    circular_type = "ECFP6",
    method = "tanimoto"
  ) %>%
    dplyr::left_join(
      x = .,
      y = anno,
      by = "feature"
    ) %>%
    # Arbitrary for now
    dplyr::filter(sim > 0.1)

  cli::cli_progress_done()

  # Get names for individual inchikey
  inchi_ks_split <- split(
    anno_sims$target_inchikey,
    ceiling(seq_along(anno_sims$target_inchikey) / 500)
  )
  inchi_ks <- lapply(
    X = inchi_ks_split,
    FUN = function(x) {
      paste0(x, collapse = ",")
    }
  )

  cli::cli_progress_step(
    paste0(
      "Annotating molecular similarity dataframes with info from pubchem"
    )
  )
  full_inchikey_map <- tibble::tibble()
  for (i in seq_along(inchi_ks)) {
    cmd <- paste0(
      "curl -s -d \"inchikey=", inchi_ks[[i]], "\" ",
      "\"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/property/",
      "Title,MolecularFormula,InChIKey,ExactMass,XLogP,TPSA/",
      "CSV\""
    )
    results <- system(cmd, intern = TRUE)
    inchikey_map <- readr::read_csv(
      file = I(results),
      show_col_types = FALSE,
      progress = FALSE
    )

    full_inchikey_map <- dplyr::bind_rows(full_inchikey_map, inchikey_map)
  }
  cli::cli_progress_done()

  anno_sims_final <- anno_sims %>%
    dplyr::left_join(
      x = .,
      y = full_inchikey_map,
      by = c("target_inchikey" = "InChIKey")
    )

  readr::write_csv(
    x = anno_sims_final,
    file = anno_sims_final_path
  )
  cli::cli_alert_success(
    paste0(
      "Saved molecular similarity of databased annotated features to: ",
      "{.path {anno_sims_final_path}}"
    )
  )
}

# ==============================================================================
# Molecular similarity biotransformer ------------------------------------------
# ==============================================================================
chem_pred_feats_path <- file.path(
  opt$output,
  "tables",
  "biotransformer_similarities.csv"
)

if (file.exists(chem_pred_feats_path)) {
  chem_pred_feats <- readr::read_csv(
    file = chem_pred_feats_path,
    show_col_types = FALSE,
    progress = FALSE
  )
  cli::cli_alert_success(
    paste0(
      "{.path {chem_pred_feats_path}} already exists, ",
      "imported molecular similarities of biotransformer predicted features"
    )
  )
} else {
  chem_pred_feats <- predicted_feats %>%
    # This is needed since one inchikey can be annotated to
    # several features because they can have different adducts
    dplyr::mutate(met_id = as.character(dplyr::row_number())) %>%
    dplyr::relocate(met_id, .before = "InChIKey") %>%
    # Here don't need unique since it's not matched against a db
    # so shouldnt have duplicates for the same feature and adduct
    dplyr::group_by(feature, adduct) %>%
    dplyr::distinct(InChIKey, .keep_all = TRUE)

  pred_smiles <- chem_pred_feats$SMILES
  names(pred_smiles) <- chem_pred_feats$met_id

  cli::cli_progress_step(
    paste0(
      "Calculate molecular similarity of biotransformer predicted ",
      "features to: {.val {opt$smiles}}"
    )
  )
  pred_sims <- mol_similarity(
    query_smiles = opt$smiles,
    target_smiles = pred_smiles,
    kekulise = TRUE,
    omit_nulls = TRUE,
    fingerprint = "circular",
    circular_type = "ECFP6",
    method = "tanimoto"
  ) %>%
    # rename to not overlap in the dplyr::left_join()
    dplyr::rename("met_id" = "feature") %>%
    dplyr::left_join(
      x = .,
      y = chem_pred_feats,
      by = c("met_id")
    )

  cli::cli_progress_done()

  readr::write_csv(
    x = pred_sims,
    file = chem_pred_feats_path
  )

  cli::cli_alert_success(
    paste0(
      "Saved molecular similarity of biotransformer predicted features to: ",
      "{.path {chem_pred_feats_path}}"
    )
  )
}
