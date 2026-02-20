# ==============================================================================
# Molecular similarity m/z matching --------------------------------------------
# ==============================================================================
anno_filt <- anno %>%
  # TODO
  # DO A more intelligent filtering than this based on
  # how much information is available in all the rows
  dplyr::group_by(adduct) %>% # Think this might have fixed the loss of same
  dplyr::distinct(target_inchikey, .keep_all = TRUE)

# TODO
# ==============================================================================
# # For checking
# # These are all the base feature without "." added
# anno %>%
#   dplyr::filter(is.na(adduct)) %>%
#   print(n = 300)

# # looks good
# anno %>%
#   dplyr::group_by(peak_id, adduct) %>%
#   dplyr::distinct(target_inchikey, .keep_all = TRUE) %>%
#   dplyr::filter(!is.na(adduct)) %>%
#   dplyr::group_by(feature, adduct, target_inchikey) %>%
#   dplyr::summarize(n = n()) %>%
#   dplyr::filter(n != 1)
# ==============================================================================

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
  # TODO
  # Arbitrary for now
  dplyr::filter(sim > 0.1)

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

full_inchikey_map <- tibble::tibble()
for (i in seq_along(inchi_ks)) {
  cmd <- paste0(
    "curl -s -d \"inchikey=", inchi_ks[[i]], "\" ",
    "\"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/property/",
    "Title,MolecularFormula,InChIKey,ExactMass,XLogP,TPSA/",
    "CSV\""
  )
  results <- system(cmd, intern = TRUE)
  inchikey_map <- readr::read_csv(file = I(results))
  full_inchikey_map <- dplyr::bind_rows(full_inchikey_map, inchikey_map)
}

anno_sims2 <- anno_sims %>%
  dplyr::left_join(
    x = .,
    y = full_inchikey_map,
    by = c("target_inchikey" = "InChIKey")
  )

# ==============================================================================
# Molecular similarity biotransformer ------------------------------------------
# ==============================================================================
chem_pred_feats <- predicted_feats %>%
  # This is needed since one inchikey can be annotated to
  # several features because they can have different adducts
  # ---------------------------------------------- #
  dplyr::mutate(met_id = as.character(dplyr::row_number())) %>%
  dplyr::relocate(met_id, .before = "InChIKey") %>%
  # ---------------------------------------------- #

  # Here don't need unique since it's not matched against a db
  # so shouldnt have duplicates for the same feature and adduct
  # ---------------------------------------------- #
  dplyr::group_by(feature, adduct) %>%
  dplyr::distinct(InChIKey, .keep_all = TRUE)
  # ---------------------------------------------- #

pred_smiles <- chem_pred_feats$SMILES
names(pred_smiles) <- chem_pred_feats$met_id

pred_sims <- mol_similarity(
  query_smiles = opt$smiles,
  target_smiles = pred_smiles,
  kekulise = TRUE, # allow for parsing incorrect smiles with electrons
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