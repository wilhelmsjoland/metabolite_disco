# ==============================================================================
# Matching m/z's against databases -------------------------------------------
# ==============================================================================
peaks_used <- full_norm_filled %>%
  dplyr::select(feature, mzmed, rtmed) %>%
  dplyr::rename("mz" = "mzmed", "rtime" = "rtmed") %>%
  # All comps with at least one significant difference
  dplyr::filter(feature %in% all_sig_diff) %>%
  tibble::column_to_rownames(var = "feature")

peaks_used$peak_id <- rownames(peaks_used) # keep XCMS peak IDs

annotation_hub <- AnnotationHub::AnnotationHub()
cdb <- annotation_hub[["AH119519"]]

target_df <- ProtGenerics::compounds(
  cdb,
  columns = c(
    "compound_id",
    "name",
    "formula",
    "exactmass",
    "smiles",
    "inchi",
    "inchikey",
    "cas",
    "pubchem"
  )
)

# parameters to match by
mz_match_param <- MetaboAnnotation::Mass2MzParam(
  adducts = c(MetaboCoreUtils::adductNames(polarity = opt$polarity)),
  ppm = opt$ppm_match
)

matches <- MetaboAnnotation::matchValues(
  query = peaks_used,
  target = target_df,
  param = mz_match_param
)

anno <- MetaboAnnotation::matchedData(matches) %>%
  tibble::as_tibble(., rownames = "feature") %>%
  dplyr::mutate(abs_score = abs(score)) %>%
  dplyr::arrange(abs_score)