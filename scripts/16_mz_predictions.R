# ==============================================================================
# m/z predictions subset -------------------------------------------------------
# ==============================================================================
message(
  "Predicting potential biotransformations based on:\n\t",
  "Biotransformation database: ", opt$biotransf_file, # + the other kegg stuff
  "\n\tppm: ", opt$ppm_match,
  sep = ""
)

# Checking specifically for the glycoside anad aglycone m/zs
glycoside <- MetaboCoreUtils::mass2mz(
  MetaboCoreUtils::calculateMass(opt$glycoside)[[1]],
  adduct = MetaboCoreUtils::adducts(polarity = opt$polarity)) %>%
  t() %>%
  tibble::as_tibble(., rownames = "adduct") %>%
  dplyr::rename("glycoside" = V1)

aglycone <- MetaboCoreUtils::mass2mz(
  MetaboCoreUtils::calculateMass(opt$aglycone)[[1]],
  adduct = MetaboCoreUtils::adducts(polarity = opt$polarity)) %>%
  t() %>%
  tibble::as_tibble(., rownames = "adduct") %>%
  dplyr::rename("aglycone" = V1)

# This is okay for now since it's only looking for the glycone and aglycone
# this functions sucks though
range_tol <- ppm_to_num(2000) # used to be glycoside_ppm at 2000

gly_agly_adducts <- glycoside %>%
  dplyr::left_join(
    x = .,
    y = aglycone,
    by = "adduct"
  ) %>%
  dplyr::mutate(
    glycoside_min = glycoside - range_tol,
    glycoside_max = glycoside + range_tol,
    aglycone_min = aglycone - range_tol,
    aglycone_max = aglycone + range_tol
  )

gly_agly <- tibble::tibble()
for (i in seq_along(gly_agly_adducts)) {
  tmp <- full_raw_filled %>%
    dplyr::filter(
      dplyr::between(
        mzmed,
        gly_agly_adducts[i, ]$glycoside_min,
        gly_agly_adducts[i, ]$glycoside_max
      ) |
        dplyr::between(
          mzmed,
          gly_agly_adducts[i, ]$aglycone_min,
          gly_agly_adducts[i, ]$aglycone_max
        )
    ) %>%
    dplyr::mutate(adduct = gly_agly_adducts[i, ]$adduct) %>%
    dplyr::relocate(adduct, .after = "feature")

  gly_agly <- bind_rows(gly_agly, tmp)
}

pot_glycosides <- unique(gly_agly$feature)
subset_matched_diffs <- pred_biot(
  data = possible_adducts_signif,
  biotransf_data = bio_transf2, # bio_transf
  tolerance_ppm = opt$ppm_match,
  features_of_interest = pot_glycosides,
  parallel = TRUE
) %>%
  dplyr::mutate(
    pair = purrr::map2(
      .x = feat1,
      .y = feat2,
      .f = c
    ),
    obs_diff = abs(obs_delta_mass - delta_mass)
  )

glycoside_pairs <- unique(
  c(
    subset_matched_diffs$feat1,
    subset_matched_diffs$feat2
  )
)

# ==============================================================================
# m/z predictions all ----------------------------------------------------------
# ==============================================================================

# Turn this around
if (file.exists(file.path(opt$output, "objects", "matched_diffs.rds"))) {
  matched_diffs <- readRDS(
    file = file.path(opt$output, "objects", "matched_diffs.rds")
  )
} else {
  matched_diffs <- pred_biot(
    data = possible_adducts_signif,
    biotransf_data = bio_transf2,
    tolerance_ppm = opt$ppm_match, # try 5 and 10, # 15 too much
    parallel = TRUE
  )
  saveRDS(
    object = matched_diffs,
    file = paste0(opt$output, "/objects/matched_diffs.rds")
  )
}

matched_diffs2 <- matched_diffs %>%
  dplyr::filter(
    dplyr::if_all(
      .cols = dplyr::all_of(c("mass1", "mass2")),
      .fns = ~ . > 0
    )
  ) %>%
  dplyr::filter(feat1 != feat2)
  # Too slow just do for a few when filtered
  # dplyr::mutate(
  #   pair = purrr::map2(
  #     .x = feat1,
  #     .y = feat2,
  #     .f = c
  #   ),
  #   obs_diff = abs(obs_delta_mass - delta_mass)
  # )

message("Writing predictions to table...")
# TODO This needs filtering first
# readr::write_csv(
#   x = matched_diffs2,
#   file = file.path(opt$output, "tables", "matched_diffs.csv")
# )