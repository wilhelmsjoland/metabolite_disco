cli::cli_h1(basename(this.path::this.path()))
# ==============================================================================
# Fetching potential glycoside/aglycone features -------------------------------
# ==============================================================================
cli::cli_h3(
  paste0(
    "Generating biotransformation predictions from ",
    "mass-to-charge ratio derived mass estimations"
  )
)

# Checking specifically for the glycoside anad aglycone m/zs
glycoside <- MetaboCoreUtils::mass2mz(
  x = MetaboCoreUtils::calculateMass(opt$glycoside)[[1]],
  adduct = MetaboCoreUtils::adducts(polarity = opt$polarity)
) %>%
  t() %>%
  as.data.frame() %>%
  tibble::as_tibble(
    x = .,
    rownames = "adduct",
    .name_repair = "universal_quiet"
  ) %>%
  dplyr::rename("glycoside" = V1)

aglycone <- MetaboCoreUtils::mass2mz(
  x = MetaboCoreUtils::calculateMass(opt$aglycone)[[1]],
  adduct = MetaboCoreUtils::adducts(polarity = opt$polarity)
) %>%
  t() %>%
  as.data.frame() %>%
  tibble::as_tibble(
    x = .,
    rownames = "adduct",
    .name_repair = "universal_quiet"
  ) %>%
  dplyr::rename("aglycone" = V1)

# This is okay for now since it's only looking for the glycone and aglycone
range_tol <- 0.002 # used to be glycoside_ppm at 2000

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
  tmp <- intensities$raw_fill_imp$untransformed %>%
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
cli::cli_alert_info(
  paste0(
    "Potential glycosides/aglycones",
    "{.val {pot_glycosides}}"
  )
)

# ==============================================================================
# m/z predictions subset -------------------------------------------------------
# ==============================================================================
cli::cli_bullets(
  c(
    "i" = "Predicting subset of potential biotransformations based on: ",
    "i" = "\tdatabase: {.val {opt$biotransf_file}}", # + the other kegg stuff
    "i" = "\tppm: {.val {opt$ppm_match}}"
  )
)

subset_matched_diffs_path <- file.path(
  opt$output,
  "objects",
  "subset_matched_diffs.rds"
)

if (file.exists(subset_matched_diffs_path)) {
  subset_matched_diffs <- readRDS(subset_matched_diffs_path)
  cli::cli_alert_success(
    paste0(
      "Imported subsetted m/z predictions object from: ",
      "{.path {subset_matched_diffs_path}}"
    )
  )
} else {
  cli::cli_progress_step("Predicting subsetted m/zs")
  subset_matched_diffs <- pred_biot(
    data = possible_adducts_signif,
    biotransf_data = bio_transf2,
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
      obs_diff = abs(obs_delta_mass - delta_mass),
      ppm_diff = num_to_ppm(mz = delta_mass, diff = obs_diff),
      ppm_mz = num_to_ppm(mz = pmax(mz1, mz2), diff = obs_diff)
    )

  cli::cli_progress_done()
  saveRDS(subset_matched_diffs, subset_matched_diffs_path)

  cli::cli_alert_success(
    paste0(
      "Saved subsetted m/z predictions to: ",
      "{.path {subset_matched_diffs_path}}"
    )
  )
}

glycoside_pairs <- unique(
  c(
    subset_matched_diffs$feat1,
    subset_matched_diffs$feat2
  )
)

# ==============================================================================
# m/z predictions all ----------------------------------------------------------
# ==============================================================================

cli::cli_bullets(
  c(
    "i" = "Predicting all potential biotransformations based on: ",
    "i" = "\tdatabase: {.val {opt$biotransf_file}}", # + the other kegg stuff
    "i" = "\tppm: {.val {opt$ppm_match}}"
  )
)
# Turn this around
matched_diffs_path <- file.path(
  opt$output,
  "objects",
  "matched_diffs.rds"
)

if (file.exists(matched_diffs_path)) {
  matched_diffs <- readRDS(file = matched_diffs_path)
  cli::cli_alert_success(
    paste0(
      "Imported m/z predictions object from: ",
      "{.path {matched_diffs_path}}"
    )
  )
} else {
  cli::cli_alert_info("Predicting m/zs")
  matched_diffs <- pred_biot(
    data = possible_adducts_signif,
    biotransf_data = bio_transf2,
    tolerance_ppm = opt$ppm_match,
    parallel = TRUE
  )
  saveRDS(
    object = matched_diffs,
    file = paste0(opt$output, "/objects/matched_diffs.rds")
  )
  cli::cli_alert_success(
    paste0(
      "Saved subsetted m/z predictions to: ",
      "{.path {matched_diffs_path}}"
    )
  )
}

matched_diffs <- matched_diffs %>%
  dplyr::filter(
    dplyr::if_all(
      .cols = dplyr::all_of(c("mass1", "mass2")),
      .fns = ~ . > 0
    )
  ) %>%
  dplyr::filter(feat1 != feat2)

cli::cli_alert_success(
  paste0(
    "Filtered full m/z predictions"
  )
)

# TODO This needs filtering first
# cli::cli_alert_info("Writing predictions to table")
# match_diffs_table_path <- file.path(
#   opt$output,
#   "tables",
#   "matched_diffs.csv"
# )
# readr::write_csv(x = matched_diffs, file = match_diffs_table_path)
# cli::cli_alert_success(
#   paste0(
#     "Wrote predictions to: ",
#     "{.val {match_diffs_table_path}}"
#   )
# )
