# ==============================================================================
# Source functions and minimal startup parameters ------------------------------
# ==============================================================================
source("scripts/functions.R")
start_log(snakemake@params$output)
script_header()

set.seed(snakemake@params$seed)
register_parallel(snakemake@params$cores)
suppressWarnings(
  suppressPackageStartupMessages({
    library(cli)
    library(BiocParallel)
    library(dplyr)
    library(tibble)
    library(tidyr)
    library(purrr)
    library(readr)
    library(MetaboCoreUtils)
    library(MsCoreUtils)
    library(future)
    library(future.apply)
    library(RSQLite)
    library(MsBackendSql)
  })
)

limma_data <- readRDS(snakemake@input[["limma"]])
prep_data <- readRDS(snakemake@input[["prep_biot"]])

intensities <- limma_data$intensities
all_sig_diff <- prep_data$all_sig_diff
possible_adducts_signif <- prep_data$possible_adducts_signif
bio_transf2 <- prep_data$bio_transf2

# ==============================================================================
# Fetching potential glycoside/aglycone features -------------------------------
# ==============================================================================
cli::cli_h3(
  paste0(
    "Generating biotransformation predictions from ",
    "mass-to-charge ratio derived mass estimations"
  )
)

# "C21H20O10;C15H10O6" - any number of formulas to search for, each checked
# independently (OR'd together) rather than as a fixed pair
metabolite_formulas <- unlist(strsplit(snakemake@params$metabolite_search, ";"))

# This is okay for now since it's only looking for the glycone and aglycone
range_tol <- 0.002 # used to be glycoside_ppm at 2000

metabolite_mz <- purrr::map_dfr(
  .x = metabolite_formulas,
  .f = ~ MetaboCoreUtils::mass2mz(
    x = MetaboCoreUtils::calculateMass(.x)[[1]],
    adduct = MetaboCoreUtils::adducts(polarity = snakemake@params$polarity)
  ) %>%
    t() %>%
    as.data.frame() %>%
    tibble::as_tibble(
      x = .,
      rownames = "adduct",
      .name_repair = "universal_quiet"
    ) %>%
    dplyr::rename("mz" = V1) %>%
    dplyr::mutate(formula = .x, .before = "adduct")
) %>%
  dplyr::mutate(
    mz_min = mz - range_tol,
    mz_max = mz + range_tol
  )

metabolite_search <- tibble::tibble()
for (i in seq_len(nrow(metabolite_mz))) {
  tmp <- intensities$raw_fill_imp$untransformed %>%
    dplyr::filter(
      dplyr::between(
        mzmed,
        metabolite_mz[i, ]$mz_min,
        metabolite_mz[i, ]$mz_max
      )
    ) %>%
    dplyr::mutate(
      formula = metabolite_mz[i, ]$formula,
      adduct = metabolite_mz[i, ]$adduct
    ) %>%
    dplyr::relocate(c(formula, adduct), .after = "feature")

  metabolite_search <- bind_rows(metabolite_search, tmp)
}

metabolite_search_path <- file.path(
  snakemake@params$output,
  "tables",
  "metabolite_search.csv"
)

if (interactive() && file.exists(metabolite_search_path)) {
  cli::cli_alert_info("Metabolite search table already exists, skipping")
} else {
  readr::write_csv(metabolite_search, metabolite_search_path)
  cli::cli_alert_info(
    paste0(
      "Metabolite search table saved to {.path {metabolite_search_path}}"
    )
  )
}

pot_metabolite_hits <- unique(metabolite_search$feature)
cli::cli_alert_info(
  paste0(
    "Potential metabolite hits: ",
    "{.val {pot_metabolite_hits}}"
  )
)

# ==============================================================================
# m/z predictions subset -------------------------------------------------------
# ==============================================================================
cli::cli_bullets(
  c(
    "i" = "Predicting subset of potential biotransformations based on: ",
    # below: + the other kegg stuff
    "i" = "\tdatabase: {.val {snakemake@params$mass_shift_path}}",
    "i" = "\tppm: {.val {snakemake@params$ppm_match}}"
  )
)

subset_matched_diffs_path <- file.path(
  snakemake@params$output,
  "objects",
  "subset_matched_diffs.rds"
)

if (interactive() && file.exists(subset_matched_diffs_path)) {
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
    tolerance_ppm = snakemake@params$ppm_match,
    features_of_interest = pot_metabolite_hits,
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
if (snakemake@params$all_vs_all) {
  cli::cli_bullets(
    c(
      "i" = "Predicting all potential biotransformations based on: ",
      # + the other kegg stuff
      "i" = "\tdatabase: {.val {snakemake@params$mass_shift_path}}",
      "i" = "\tppm: {.val {snakemake@params$ppm_match}}"
    )
  )
  # Turn this around
  matched_diffs_path <- file.path(
    snakemake@params$output,
    "objects",
    "matched_diffs.rds"
  )

  if (interactive() && file.exists(matched_diffs_path)) {
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
      tolerance_ppm = snakemake@params$ppm_match,
      parallel = TRUE
    )
    saveRDS(
      object = matched_diffs,
      file = paste0(snakemake@params$output, "/objects/matched_diffs.rds")
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
  #   snakemake@params$output,
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
} else {
  cli::cli_alert_info(
    paste0(
      "all_vs_all is {.val {snakemake@params$all_vs_all}}",
      ", skipping all-vs-all biotransformation predictions"
    )
  )
}

# ==============================================================================
# Snakesave --------------------------------------------------------------------
# ==============================================================================
saveRDS(
  object = list(
    subset_matched_diffs = subset_matched_diffs,
    glycoside_pairs = glycoside_pairs,
    matched_diffs = if (exists("matched_diffs")) matched_diffs else NULL,
    metabolite_search = metabolite_search
  ),
  file = snakemake@output[[1]]
)
end_log()