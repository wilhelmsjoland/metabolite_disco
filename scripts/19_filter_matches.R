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
    library(dplyr)
    library(tibble)
    library(readr)
    library(purrr)
    library(RSQLite)
    library(MsBackendSql)
    library(httr2)
  })
)

anno_data <- readRDS(snakemake@input[["annotation"]])
biot_data <- readRDS(snakemake@input[["biotransformer"]])

anno <- anno_data$anno
predicted_feats <- biot_data$predicted_feats

# ==============================================================================
# Molecular similarity m/z matching --------------------------------------------
# ==============================================================================
anno_sims_final_path <- file.path(
  snakemake@params$output,
  "tables",
  "annotation_predictions.csv"
)

if (interactive() && file.exists(anno_sims_final_path)) {
  anno_sims_final <- readr::read_csv(
    file = anno_sims_final_path,
    show_col_types = FALSE,
    progress = FALSE
  )
  cli::cli_alert_success(
    paste0(
      "{.path {anno_sims_final_path}} already exists, ",
      "imported database annotated features"
    )
  )
} else {
  anno_sims <- anno %>%
    dplyr::group_by(adduct) %>%
    dplyr::distinct(target_inchikey, .keep_all = TRUE) %>%
    dplyr::ungroup()

  cli::cli_progress_step("Annotating dataframes with info from pubchem")
  if (!is.null(snakemake@params$annotate_path)) {
    pubchem_con <- RSQLite::dbConnect(
      drv = RSQLite::SQLite(),
      snakemake@params$annotate_path
    )

    full_inchikey_map <- dplyr::tbl(pubchem_con, "compounds") %>%
      dplyr::filter(inchikey %in% anno_sims$target_inchikey) %>%
      dplyr::collect()

    RSQLite::dbDisconnect(pubchem_con)
  } else {
    full_inchikey_map <- fetch_pubchem_annotations(
      inchikeys = unique(anno_sims$target_inchikey)
    )
  }
  cli::cli_progress_done()

  anno_sims_final <- anno_sims %>%
    dplyr::left_join(
      x = .,
      y = full_inchikey_map,
      by = c("target_inchikey" = "inchikey"),
      # Expected: same inchikey can repeat across adducts (x), and PubChem
      # can have multiple CIDs for one inchikey (y) - keeping all of them
      # since there's no principled way to pick a single "best" CID.
      relationship = "many-to-many"
    )

  readr::write_csv(
    x = anno_sims_final,
    file = anno_sims_final_path
  )
  cli::cli_alert_success(
    paste0(
      "Saved filtered database annotated features to: ",
      "{.path {anno_sims_final_path}}"
    )
  )
}

# ==============================================================================
# Molecular similarity biotransformer ------------------------------------------
# ==============================================================================
chem_pred_feats_path <- file.path(
  snakemake@params$output,
  "tables",
  "biotransformer_predictions.csv"
)

if (interactive() && file.exists(chem_pred_feats_path)) {
  chem_pred_feats <- readr::read_csv(
    file = chem_pred_feats_path,
    show_col_types = FALSE,
    progress = FALSE
  )
  cli::cli_alert_success(
    paste0(
      "{.path {chem_pred_feats_path}} already exists, ",
      "imported biotransformer predicted features"
    )
  )
} else {
  if (!is.null(snakemake@params$annotate_path)) {
    pubchem_con <- RSQLite::dbConnect(
      drv = RSQLite::SQLite(),
      snakemake@params$annotate_path
    )

    predicted_inchikeys <- unique(predicted_feats$inchikey)
    biotransformer_inchi_map <- dplyr::tbl(pubchem_con, "compounds") %>%
      dplyr::filter(inchikey %in% predicted_inchikeys) %>%
      dplyr::collect()

    RSQLite::dbDisconnect(pubchem_con)
  } else {
    biotransformer_inchi_map <- NULL
  }
  chem_pred_feats <- predicted_feats %>%
    {
      if (!is.null(biotransformer_inchi_map)) {
        dplyr::left_join(
          x = .,
          y = biotransformer_inchi_map,
          by = "inchikey"
        )
      } else {
        .
      }
    } %>%
    # This is needed since one inchikey can be annotated to
    # several features because they can have different adducts
    dplyr::mutate(met_id = as.character(dplyr::row_number())) %>%
    dplyr::relocate(met_id, .before = "inchikey") %>%
    # Here don't need unique since it's not matched against a db
    # so shouldnt have duplicates for the same feature and adduct
    dplyr::group_by(feature, adduct) %>%
    dplyr::distinct(inchikey, .keep_all = TRUE)

  readr::write_csv(
    x = chem_pred_feats,
    file = chem_pred_feats_path
  )

  cli::cli_alert_success(
    paste0(
      "Saved filtered biotransformer predicted features to: ",
      "{.path {chem_pred_feats_path}}"
    )
  )
}

# ==============================================================================
# Snakesave --------------------------------------------------------------------
# ==============================================================================
saveRDS(
  object = list(
    anno_sims_final = anno_sims_final,
    chem_pred_feats = chem_pred_feats
  ),
  file = snakemake@output[[1]]
)

script_footer()
end_log()