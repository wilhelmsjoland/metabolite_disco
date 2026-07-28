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
    library(readr)
    library(purrr)
    library(RSQLite)
    library(MsBackendSql)
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
  "anno_similarities.csv"
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
      "imported molecular similarities of database annotated features"
    )
  )
} else {
  anno_sims <- anno %>%
    dplyr::group_by(adduct) %>%
    dplyr::distinct(target_inchikey, .keep_all = TRUE) %>%
    dplyr::ungroup()

  # Get names for individual inchikey
  unique_inchikeys <- unique(anno_sims$target_inchikey) # check each key once
  inchi_ks_split <- split(
    unique_inchikeys,
    ceiling(seq_along(unique_inchikeys) / 500)
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
      by = c("target_inchikey" = "InChIKey"),
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
      "Saved molecular similarity of databased annotated features to: ",
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
  "biotransformer_similarities.csv"
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

  readr::write_csv(
    x = chem_pred_feats,
    file = chem_pred_feats_path
  )

  cli::cli_alert_success(
    paste0(
      "Saved molecular similarity of biotransformer predicted features to: ",
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
end_log()