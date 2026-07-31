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
    library(dplyr)
    library(readr)
    library(MetaboCoreUtils)
    library(MsCoreUtils)
    library(xcms)
    library(RSQLite) # for sql backend
    library(MsBackendSql) # for sql backend
  })
)

filter_data <- readRDS(snakemake@input[["filter_features"]])
xchr9 <- filter_data$xchr9

# ==============================================================================
# Biotransformer.jar -----------------------------------------------------------
# ==============================================================================
cli::cli_h3("Predicting biotransformations from SMILES")

# "name,smiles;name,smiles" - ";" and "," never appear in plain SMILES, unlike
# "=" (double bonds) or ":" (aromatic bonds/atom maps), and neither collides
# with this project's underscore-heavy naming convention.
smiles_pairs <- unlist(strsplit(snakemake@params$smiles, ";"))
smiles_vec <- setNames(
  sub("^[^,]+,", "", smiles_pairs),
  sub(",.*$", "", smiles_pairs)
)

output_dir <- snakemake@params$output

prediction_path <- function(mol_name, output_dir) {
  file.path(output_dir, "tables", paste0("prediction_", mol_name, ".csv"))
}

cli::cli_alert_info(
  paste0(
    "Predicting biotransformations for {.val {names(smiles_vec)}} ",
    "with biotransformer.jar (in parallel, one process per molecule)"
  )
)
# One process per molecule -- more workers than SMILES strings just spawns
# processes with nothing to do
BiocParallel::bpworkers(bp) <- min(
  BiocParallel::bpworkers(bp),
  length(smiles_vec)
)
biot_pred <- BiocParallel::bplapply(
  X = names(smiles_vec),
  FUN = function(mol_name) {
    run_biotransformer(
      bt_dir = "biotransformer3.0jar",
      smiles = smiles_vec[[mol_name]],
      b_type = "superbio",
      k_task = "pred",
      output_file = paste0("prediction_", mol_name),
      results_path = output_dir,
      annotate_path = snakemake@params$annotate_path
    )

    readr::read_csv(
      file = prediction_path(mol_name, output_dir),
      show_col_types = FALSE,
      progress = FALSE
    ) %>%
      dplyr::mutate(input = mol_name) %>%
      dplyr::rename("inchikey" = "InChIKey")
  },
  BPPARAM = bp
) %>%
  dplyr::bind_rows()

cli::cli_alert_success(
  paste0(
    "Saved biotransformer predictions for {.val {names(smiles_vec)}} to: ",
    "{.path {file.path(output_dir, 'tables')}}"
  )
)

cli::cli_alert_success("Merging predicted features with feature definitions")
predicted_feats_path <- file.path(
  snakemake@params$output,
  "tables",
  "predicted_annotated_feats.csv"
)

biot_dedup <- biot_pred %>%
  dplyr::group_by(inchikey, input) %>%
  dplyr::summarize(
    dplyr::across(
      .cols = setdiff(colnames(.), c("inchikey", "input")),
      .fns  = ~ paste(unique(.x), collapse = ", ")
    ),
    .groups = "keep"
  )

biot_mass <- biot_dedup %>%
  dplyr::ungroup() %>% # changed so row_number is consistent - check!!!!!
  dplyr::mutate(
    mass = MetaboCoreUtils::calculateMass(`Molecular formula`),
    row_id = dplyr::row_number()
  ) %>%
  dplyr::relocate(mass, .before = "InChI")

# row_id (plain, unique integers) is the join key back onto biot_mass below -
# unambiguous, unlike trying to key on InChIKey/input directly (InChIKey
# alone isn't unique anymore now that the same metabolite can come from
# multiple inputs).
biot_mets <- biot_mass$mass
names(biot_mets) <- biot_mass$row_id

biot_final <- MetaboCoreUtils::mass2mz(
  x = biot_mets,
  adduct = MetaboCoreUtils::adducts(polarity = snakemake@params$polarity)
) %>%
  tibble::as_tibble(., rownames = "row_id") %>%
  dplyr::mutate(row_id = as.integer(row_id)) %>%
  tidyr::pivot_longer(
    cols = -row_id,
    names_to = "adduct",
    values_to = "mz"
  ) %>%
  dplyr::left_join(
    x = .,
    y = biot_mass,
    by = "row_id",
    relationship = "many-to-one"
  ) %>%
  dplyr::select(-row_id) %>%
  dplyr::arrange(inchikey)

biot_mass_len <- length(biot_mass$inchikey) *
  nrow(MetaboCoreUtils::adducts(polarity = snakemake@params$polarity))

if (biot_mass_len != nrow(biot_final)) {
  cli::cli_alert_danger(
    "The transformation prediction dataframes are not the same length."
  )
}

def_tib <- xcms::featureDefinitions(xchr9) %>%
  tibble::as_tibble(., rownames = "feature")
# biot_final = mass to mzs - > match the m/zs to the m/zs in the data
predicted_feats <- biot_final %>%
  dplyr::inner_join(
    x = .,
    y = def_tib %>%
      dplyr::mutate(
        tol = MsCoreUtils::ppm(mzmed, snakemake@params$ppm_match),
        mz_lo = mzmed - tol,
        mz_hi = mzmed + tol
      ) %>%
      dplyr::relocate(c("mz_lo", "mz_hi"), .after = "feature"),
    by = dplyr::join_by(dplyr::between(mz, mz_lo, mz_hi))
  )
cli::cli_alert_success("Merged predicted features with feature definitions")

readr::write_csv(
  x = predicted_feats,
  file = predicted_feats_path
)
cli::cli_alert_success(
  paste0(
    "Saved biotransformer predicted features to: ",
    "{.path {predicted_feats_path}}"
  )
)

pred_peak_ids <- sort(unique(predicted_feats$feature))

pred_chrs_path <- file.path(
  snakemake@params$output,
  "objects",
  "pred_chrs.rds"
)

cli::cli_alert_info(
  paste0(
    "Generating chromatograms for annotated ",
    "biotransformer predicted features"
  )
)
pred_chrs <- xcms::featureChromatograms(
  BPPARAM = BiocParallel::SerialParam(),
  chunkSize = 1L,
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
  file = pred_chrs_path
)
cli::cli_alert_success(
  paste0(
    "Saved biotransformer predicted feature chromatograms: ",
    "{.path {pred_chrs_path}}"
  )
)

# ==============================================================================
# Snakesave --------------------------------------------------------------------
# ==============================================================================
saveRDS(
  object = list(
    predicted_feats = predicted_feats,
    pred_chrs = pred_chrs,
    pred_peak_ids = pred_peak_ids
  ),
  file = snakemake@output[[1]]
)

script_footer()
end_log()