# ==============================================================================
# Matching m/z's against databases -------------------------------------------
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
    library(AnnotationHub)
    library(MetaboAnnotation)
    library(MetaboCoreUtils)
    library(ProtGenerics)
  })
)

limma_data <- readRDS(snakemake@input[["limma"]])
prep_data <- readRDS(snakemake@input[["prep_biot"]])

intensities <- limma_data$intensities
all_sig_diff <- prep_data$all_sig_diff


# ==============================================================================
# Matching m/z's against databases -------------------------------------------
# ==============================================================================
cli::cli_h3(
  paste0(
    "Annotating features by matching m/z's against databases"
  )
)

peaks_used <- intensities$norm_fill_imp$untransformed %>%
  dplyr::select(feature, mzmed, rtmed) %>%
  dplyr::rename("mz" = "mzmed", "rtime" = "rtmed") %>%
  # All comps with at least one significant difference
  dplyr::filter(feature %in% all_sig_diff) %>%
  tibble::column_to_rownames(var = "feature")

peaks_used$peak_id <- rownames(peaks_used) # keep XCMS peak IDs

anno_path <- file.path(
  snakemake@params$output,
  "tables",
  "mz_annotations.csv"
)
if (interactive() && file.exists(anno_path)) {
  anno <- readr::read_csv(
    file = anno_path,
    show_col_types = FALSE,
    progress = FALSE
  )
  cli::cli_bullets(
    c(
      "i" = "Feature annotations already exist",
      "v" = "Imported feature annotatations from: {.path {anno_path}}"
    )
  )
} else {
  annotation_hub <- AnnotationHub::AnnotationHub()
  ah_id <- "AH119519"
  ah_metadata <- mcols(annotation_hub[ah_id]) %>%
    as.list() %>%
    purrr::map(
      .x = .,
      .f = ~ as.character(.x)
    ) %>%
    tibble::enframe() %>%
    tidyr::unnest(value)

  cli::cli_alert_info("Metadata for: {.val {ah_id}}")
  purrr::walk2(
    .x = ah_metadata$name,
    .y = ah_metadata$value,
    .f = ~ cli::cli_alert_info("{.x}: {.y}")
  )

  cli::cli_progress_step("Importing annotation database")
  cdb <- annotation_hub[[ah_id]]
  cli::cli_progress_done("Imported annotation database")

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
    adducts = c(
      MetaboCoreUtils::adductNames(polarity = snakemake@params$polarity)
    ),
    ppm = snakemake@params$ppm_match
  )

  cli::cli_bullets(
    c(
      "i" = "Matching features to database with: ",
      "i" = "polarity: {.val {snakemake@params$polarity}}",
      "i" = "ppm: {.val {snakemake@params$ppm_match}}"
    )
  )

  matches <- MetaboAnnotation::matchValues(
    query = peaks_used,
    target = target_df,
    param = mz_match_param
  )

  anno <- MetaboAnnotation::matchedData(matches) %>%
    tibble::as_tibble(x = ., rownames = "feature") %>%
    dplyr::mutate(abs_score = abs(score)) %>%
    dplyr::arrange(abs_score)

  readr::write_csv(x = anno, file = anno_path)
  cli::cli_alert_success(
    paste0(
      "Feature annotatations saved to :",
      "{.path {anno_path}}"
    )
  )
}

# ==============================================================================
# Snakesave --------------------------------------------------------------------
# ==============================================================================
saveRDS(
  object = list(
    anno = anno
  ),
  file = snakemake@output[[1]]
)

end_log()