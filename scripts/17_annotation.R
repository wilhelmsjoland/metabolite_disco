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
    library(AnnotationHub)
    library(MetaboAnnotation)
    library(xcms)
    library(MetaboCoreUtils)
    library(ProtGenerics)
  })
)

limma_data <- readRDS(snakemake@input[["limma"]])
prep_data <- readRDS(snakemake@input[["prep_biot"]])
filter_data <- readRDS(snakemake@input[["filter_features"]])

intensities <- limma_data$intensities
all_sig_diff <- prep_data$all_sig_diff
xchr9 <- filter_data$xchr9

# ==============================================================================
# Match m/z's against databases ------------------------------------------------
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
  cli::cli_progress_done()

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

  cli::cli_progress_step("Matching features")
  matches <- MetaboAnnotation::matchValues(
    query = peaks_used,
    target = target_df,
    param = mz_match_param
  )
  cli::cli_progress_done()

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

anno_peak_ids <- sort(unique(anno$peak_id))

anno_chrs_path <- file.path(
  snakemake@params$output,
  "objects",
  "anno_chrs.rds"
)
if (interactive() && file.exists(anno_chrs_path)) {
  anno_chrs <- readRDS(anno_chrs_path)
  cli::cli_alert_success(
    paste0(
      "Imported database-matched feature chromatograms from: ",
      "{.path {anno_chrs_path}}"
    )
  )
} else {
  cli::cli_alert_info(
    paste0(
      "Generating chromatograms for annotated ",
      "database-matched features"
    )
  )
  anno_chrs <- xcms::featureChromatograms(
    BPPARAM = BiocParallel::SerialParam(),
    chunkSize = 1L,
    object = xchr9,
    expandRt = 0,
    expandMz = 0,
    aggregationFun = "sum",
    filled = TRUE,
    features = anno_peak_ids,
    missing = 0,
    return.type = "XChromatograms"
  )
  saveRDS(
    object = anno_chrs,
    file = anno_chrs_path
  )
  cli::cli_alert_success(
    paste0(
      "Saved annotated database-matched feature chromatograms: ",
      "{.path {anno_chrs_path}}"
    )
  )
}

# ==============================================================================
# Snakesave --------------------------------------------------------------------
# ==============================================================================
saveRDS(
  object = list(
    anno = anno,
    anno_chrs = anno_chrs,
    anno_peak_ids = anno_peak_ids
  ),
  file = snakemake@output[[1]]
)
end_log()