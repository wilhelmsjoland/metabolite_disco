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
    library(xcms)
    library(tibble)
    library(dplyr)
    library(tidyr)
    library(RSQLite) # for sql backend
    library(MsBackendSql) # for sql backend
  })
)

correspondence <- readRDS(snakemake@input[["correspondence"]])
xchr7 <- correspondence$xchr7


# ==============================================================================
# Gap filling ------------------------------------------------------------------
# ==============================================================================
cli::cli_h3("Filling gaps")

xchr8_path <- file.path(
  snakemake@params$output,
  "objects",
  "xchr8.rds"
)

if (interactive() && file.exists(xchr8_path)) {
  xchr8 <- readRDS(file = xchr8_path)
  cli::cli_alert_success(
    paste0(
      "Imported gap filled object from: ",
      "{.path {xchr8_path}}"
    )
  )
} else {
  cli::cli_alert_info("Filling gaps")
  xchr8 <- xcms::fillChromPeaks(
    object = xchr7,
    BPPARAM = BiocParallel::bpparam(),
    chunkSize = snakemake@params$cores,
    param = xcms::ChromPeakAreaParam()
  )
  saveRDS(object = xchr8, file = xchr8_path)
  cli::cli_alert_success(
    paste0(
      "Saved gap filled object to: ",
      "{.path {xchr8_path}}"
    )
  )
}

# ==============================================================================
# Extract information on features with NAs prior to gap filling ----------------
# ==============================================================================
na_feat_samp_n <- featureValues(xchr7) %>%
  tibble::as_tibble(rownames = "feature") %>%
  dplyr::filter(
    dplyr::if_any(
      .cols = dplyr::contains(".mzML"),
      .fns = ~ is.na(.)
    )
  ) %>%
  nrow()
total_feat_samp_n <- nrow(featureValues(xchr7))
cli::cli_alert_warning(
  paste0(
    "Features with at least one NA across samples prior to gap filling: ",
    "{.val {na_feat_samp_n}}, total: {.val {total_feat_samp_n}}"
  )
)

# Number of missing entries
na_feat_n <- sum(is.na(featureValues(xchr7)))
# Total number of entries
total_feat_n <- nrow(featureValues(xchr7)) * ncol(featureValues(xchr7))
cli::cli_alert_warning(
  paste0(
    "Total number of NA entries across samples prior to gap filling: ",
    "{.val {na_feat_n}}, total: {.val {total_feat_n}}"
  )
)

# ==============================================================================
# Extract information on features with NAs after gap filling -------------------
# ==============================================================================
# Number of missing entries after gap filling
filled_na_feat_n <- sum(is.na(featureValues(xchr8)))
# Total number of entries after gap filling
filled_total_feat_n <- nrow(featureValues(xchr8)) * ncol(featureValues(xchr8))
# Total number of filled peaks after gap filling
filled_peak_n <- abs(filled_na_feat_n - na_feat_n)
cli::cli_alert_success("Filled {.val {filled_peak_n}} peaks")
cli::cli_alert_info(
  paste0(
    "Total number of NA entries across samples after gap filling: ",
    "{.val {filled_na_feat_n}}, total: {.val {filled_total_feat_n}}"
  )
)

# ==============================================================================
# Create chromatograms for all features with NA after gap filling ..............
# ==============================================================================
# Extract the m/z - rt regions for these features
# Extract features with nas for peak filling
# feat_with_na_after <- xcms::featureValues(
#   object = xchr8,
#   method = "sum",
#   value = "into",
#   intensity = "into",
#   filled = TRUE,
#   missing = NA,
#   msLevel = 1
# ) %>%
#   tibble::as_tibble(., rownames = "feature") %>%
#   tidyr::pivot_longer(cols = contains(".mzML")) %>%
#   dplyr::filter(is.na(value)) %>%
#   dplyr::pull(feature) %>%
#   unique(.)

# chrs_na_feat <- xcms::featureArea(
#   object = xchr8,
#   features = feat_with_na_after
# )

# # Expand the retention time by 1 second on both sides
# chrs_na_feat[, "rtmin"] <- chrs_na_feat[, "rtmin"] - 1
# chrs_na_feat[, "rtmax"] <- chrs_na_feat[, "rtmax"] + 1

# # For later plotting of non-filled peaks
# chrs_na_path <- file.path(
#   snakemake@params$output,
#   "objects",
#   "chrs_na.rds"
# )

# if (interactive() && file.exists(chrs_na_path)) {
#   chrs_na <- readRDS(file = chrs_na_path)
#   cli::cli_alert_success(
#     paste0(
#       "Imported saved chromatograms for NAs from: ",
#       "{.path {chrs_na_path}}"
#     )
#   )
# } else {
#   cli::cli_alert_info("Generating chromatograms for NAs")
#   chrs_na <- xcms::chromatogram(
#     BPPARAM = BiocParallel::bpparam(),
#     chunkSize = snakemake@params$cores,
#     object = xchr8,
#     mz = chrs_na_feat[, c("mzmin", "mzmax")],
#     # mabye increase this a little?
#     rt = chrs_na_feat[, c("rtmin", "rtmax")]
#   )
#   saveRDS(object = chrs_na, file = chrs_na_path)
#   cli::cli_alert_success(
#     paste0(
#       "Saved chromatograms for NAs to: ",
#       "{.path {chrs_na_path}}"
#     )
#   )
# }

# ==============================================================================
# Snakesave --------------------------------------------------------------------
# ==============================================================================
saveRDS(
  object = list(xchr8 = xchr8),
  file = snakemake@output[[1]]
)

script_footer()
end_log()