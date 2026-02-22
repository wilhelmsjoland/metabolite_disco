cli::cli_h1(basename(this.path::this.path()))
# ==============================================================================
# Median scaling ---------------------------------------------------------------
# ==============================================================================
cli::cli_h3("Median scaling data")

res <- xcms::quantify(
  xchr9,
  method = "sum",
  value = "into",
  filled = FALSE,
  missing = "rowmin_half" # 0 ? # dont impute
)

SummarizedExperiment::assays(res)$raw_filled <- xcms::featureValues(
  xchr9,
  method = "sum",
  value = "into",
  filled = TRUE,
  missing = "rowmin_half" # 0 ? # dont impute
)

# Compute median and generate normalization factor
mdns <- apply(
  SummarizedExperiment::assay(res, "raw"),
  MARGIN = 2,
  median,
  na.rm = TRUE
)
nf_mdn <- mdns / median(mdns)

# Dividing dataset by median of median and creating a new assay
SummarizedExperiment::assays(res)$norm <- sweep(
  SummarizedExperiment::assay(res, "raw"),
  MARGIN = 2,
  nf_mdn,
  "/"
)

# Compute median and generate normalization factor
mdns <- apply(
  X = SummarizedExperiment::assay(res, "raw_filled"),
  MARGIN = 2,
  FUN = function(x) {
    median(x, na.rm = TRUE)
  }
)
nf_mdn <- mdns / median(mdns)

# Dividing dataset by median of median and creating a new assay
SummarizedExperiment::assays(res)$norm_filled <- sweep(
  x = SummarizedExperiment::assay(res, "raw_filled"),
  MARGIN = 2,
  STATS = nf_mdn,
  FUN = "/"
)

cli::cli_alert_success(
  paste0(
    "Data median scaled and stored in: ",
    "{.val {paste0('res$', names(SummarizedExperiment::assays(res)))}}"
  )
)
