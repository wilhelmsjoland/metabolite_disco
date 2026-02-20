# ==============================================================================
# Median scaling & PCA ------------------------------------------------------
# ==============================================================================
message(
  "===========================================================================",
  "\n",
  "Generating PCAs -----------------------------------------------------------",
  "\n",
  "==========================================================================="
)

res <- xcms::quantify(
  xchr9,
  method = "sum",
  value = "into",
  filled = FALSE,
  missing = "rowmin_half" # 0 ?
)

SummarizedExperiment::assays(res)$raw_filled <- xcms::featureValues(
  xchr9,
  method = "sum",
  value = "into",
  filled = TRUE,
  missing = "rowmin_half" # 0 ?
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
  SummarizedExperiment::assay(res, "raw_filled"),
  MARGIN = 2,
  median,
  na.rm = TRUE
)
nf_mdn <- mdns / median(mdns)

# Dividing dataset by median of median and creating a new assay
SummarizedExperiment::assays(res)$norm_filled <- sweep(
  SummarizedExperiment::assay(res, "raw_filled"),
  MARGIN = 2,
  nf_mdn,
  "/"
)