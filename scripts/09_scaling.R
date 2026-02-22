cli::cli_h1(basename(this.path::this.path()))
# ==============================================================================
# Median scaling ---------------------------------------------------------------
# ==============================================================================
cli::cli_h3("Median scaling data")

res <- xcms::quantify(
  object = xchr9,
  method = "sum",
  value = "into",
  intensity = "into",
  filled = FALSE,
  missing = NA,
  msLevel = 1
)

SummarizedExperiment::assays(res)$raw_fill <- xcms::featureValues(
  object = xchr9,
  method = "sum",
  value = "into",
  intensity = "into",
  filled = TRUE,
  missing = NA,
  msLevel = 1
)

SummarizedExperiment::assays(res)$raw_fill_imp <- xcms::featureValues(
  object = xchr9,
  method = "sum",
  value = "into",
  intensity = "into",
  filled = TRUE,
  missing = "rowmin_half",
  msLevel = 1
)

SummarizedExperiment::assays(res)$norm <- median_scale_base(
  res_obj = res,
  assay = "raw"
)
SummarizedExperiment::assays(res)$norm_fill <- median_scale_base(
  res_obj = res,
  assay = "raw_fill"
)
SummarizedExperiment::assays(res)$norm_fill_imp <- median_scale_base(
  res_obj = res,
  assay = "raw_fill_imp"
)

cli::cli_alert_success(
  paste0(
    "Data median scaled and stored in: ",
    "{.val {paste0('res$', SummarizedExperiment::assayNames(res))}}"
  )
)
