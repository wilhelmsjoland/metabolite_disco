# ==============================================================================
# Median scaling ---------------------------------------------------------------
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
    library(SummarizedExperiment)
  })
)

filter_features <- readRDS(snakemake@input[["filter_features"]])
xchr9 <- filter_features$xchr9

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

assay_names <- SummarizedExperiment::assayNames(res)
cli::cli_bullets(
  c(
    "i" = paste0(
      "Unscaled dataframes stored in",
      "{.val {paste0('res$', assay_names[!grepl('norm', assay_names)])}}"
    ),
    "v" = paste0(
      "Median scaled dataframes stored in",
      "{.val {paste0('res$', assay_names[grepl('norm', assay_names)])}}"
    )
  )
)

# ==============================================================================
# Snakesave --------------------------------------------------------------------
# ==============================================================================
saveRDS(
  object = list(res = res),
  file = snakemake@output[[1]]
)

end_log()