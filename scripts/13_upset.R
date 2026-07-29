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
    library(tidyr)
    library(tibble)
    library(ComplexHeatmap)
    library(RSQLite)
    library(MsBackendSql)
  })
)

limma_data <- readRDS(snakemake@input[["limma"]])
full_limma <- limma_data$full_limma
comparisons <- limma_data$comparisons


# ==============================================================================
# Generating upset plots -------------------------------------------------------
# ==============================================================================
cli::cli_h3("Generating upset plots")

upset_tib <- full_limma %>%
  dplyr::select(feature, contrast, adj.P.Val) %>%
  tidyr::pivot_wider(
    names_from = "contrast",
    values_from = "adj.P.Val"
  ) %>%
  dplyr::mutate(
    dplyr::across(
      .cols = 2:ncol(.),
      .fns = ~ dplyr::if_else(
        . < snakemake@params$qvalue,
        1, # TRUE
        0 # FALSE
      )
    )
  )

# ==============================================================================
# Generating distinct upset plot -----------------------------------------------
# ==============================================================================
# Distinct is all that are distinctly significant in this set
cli::cli_alert_info(
  paste0(
    "Generating intersecting upset plot"
  )
)
upset_intersect <- upset_tib %>%
  tibble::column_to_rownames(var = "feature") %>%
  ComplexHeatmap::make_comb_mat(mode = "intersect")

upset_intersect_p_path <- file.path(
  snakemake@params$output,
  "graphs",
  "upset",
  "upset_intersection.pdf"
)
if (interactive() && file.exists(upset_intersect_p_path)) {
  cli::cli_alert_info(
    paste0(
      "{.path {upset_intersect_p_path} already exists. Not overwriting.}"
    )
  )
} else {
  pdf(file = upset_intersect_p_path, width = 12, height = 8)
  produce_complex_upset(
    input = upset_intersect,
    comps = comparisons,
    qvalue = snakemake@params$qvalue
  )
  invisible(dev.off())
  cli::cli_alert_success(
    paste0(
      "Intersection upset plot for {.val {snakemake@params$gap_filling}}",
      " saved to: {.path {upset_intersect_p_path}}"
    )
  )
}

# ==============================================================================
# Generating distinct upset plot -----------------------------------------------
# ==============================================================================
# Intersect are all that are significant alone and together
cli::cli_alert_info(
  paste0(
    "Generating distinct upset plot"
  )
)
upset_distinct <- upset_tib %>%
  tibble::column_to_rownames(var = "feature") %>%
  ComplexHeatmap::make_comb_mat(mode = "distinct")

upset_distinct_p_path <- file.path(
  snakemake@params$output,
  "graphs",
  "upset",
  "upset_distinct.pdf"
)
if (interactive() && file.exists(upset_distinct_p_path)) {
  cli::cli_alert_info(
    paste0(
      "{.path {upset_distinct_p_path} already exists. Not overwriting.}"
    )
  )
} else {
  pdf(file = upset_distinct_p_path, width = 12, height = 8)
  produce_complex_upset(
    input = upset_distinct,
    comps = comparisons,
    qvalue = snakemake@params$qvalue
  )
  invisible(dev.off())
  cli::cli_alert_success(
    paste0(
      "Distinct upset plot for {.val {snakemake@params$gap_filling}}",
      " saved to: {.path {upset_distinct_p_path}}"
    )
  )
}

# ==============================================================================
# Snakesave --------------------------------------------------------------------
# ==============================================================================
saveRDS(
  object = list(
    upset_distinct = upset_distinct,
    upset_intersect = upset_intersect
  ),
  file = snakemake@output[[1]]
)
end_log()