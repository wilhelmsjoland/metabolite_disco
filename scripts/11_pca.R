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
    library(ggplot2)
    library(dplyr)
    library(patchwork)
    library(RSQLite)
    library(MsBackendSql)
  })
)

limma_data <- readRDS(snakemake@input[["limma"]])
setup <- readRDS(snakemake@input[["setup"]])

intensities_mat <- limma_data$intensities_mat
meta <- setup$meta
group_colors <- setup$group_colors

# ==============================================================================
# PCA - PC1 & PC2 - before and after technical normalization -------------------
# ==============================================================================
cli::cli_h3("Generating principal components analyses")

# Data before median normalization
pca_res <- prcomp(
  x = t(intensities_mat$raw_fill_imp$log2_scale),
  scale = FALSE,
  center = FALSE
)

# Data after median normalization
pca_res_norm <- prcomp(
  x = t(intensities_mat$norm_fill_imp$log2_scale),
  scale = FALSE,
  center = FALSE
)
# ==============================================================================
# PCA - PC1 & PC2 - before and after technical normalization -------------------
# ==============================================================================
pca_raw <- plot_pca(
  prcomp_res = pca_res,
  metad = meta,
  x = PC1,
  y = PC2
) +
  ggplot2::labs(title = "Before median scaling - log2 transformed")

pca_adj <- plot_pca(
  prcomp_res = pca_res_norm,
  metad = meta,
  x = PC1,
  y = PC2
) +
  ggplot2::labs(title = "After median scaling - log2 transformed")

norm_filled_12_pca_p <- pca_raw / pca_adj +
  patchwork::plot_layout(guides = "collect")

norm_filled_12_pca_p_path <- file.path(
  snakemake@params$output,
  "graphs",
  "pca",
  "norm_filled_pca_1_2.pdf"
)

ggplot2::ggsave(
  filename = norm_filled_12_pca_p_path,
  plot = norm_filled_12_pca_p,
  device = "pdf",
  height = 6,
  width = 6,
  units = "in"
)

cli::cli_alert_success(
  paste0(
    "Saved PCA of PC1 & PC2 before and after median scaling to: ",
    "{.path {norm_filled_12_pca_p_path}}"
  )
)

# ==============================================================================
# PCA - PC3 & PC4 - before and after technical normalization -------------------
# ==============================================================================

# PC 3-4 before & after normalization
pca_raw <- plot_pca(
  prcomp_res = pca_res,
  metad = meta,
  x = PC3,
  y = PC4
) +
  ggplot2::labs(title = "Before median scaling - log2 transformed")

pca_adj <- plot_pca(
  prcomp_res = pca_res_norm,
  metad = meta,
  x = PC3,
  y = PC4
) +
  ggplot2::labs(title = "After median scaling - log2 transformed")

norm_filled_34_pca_p  <- pca_raw / pca_adj +
  patchwork::plot_layout(guides = "collect")

norm_filled_34_pca_p_path <- file.path(
  snakemake@params$output,
  "graphs",
  "pca",
  "norm_filled_pca_3_4.pdf"
)

ggplot2::ggsave(
  filename = norm_filled_34_pca_p_path,
  plot = norm_filled_34_pca_p,
  device = "pdf",
  height = 6,
  width = 6,
  units = "in"
)

cli::cli_alert_success(
  paste0(
    "Saved PCA of PC3 & PC4 before and after median scaling to: ",
    "{.path {norm_filled_34_pca_p_path}}"
  )
)

# ==============================================================================
# Snakesave --------------------------------------------------------------------
# ==============================================================================
saveRDS(
  object = list(),
  file = snakemake@output[[1]]
)

script_footer()
end_log()