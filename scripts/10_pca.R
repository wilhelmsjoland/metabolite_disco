cli::cli_h1(basename(this.path::this.path()))
# ==============================================================================
# PCA - PC1 & PC2 - before and after technical normalization -------------------
# ==============================================================================
cli::cli_h3("Generating principal components analyses")

# Data before normalization
vals <- SummarizedExperiment::assay(res, "raw_fill_imp") %>%
  log2() %>%
  t() %>%
  scale(center = TRUE, scale = TRUE) %>%
  as.matrix(.)

# Data after normalization
vals_norm <- SummarizedExperiment::assay(res, "norm_fill_imp") %>%
  log2() %>%
  t() %>%
  scale(center = TRUE, scale = TRUE) %>%
  as.matrix(.)
# ==============================================================================
# PCA - PC1 & PC2 - before and after technical normalization -------------------
# ==============================================================================
pca_res <- prcomp(vals, scale = FALSE, center = FALSE)
pca_raw <- plot_pca(
  prcomp_res = pca_res,
  metad = meta,
  x = PC1,
  y = PC2
) +
  ggplot2::labs(title = "Before median scaling")

pca_res_norm <- prcomp(vals_norm, scale = FALSE, center = FALSE)
pca_adj <- plot_pca(
  prcomp_res = pca_res_norm,
  metad = meta,
  x = PC1,
  y = PC2
) +
  ggplot2::labs(title = "After median scaling")

norm_filled_12_pca_p <- pca_raw / pca_adj +
  patchwork::plot_layout(guides = "collect")

norm_filled_12_pca_p_path <- file.path(
  opt$output,
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
  ggplot2::labs(title = "Before median scaling")

pca_adj <- plot_pca(
  prcomp_res = pca_res_norm,
  metad = meta,
  x = PC3,
  y = PC4
) +
  ggplot2::labs(title = "After median scaling")

norm_filled_34_pca_p  <- pca_raw / pca_adj +
  patchwork::plot_layout(guides = "collect")

norm_filled_34_pca_p_path <- file.path(
  opt$output,
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