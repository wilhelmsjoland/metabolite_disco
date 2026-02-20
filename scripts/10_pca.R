# ==============================================================================
# PCA ------------------------------------------------------
# ==============================================================================
message(
  "===========================================================================",
  "\n",
  "Generating PCAs -----------------------------------------------------------",
  "\n",
  "==========================================================================="
)


# Data before normalization
vals <- SummarizedExperiment::assay(res, "raw_filled") %>%
  log2() %>%
  t() %>%
  scale(center = TRUE, scale = TRUE) %>%
  as.matrix(.)

pca_res <- prcomp(vals, scale = FALSE, center = FALSE)
pca_raw <- plot_pca(
  prcomp_res = pca_res,
  metad = meta,
  x = PC1,
  y = PC2
) +
  ggplot2::labs(title = "Before normalization")

# Data after normalization
vals_norm <- SummarizedExperiment::assay(res, "norm_filled") %>%
  log2() %>%
  t() %>%
  scale(center = TRUE, scale = TRUE) %>%
  as.matrix(.)

pca_res_norm <- prcomp(vals_norm, scale = FALSE, center = FALSE)
pca_adj <- plot_pca(
  prcomp_res = pca_res_norm,
  metad = meta,
  x = PC1,
  y = PC2
) +
  ggplot2::labs(title = "After normalization")

norm_filled_12_pca_p <- pca_raw / pca_adj +
  patchwork::plot_layout(guides = "collect")

ggplot2::ggsave(
  filename = paste0(opt$output, "/graphs/pca/norm_filled_pca_1_2.pdf"),
  plot = norm_filled_12_pca_p,
  device = "pdf",
  height = 6,
  width = 6,
  units = "in"
)

# PC 3-4 before & after normalization
pca_raw <- plot_pca(
  prcomp_res = pca_res,
  metad = meta,
  x = PC3,
  y = PC4
) +
  ggplot2::labs(title = "Before normalization")

pca_adj <- plot_pca(
  prcomp_res = pca_res_norm,
  metad = meta,
  x = PC3,
  y = PC4
) +
  ggplot2::labs(title = "After normalization")

norm_filled_34_pca_p  <- pca_raw / pca_adj +
  patchwork::plot_layout(guides = "collect")

ggplot2::ggsave(
  filename = paste0(opt$output, "/graphs/pca/norm_filled_pca_3_4.pdf"),
  plot = norm_filled_34_pca_p,
  device = "pdf",
  height = 6,
  width = 6,
  units = "in"
)