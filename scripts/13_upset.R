# ==============================================================================
# Produce upset plots ----------------------------------------------------------
# ==============================================================================
message(
  "===========================================================================",
  "\n",
  "Generating upset plots ----------------------------------------------------",
  "\n",
  "==========================================================================="
)

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
        . < opt$qvalue,
        TRUE,
        FALSE
      )
    )
  )

upset_p <- upset_tib %>%
  ComplexUpset::upset(
    intersect = comparisons,
    name = paste0("Features with p adjusted < ", opt$qvalue),
    width_ratio = 0.15,
    base_annotations = list(
      "Intersecting features" = ComplexUpset::intersection_size()
    )
  )

ggplot2::ggsave(
  filename = paste0(opt$output, "/graphs/upset/upset_pdf"),
  plot = upset_p,
  device = "pdf",
  height = 7,
  width = 22,
  units = "in"
)