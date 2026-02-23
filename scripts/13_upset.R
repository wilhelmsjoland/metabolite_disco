cli::cli_h1(basename(this.path::this.path()))
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

upset_p_path <- file.path(
  opt$output,
  "graphs",
  "upset",
  "upset.pdf"
)

if (file.exists(upset_p_path)) {
  cli::cli_alert_danger(
    paste0(
      "{.path {upset_p_path} already exists. Not overwriting.}"
    )
  )
} else {
  ggplot2::ggsave(
    filename = upset_p_path,
    plot = upset_p,
    device = "pdf",
    height = 7,
    width = 22,
    units = "in"
  )
}

cli::cli_alert_success(
  paste0(
    "Upset plot for {.val {opt$gap_filling}} saved to: ",
    "{.path {upset_p_path}}"
  )
)
