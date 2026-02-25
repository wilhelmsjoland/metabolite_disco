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
  opt$output,
  "graphs",
  "upset",
  "upset_intersection.pdf"
)
if (file.exists(upset_intersect_p_path)) {
  cli::cli_alert_info(
    paste0(
      "{.path {upset_intersect_p_path} already exists. Not overwriting.}"
    )
  )
} else {
  pdf(file = upset_intersect_p_path, width = 12, height = 8)
  produce_complex_upset(
    input = upset_intersect,
    comps = comparisons
  )
  invisible(dev.off())
  cli::cli_alert_success(
    paste0(
      "Upset plot for {.val {opt$gap_filling}} saved to: ",
      "{.path {upset_intersect_p_path}}"
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
  opt$output,
  "graphs",
  "upset",
  "upset_distinct.pdf"
)
if (file.exists(upset_distinct_p_path)) {
  cli::cli_alert_info(
    paste0(
      "{.path {upset_distinct_p_path} already exists. Not overwriting.}"
    )
  )
} else {
  pdf(file = upset_distinct_p_path, width = 12, height = 8)
  produce_complex_upset(
    input = upset_distinct,
    comps = comparisons
  )
  invisible(dev.off())
  cli::cli_alert_success(
    paste0(
      "Upset plot for {.val {opt$gap_filling}} saved to: ",
      "{.path {upset_distinct_p_path}}"
    )
  )
}


# Not correct -> upset_p changes order of stuff
extract_comb(m, "001110") %>% length
extract_comb(m, "010101") %>% length