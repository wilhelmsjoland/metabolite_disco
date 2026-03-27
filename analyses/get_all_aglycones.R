gly_agly <- readr::read_csv(
  paste0(
    output_path, "/",
    "/afzelin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon",
    "/tables/gly_agly.csv"
  ),
  progress = FALSE,
  show_col_types = FALSE
)

chrs <- xcms::featureChromatograms(
  BPPARAM = MulticoreParam(workers = 4),
  chunkSize = 4L,
  object = xchr9,
  expandRt = 0,
  expandMz = 0,
  aggregationFun = "sum",
  filled = TRUE,
  features = gly_agly$feature,
  missing = 0,
  msLevel = 1L,
  return.type = "XChromatograms"
)
colnames(chrs) <- sub(".*[/\\\\]", "", colnames(chrs))

# Gly agly plotting
chr_data <- get_chr_data(chrs, gly_agly$feature, meta) %>%
  dplyr::left_join(
    x = .,
    y = meta,
    by = c("sample" = "sample")
  )
for (i in gly_agly$feature) {
  tmp_p <- chr_data %>%
    dplyr::filter(feature == i) %>%
    dplyr::mutate(intensity = tidyr::replace_na(intensity, 0)) %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = rtime,
        y = intensity,
        color = group
      )
    ) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~sample, scales = "fixed") +
    ggplot2::scale_color_manual(values = group_colors) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_blank(),
      legend.position = "bottom"
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 2)) +
    labs(title = i)

  print(tmp_p)
  readline("Enter for next: ")
}

# Extract the glycoside and aglycone from each gly_agly
result <- plot_feature(
  feature_chrom = chrs,
  feature = "FT09024",
  meta = meta,
  limma_results = full_limma,
  method = "sum",
  value = "into",
  filled = TRUE,
  missing = 0,
  ms_level = 1L,
  save_loc = NULL,
  device = "pdf",
  overwrite = FALSE
)
result$full
result$half
result$overlay