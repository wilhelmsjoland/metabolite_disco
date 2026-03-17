cli::cli_h1(basename(this.path::this.path()))
# ==============================================================================
# Plotting features ------------------------------------------------------------
# ==============================================================================
cli::cli_progress_step("Writing feature chromatograms to plots...")
# Features
plot_twenty_feats(
  chromatogram = feature_chrs,
  save_loc = file.path(opt$output, "graphs", "features")
)
cli::cli_progress_done()

# TODO
# Choose features to plot from:
# 1. all the chem sim 1 & 2
# 2. biotransformer.jar,
# 3. biotransformations
# 4. anno
# placeholder
# feats_to_plot <- sort(unique(all_sig_diff))
# placeholder
cli::cli_progress_step(
  "Writing feature chromatograms and intensity boxplots"
)
feats_to_plot <- rownames(xcms::featureDefinitions(feature_chrs))[1:10]
for (i in feats_to_plot) {
  ft_p <- plot_feat_chrom_int(
    feature_chrom = feature_chrs,
    feature = i,
    method = "sum",
    value = "into",
    filled = TRUE,
    missing = 0,
    ms_level = 1,
    save_loc = "/graphs/feature_chromatogram_intensity/",
    device = "pdf",
    feat_pairs = FALSE,
    overwrite = FALSE
  )
}
cli::cli_progress_done()

cli::cli_progress_step("Plotting feature pairs in filtered biotransformations")
# TODO
# FIX the feature vector to plot with
# Also make sure that i dont plot the same features but with
# different adducts or anything like that
# for (i in seq_len(nrow(xchr9_filt$biot_filt_sig_features_tib))) <- old

test <- subset_matched_diffs %>%
  dplyr::filter(
    dplyr::if_all(
      .cols = dplyr::all_of(c("feat1", "feat2")),
      .fns = ~ .x %in% c(
        "FT02088", "FT02089", "FT08181", "FT02223",
        "FT02409", "FT02839", "FT02925", "FT03582",
        "FT04452", "FT06025"
      )
    )
  )

for (i in seq_len(nrow(test[1:5, ]))) {
  ft_pair_p <- plot_feature_pairs(
    feature_chrom = feature_chrs,
    filt_match_row = test[i, ],
    method = "sum",
    value = "into",
    filled = TRUE,
    missing = 0,
    ms_level = 1,
    save_pairs_loc = "/graphs/feature_pairs/",
    device = "pdf",
    overwrite = FALSE
  )
}
cli::cli_progress_done()

cli::cli_progress_step(
  paste0(
    "Writing feature chromatograms and intensity boxplots ",
    "for glycosides/aglycones..."
  )
)
for (i in pot_glycosides) {
  ft_p <- plot_feat_chrom_int(
    feature_chrom = feature_chrs,
    feature = i,
    method = "sum",
    value = "into",
    filled = TRUE,
    missing = 0,
    ms_level = 1,
    save_loc = "/graphs/glycoside/",
    device = "pdf",
    feat_pairs = FALSE,
    overwrite = FALSE
  )
}
cli::cli_progress_done()

cli::cli_progress_step(
  "Plotting glycoside/aglycone feature pair biotransformations"
)
for (i in seq_len(nrow(test[1:5, ]))) {
  ft_pair_p <- plot_feature_pairs(
    feature_chrom = feature_chrs,
    filt_match_row = test[i, ],
    method = "sum",
    value = "into",
    filled = TRUE,
    missing = 0,
    ms_level = 1,
    save_pairs_loc = "/graphs/glycoside_feature_pairs/",
    device = "pdf",
    overwrite = FALSE
  )
}
cli::cli_progress_done()

cli::cli_progress_step("Producing significant intersecting feature boxplots...")
# Should work in the future but feature_chrs is wrong for now
for (i in upset_comp) {
  ft_p <- plot_feat_chrom_int(
    feature_chrom = feature_chrs,
    feature = i,
    method = "sum",
    value = "into",
    filled = TRUE,
    missing = 0,
    ms_level = 1,
    save_loc = "/graphs/feature_chromatogram_intensity/",
    device = "pdf",
    feat_pairs = FALSE,
    overwrite = FALSE
  )
}
cli::cli_progress_done()

message("Plotting annotated features after molecular similarity")
anno_sims_final %>%
  dplyr::filter(feature %in% unname(unlist(upset_distinct_comp))) %>%
  dplyr::filter(sim > 0.2) %>%
  dplyr::mutate(
    Title = forcats::fct_reorder( # target_name
      .f = Title, # target_name
      .x = sim,
      .fun = "mean",
      .desc = TRUE
    )
  ) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = Title, # target_name
      y = sim
    )
  ) +
  # because there are duplicates - just choose the best one
  ggplot2::geom_col(
    ggplot2::aes(fill = peak_id),
    stat = "summary",
    fun = "max",
    color = "black",
    position = ggplot2::position_dodge()
  ) +
  ggplot2::guides(x = ggplot2::guide_axis(angle = -45)) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(c(0, 0)),
    limits = c(0, 1)
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_text(angle = -90),
    legend.title = ggplot2::element_blank()
  ) +
  ggplot2::labs(y = "Tanimoto similarity")

message("Plotting predicted features after molecular similarity")
chem_pred_feats %>%
  dplyr::mutate(met_id = as.factor(met_id)) %>%
  dplyr::filter(feature %in% unname(unlist(upset_distinct_comp))) %>%
  dplyr::filter(sim > 0.2) %>%
  dplyr::mutate(
    met_id = forcats::fct_reorder( # target_name
      .f = met_id, # target_name
      .x = sim,
      .fun = "mean",
      .desc = TRUE
    )
  ) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = met_id, # target_name
      y = sim
    )
  ) +
  # because there are duplicates - just choose the best one
  ggplot2::geom_col(
    ggplot2::aes(fill = feature),
    stat = "summary",
    fun = "max",
    color = "black",
    position = ggplot2::position_dodge()
  ) +
  ggplot2::guides(x = ggplot2::guide_axis(angle = -45)) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(c(0, 0)),
    limits = c(0, 1)
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_text(angle = -90),
    legend.title = ggplot2::element_blank()
  ) +
  ggplot2::labs(y = "Tanimoto similarity")
