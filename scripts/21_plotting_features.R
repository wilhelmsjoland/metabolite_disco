cli::cli_h1(basename(this.path::this.path()))
# ==============================================================================
# Plotting features ------------------------------------------------------------
# ==============================================================================
cli::cli_progress_step("Writing feature chromatograms to plots...")
# Features
plot_twenty_feats(
  chromatogram = feature_chrs,
  save_loc = "/graphs/features/"
)
cli::cli_progress_done()



message("Writing feature chromatograms and intensity boxplots...")
# feats_to_plot <- sort(unique(all_sig_diff))
# for (i in feats_to_plot) {
#   ft_p <- plot_feat_chrom_int(
#     feature_chrom = feature_chrs,
#     feature = i,
#     method = "sum",
#     value = "into",
#     filled = TRUE,
#     missing = 0,
#     ms_level = 1,
#     save_loc = "/graphs/feature_chromatogram_intensity/",
#     device = "pdf",
#     feat_pairs = FALSE
#   )
# }

message("Plotting feature pairs in filtered biotransformations...")
# TODO
# FIX this
# for (i in seq_len(nrow(xchr9_filt$biot_filt_sig_features_tib))) {
#   ft_pair_p <- plot_feature_pairs(
#     feature_chrom = feature_chrs,
#     filt_match_row = xchr9_filt$biot_filt_sig_features_tib[i, ],
#     method = "sum",
#     value = "into",
#     filled = TRUE,
#     missing = 0,
#     ms_level = 1,
#     save_pairs_loc = "/graphs/feature_pairs/",
#     device = "pdf"
#   )
# }

message(
  "Writing feature chromatograms and intensity boxplots ",
  "for glycosides/aglycones..."
)
# for (i in pot_glycosides) {
#   ft_p <- plot_feat_chrom_int(
#     feature_chrom = feature_chrs,
#     feature = i,
#     method = "sum",
#     value = "into",
#     filled = TRUE,
#     missing = 0,
#     ms_level = 1,
#     save_loc = "/graphs/glycoside/",
#     device = "pdf",
#     feat_pairs = FALSE
#   )
# }

message("Plotting glycoside/aglycone feature pairs biotransformations...")
# for (i in seq_len(nrow(filt_match_diffs2))) {
#   ft_pair_p <- plot_feature_pairs(
#     feature_chrom = feature_chrs,
#     filt_match_row = filt_match_diffs2[i, ],
#     method = "sum",
#     value = "into",
#     filled = TRUE,
#     missing = 0,
#     ms_level = 1,
#     save_pairs_loc = "/graphs/glycoside_feature_pairs/",
#     device = "pdf"
#   )
# }

message("Producing significant intersecting feature boxplots...")
# for (i in upset_comp) {
#   ft_p <- plot_feat_chrom_int(
#     feature_chrom = feature_chrs,
#     feature = i,
#     method = "sum",
#     value = "into",
#     filled = TRUE,
#     missing = 0,
#     ms_level = 1,
#     save_loc = "/graphs/feature_chromatogram_intensity/",
#     device = "pdf",
#     feat_pairs = FALSE
#   )
# }

message("Plotting annotated features after molecular similarity")
# anno_sims2 %>%
#   # TODO
#   # arbitrary for now
#   dplyr::filter(sim > 0.3) %>%
#   dplyr::mutate(
#     Title = forcats::fct_reorder( # target_name
#       .f = Title, # target_name
#       .x = sim,
#       .fun = "mean",
#       .desc = TRUE
#     )
#   ) %>%
#   ggplot2::ggplot(
#     ggplot2::aes(
#       x = Title, # target_name
#       y = sim
#     )
#   ) +
#   # because there are duplicates - just choose the best one
#   ggplot2::geom_col(
#     # ggplot2::aes(fill = peak_id),
#     stat = "summary",
#     fun = "max",
#     color = "black",
#     position = ggplot2::position_dodge()
#   ) +
#   ggplot2::guides(x = ggplot2::guide_axis(angle = -45)) +
#   ggplot2::scale_y_continuous(
#     expand = ggplot2::expansion(c(0, 0)),
#     limits = c(0, 1)
#   ) +
#   ggplot2::theme_bw() +
#   ggplot2::theme(
#     axis.title.x = ggplot2::element_blank(),
#     axis.title.y = ggplot2::element_text(angle = -90),
#     legend.title = ggplot2::element_blank()
#   ) +
#   ggplot2::labs(y = "Tanimoto similarity")

message("Plotting predicted features after molecular similarity")
# pred_sims %>%
#   dplyr::mutate(
#     met_id = forcats::fct_reorder( # target_name
#       .f = met_id, # target_name
#       .x = sim,
#       .fun = "mean",
#       .desc = TRUE
#     )
#   ) %>%
#   ggplot2::ggplot(
#     ggplot2::aes(
#       x = met_id, # target_name
#       y = sim
#     )
#   ) +
#   # because there are duplicates - just choose the best one
#   ggplot2::geom_col(
#     ggplot2::aes(fill = feature),
#     stat = "summary",
#     fun = "max",
#     color = "black",
#     position = ggplot2::position_dodge()
#   ) +
#   ggplot2::guides(x = ggplot2::guide_axis(angle = -45)) +
#   ggplot2::scale_y_continuous(
#     expand = ggplot2::expansion(c(0, 0)),
#     limits = c(0, 1)
#   ) +
#   ggplot2::theme_bw() +
#   ggplot2::theme(
#     axis.title.x = ggplot2::element_blank(),
#     axis.title.y = ggplot2::element_text(angle = -90),
#     legend.title = ggplot2::element_blank()
#   ) +
#   ggplot2::labs(y = "Tanimoto similarity")
