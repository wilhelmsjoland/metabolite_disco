
hits <- subset_matched_diffs %>%
  dplyr::filter(feat1 != feat2) %>%
  dplyr::arrange(obs_diff)

hits_feats <- sort(unique(c(hits$feat1, hits$feat2)))

if (!exists("hits_chrs")) {
  hits_chrs <- xcms::featureChromatograms(
    object = xchr9,
    expandRt = 0,
    expandMz = 0,
    aggregationFun = "sum",
    filled = TRUE,
    features = hits_feats,
    missing = 0,
    return.type = "XChromatograms"
  )
} else {
  message("'hits_chrs' already exists")
}

# Individual graphs
# for (i in unique(gly.agly$feature)) {
#   tmp_hits_chrs <- plotFeatChrInt(
#     feature_chrom = hits_chrs,
#     feature = i,
#     method = "sum",
#     value = "into",
#     filled = TRUE,
#     missing = "rowmin_half",
#     msLevel = 1,
#     save_loc = NULL, # "/graphs/saved_pairs/",
#     device = "pdf", # "pdf"
#     feat_pairs = FALSE
#   )
#   print(tmp_hits_chrs$combined)
#   stop_loop <- readline("Enter for next, break for stop: ")
#   if (stop_loop == "break") {
#     break
#   }
# }

# Feature pair graphs
# ==============================================================================
# TODO -------------------------------------------------------------------------
# Fix the chemical similarity steps
# add feat1: FT02089, feat2: FT06587 to int_pairs
# add feat1: FT08181, feat2: FT09844
# add Feat1: FT08191, Feat2: FT10025
# add Feat1: FT02089, Feat2: FT02434
# This one is really interesting Feat1: FT02089, Feat2: FT01640
# ==============================================================================

if (file.exists(
  file.path(
    res_folder,
    "graphs",
    "saved_pairs",
    "saved_pairs.tsv"
  )
)) {
  int_pairs <- read_tsv(
    file = file.path(res_folder, "graphs", "saved_pairs", "saved_pairs.tsv")
  )
} else {
  # int_pairs <- tibble::tibble()
  for (i in seq_len(nrow(hits))) {
    # break
    ft_pair_p <- plotFeatPairs(
      feature_chrom = hits_chrs,
      filt.match.row = hits[i, ],
      method = "sum",
      value = "into",
      filled = TRUE,
      missing = 0,
      msLevel = 1,
      save_pairs_loc = "/graphs/saved_pairs/",
      device = "pdf"
    )
    print(ft_pair_p)
    idx <- which(df$rpair_num == hits[i, ][["name"]])
    message(
      "Feature pair: ", i, " / ", nrow(hits),
      "\nFeature X ± ", df[idx, ]$name1, " <=> ",
      "Feature Y ± ", df[idx, ]$name2,
      "\nFeat1: ", hits[i, ][["feat1"]], ", Feat2: ", hits[i, ][["feat2"]],
      "\nDelta_form: ", hits[i, ][["chem_change"]],
      "\nAdduct1: ", hits[i, ][["adduct1"]],
      ", Adduct2: ", hits[i, ][["adduct2"]],
      "\nObs ppm: ", hits[i, ][["obs_ppm"]],
      "\nRclass: ", df[idx, ][["entry"]],
      "\nRpair: ", df[idx, ][["rpair"]]
    )

    print(hits[i, ])
    stop_loop <- readline("Enter for next, 'break' for stop, 'save' for save: ")
    if (stop_loop == "break") {
      break
    }
    if (stop_loop == "save") {
      int_pairs <- dplyr::bind_rows(int_pairs, hits[i, ])
    }

    if (i %in% seq(from = 0, to = nrow(hits),  by = 15)) {
      dev.off()
    }
  }
}


# saveRDS(
#   object = int_pairs,
#   file = file.path(res_folder, "graphs", "saved_pairs", "saved_pairs.RDS")
# )
# write_tsv(
#   x = int_pairs,
#   file = file.path(res_folder, "graphs", "saved_pairs", "saved_pairs.tsv")
# )




# FT02089 is apigenin
# FT08181 is apiin
# FT08191 is likely apiin too


test <- int_pairs %>%
  dplyr::mutate(
    check_pairs = paste0(feat1, "_", feat2),
  ) %>%
  dplyr::distinct(check_pairs, .keep_all = TRUE)

test2 <- full_raw_filled %>% # full_norm_filled
  dplyr::select(
    -c(
      "mzmed",
      "mzmin",
      "mzmax",
      "rtmed",
      "rtmin",
      "rtmax",
      "npeaks",
      "bu_mutant_apiin",
      "bu_mutant_control",
      "bu_wt_apiin",
      "bu_wt_control",
      "ms_level"
    )
  ) %>%
  tidyr::pivot_longer(cols = contains(".mzML"))


# Do correlations for all data overall to see if there any
# interesting metabolites that seem matched
int_pair_feats <- sort(unique(c(int_pairs$feat1, int_pairs$feat2)))

test3 <- test2 %>%
  tidyr::pivot_wider(
    names_from = "feature",
    values_from = "value"
  ) %>%
  dplyr::left_join(
    x = .,
    y = as_tibble(meta, rownames = "name"),
    by = "name"
  ) %>%
  dplyr::relocate(c("group", "path"), .after = "name") %>%
  # now only for one group
  dplyr::filter(group == "bu_mutant_apiin") %>%
  dplyr::select(name, dplyr::all_of(int_pair_feats)) %>%
  # group_by(group_var) %>%
  summarize(
    results = list(
      {
        # data.frame of just the cols
        dat <- dplyr::pick(dplyr::all_of(int_pair_feats))
        purrr::map_dfr(
          .x = combn(int_pair_feats, 2, simplify = FALSE),
          .f = function(p) {
            broom::tidy(
              cor.test(
                dat[[p[1]]],
                dat[[p[2]]],
                method = "pearson"
              )
            ) %>%
              dplyr::mutate(var1 = p[1], var2 = p[2], .before = 1)
          }
        )
      }
    ),
    .groups = "drop"
  ) %>%
  tidyr::unnest(results)

test4 <- test3 %>%
  dplyr::filter(
    dplyr::if_any(
      .cols = dplyr::all_of(c("var1", "var2")),
      .fns = ~ .x %in% pot_glycosides
    )
  ) %>%
  rstatix::add_significance() %>%
  dplyr::filter(p.value.signif != "ns") %>%
  tidyr::drop_na(estimate) %>%
  dplyr::mutate(p.value.signif = if_else(
    p.value.signif == "ns",
    NA,
    as.character(p.value.signif)
  ))

test4 %>%
  ggplot(.,
    aes(
      x = var1,
      y = var2,
      fill = estimate
    )
  ) +
  geom_tile() +
  geom_text(aes(label = p.value.signif)) +
  scale_fill_gradient2() +
  guides(
    x = guide_axis(angle = -45)
  )

ind_int_fts <- sort(unique(c(test4$var1, test4$var2)))

for (i in ind_int_fts) {
  tmp_feat_corr <- plotFeatChrInt(
    feature_chrom = hits_chrs,
    feature = i,
    method = "sum",
    value = "into",
    filled = TRUE,
    missing = "rowmin_half",
    msLevel = 1,
    save_loc = NULL,
    device = "pdf",
    feat_pairs = FALSE
  )
  print(tmp_feat_corr$combined)

  stop_loop <- readline("Enter for next, 'break' for stop: ")
  if (stop_loop == "break") {
    break
  }
  if (i %in% seq(from = 0, to = nrow(hits),  by = 15)) {
    dev.off()
  }
}

# ==============================================================================
# ComplexHeatmap tibble
# ==============================================================================

# TODO
# Make these heatmaps for all features!!!!
# ALL FEATURES
# Add more database matching and better databases
# mabye even combine more databases

test5 <- full_norm_filled %>% # full_norm_filled
  dplyr::select(
    -c(
      "mzmed",
      "mzmin",
      "mzmax",
      "rtmed",
      "rtmin",
      "rtmax",
      "npeaks",
      "bu_mutant_apiin",
      "bu_mutant_control",
      "bu_wt_apiin",
      "bu_wt_control",
      "ms_level"
    )
  ) %>%
  tidyr::pivot_longer(cols = contains(".mzML")) %>%
  dplyr::filter(feature %in% ind_int_fts) %>%
  # dplyr::mutate(feature_sample = paste0(feature, "_", name)) %>%
  dplyr::left_join(
    x = .,
    y = tibble::as_tibble(meta, rownames = "name"),
    by = "name"
  )

# ==============================================================================
# 1) Clustering with correlation distances
# ==============================================================================
mat_cor <- test5 %>%
  dplyr::select(feature, name, value) %>%
  tidyr::pivot_wider(
    names_from = name,
    values_from = value,
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("feature") %>%
  as.matrix()

d_samp_cor <- as.dist(
  1 - cor(mat_cor,  use = "pairwise.complete.obs", method = "pearson")
)
hc_samp_cor <- hclust(d_samp_cor, method = "average")
name_order_cor <- colnames(mat_cor)[hc_samp_cor$order]

d_feat_cor <- as.dist(
  1 - cor(t(mat_cor), use = "pairwise.complete.obs", method = "pearson")
)
hc_feat_cor <- hclust(d_feat_cor, method = "average")
feature_order_cor <- rownames(mat_cor)[hc_feat_cor$order]

# reorder matrix to match hclust orders
mat_cor_ord <- mat_cor[feature_order_cor, name_order_cor, drop = FALSE]

ann_df_cor <- test5 %>%
  dplyr::distinct(name, group) %>%
  dplyr::right_join(tibble::tibble(name = colnames(mat_cor_ord)), by = "name")

ha_cor <- HeatmapAnnotation(
  group = ann_df_cor$group,
  annotation_name_side = "left",
  col = list(group = group_colors)
)

ht_cor <- Heatmap(
  mat_cor_ord,
  name = "value",
  col = circlize::colorRamp2(
    c(min(mat_cor_ord, na.rm = TRUE), 0, max(mat_cor_ord, na.rm = TRUE)),
    c("cornflowerblue", "white", "firebrick")
  ),
  top_annotation = ha_cor,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_rot = 45
)

draw(ht_cor, heatmap_legend_side = "right", annotation_legend_side = "right")

# =============================================================================
# Clustered ComplexHeatmap -----------------------------------------------------
# =============================================================================
mat_ht <- test5 %>%
  dplyr::select(feature, name, value) %>%
  tidyr::pivot_wider(
    names_from = name,
    values_from = value,
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("feature") %>%
  as.matrix()

ann_df <- test5 %>%
  dplyr::distinct(name, group) %>%
  dplyr::right_join(tibble::tibble(name = colnames(mat_ht)), by = "name")

ha <- HeatmapAnnotation(
  group = ann_df$group,
  annotation_name_side = "left",
  col = list(group = group_colors)
)

ht <- Heatmap(
  mat_ht,
  name = "value",
  col = circlize::colorRamp2(
    c(min(mat_ht, na.rm = TRUE), 0, max(mat_ht, na.rm = TRUE)),
    c("cornflowerblue", "white", "firebrick")
  ),
  top_annotation = ha,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_rot = 45
)

draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")