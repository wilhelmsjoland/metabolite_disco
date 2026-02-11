
hits <- filt_match_diffs2 %>%
  dplyr::mutate(obs_ppm = num_to_ppm(abs(delta_mass - obs_delta_mass))) %>%
  dplyr::filter(feat1 != feat2) %>%
  dplyr::arrange(obs_ppm)

hits_feats <- unique(c(hits$feat1, hits$feat2))

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
    res.folder,
    "graphs",
    "saved_pairs",
    "saved_pairs.tsv"
  )
)) {
  int_pairs <- read_tsv(
    file = file.path(res.folder, "graphs", "saved_pairs", "saved_pairs.tsv")
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
#   file = file.path(res.folder, "graphs", "saved_pairs", "saved_pairs.RDS")
# )
# write_tsv(
#   x = int_pairs,
#   file = file.path(res.folder, "graphs", "saved_pairs", "saved_pairs.tsv")
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
      .fns = ~ .x %in% pot.glycosides
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
  col = list(group = group.colors)
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
  col = list(group = group.colors)
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


# =============================================================================
# Molecular similarity ---------------------------------------------------------
# =============================================================================
# from the heatmap
int_bio_fts <- c(
  "FT11675",
  "FT08181",
  "FT08191",
  "FT04076",
  "FT05019",
  "FT05465",
  "FT02089",
  "FT04760",
  "FT01640",
  "FT06679"
)

int_annos <- c(
  "FT02161",
  "FT00292",
  "FT03755",
  "FT00485",
  "FT01696",
  "FT06069",
  "FT01657",
  "FT02948",
  "FT03779",
  "FT02246",
  "FT03341",
  "FT02024",
  "FT09148",
  "FT04929",
  "FT03009",
  "FT01252",
  "FT00081",
  "FT05181",
  "FT04812",
  "FT04750", # int no prod in 2 groups
  "FT01901",
  "FT01345",
  "FT05216",
  "FT01428",
  "FT00099",
  "FT05532",
  "FT01840",
  "FT01655",
  "FT02539",
  "FT05753",
  "FT00762",
  "FT00049",
  "FT02718", # weird
  "FT05496",
  "FT05809",
  "FT02012",
  "FT00307",
  "FT07141",
  "FT03793",
  "FT01327",
  "FT01404",
  "FT01658",
  "FT00183",
  "FT03331", # int
  "FT00157",
  "FT03095",
  "FT00103",
  "FT05777",
  "FT08395",
  "FT01133",
  "FT02027",
  "FT03826",
  "FT06714",
  "FT00372",
  "FT08744",
  "FT06988", # only in controls
  "FT02994",
  "FT07149", # only in mutant
  "FT01879", # super low in mutant apiin
  "FT00475",
  "FT06113",
  "FT04384",
  "FT03101", # int
  "FT03999",
  "FT02109", # almost only in mutant apiin
  "FT06809",
  "FT07476",
  "FT02952",
  "FT04898", # none in mutant apiin
  "FT10427",
  "FT04154",
  "FT06532", # int, strong peak
  "FT02265",
  "FT03040",
  "FT01936",
  "FT02266",
  "FT06508",
  "FT01385",
  "FT00397",
  "FT01640",
  "FT04480",
  "FT06507", # int, big difference
  "FT03435",
  "FT06584", # int - almost only in controls
  "FT07425", # int
  "FT02898",
  "FT05339",
  "FT01151", # only in mutant apiin?
  "FT02134",
  "FT09648",
  "FT08147",
  "FT09207",
  "FT01121",
  "FT04548",
  "FT01150", # weird
  "FT05943",
  "FT01119",
  "FT04168",
  "FT07026", # int only in controls
  "FT05666", # int only in controls
  "FT05782", # int only in controls
  "FT07752",
  "FT07522",
  "FT06533",
  "FT06244", # int only in controls
  "FT05874", # much higher in mutants
  "FT10419", # much higher in only mutant apiin
  "FT09078",
  "FT11262", # int only in controls
  "FT04555",
  "FT06998", # int much higher in controls
  "FT10557", # int much higher in controls
  "FT06635", # int higher in controls
  "FT06935", # int higher in controls
  "FT07834",
  "FT09271", # int only in apiin
  "FT06876", # int only in apiin
  "FT10018", # int only in apiin
  "FT02104",
  "FT07269",
  "FT05914", # int, probably different peaks in mutants and wt
  "FT06738",
  "FT05991",
  "FT10696", # weird
  "FT06263",
  "FT08063",
  "FT09227",
  "FT07139",
  "FT06789",
  "FT03730",
  "FT10864",
  "FT07014",
  "FT07616",
  "FT08733", # int only in controls
  "FT09819", # int only in mutant control
  "FT10162",
  "FT06676",
  "FT10048",
  "FT10885",
  "FT09160",
  "FT09410", # int only in controls
  "FT11592",
  "FT08566",
  "FT09268", # int higher in controls
  "FT10851"
)

# missing smiles
anno$peak_id[is.na(anno$target_smiles)]

# 648 rows
# 1,592
anno_filt <- anno %>%
  # TODO
  # DO A more intelligent filtering than this based on
  # how much information is available in all the rows
  dplyr::group_by(adduct) %>% # Think this might have fixed the loss of same
  dplyr::distinct(target_inchikey, .keep_all = TRUE) %>%
  dplyr::filter(peak_id %in% int_annos) # int.bio.feats

smiles <- anno_filt$target_smiles
names(smiles) <- anno_filt$feature

apigenin_smiles <- "C1=CC(=CC=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)O"
apiin_smiles <- paste0(
  "C1[C@@]([C@H]([C@@H](O1)O[C@@H]2[C@H]([C@@H]([C@H]",
  "(O[C@H]2OC3=CC(=C4C(=C3)OC(=CC4=O)C5=CC=C(C=C5)O)O)CO)O)O)O)(CO)O"
)

apigenin_sims <- mol_similarity(
  query_smiles = apigenin_smiles,
  target_smiles = smiles,
  kekulise = TRUE, # parsing incorrect smiles with electrons
  omit_nulls = TRUE,
  fingerprint = "circular",
  circular_type = "ECFP6",
  method = "tanimoto"
) %>%
  dplyr::left_join(
    x = .,
    y = anno,
    by = "feature"
  ) %>%
  dplyr::filter(sim > 0.15)

# Get names for individual inchikey
inchi_ks <- paste0(apigenin_sims$target_inchikey, collapse = ",")
cmd <- paste0(
  "curl -s -d \"inchikey=", inchi_ks, "\" ",
  "\"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/property/",
  "Title,MolecularFormula,InChIKey,ExactMass,XLogP,TPSA/",
  "CSV\""
)
results <- system(cmd, intern = TRUE)
inchikey_map <- read_csv(file = paste0(results, collapse = "\n"))

apigenin_sims %>%
  dplyr::left_join(
    x = .,
    y = inchikey_map,
    by = c("target_inchikey" = "InChIKey")
  ) %>%
  dplyr::mutate(Title = forcats::fct_reorder( # target_name
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
    aes(fill = peak_id),
    stat = "summary",
    fun = "max",
    color = "black",
    position = ggplot2::position_dodge()
  ) +
  ggplot2::guides(x = ggplot2::guide_axis(angle = -45)) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(c(0, 0)),
    limits = c(0, 0.8)
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_text(angle = -90),
    legend.title = ggplot2::element_blank()
  ) +
  ggplot2::labs(y = "Tanimoto similarity")
