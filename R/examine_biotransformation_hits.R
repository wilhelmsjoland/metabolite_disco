hits <- filt.match.diffs2 %>%
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
for (i in unique(gly.agly$feature)) {
  tmp_hits_chrs <- plotFeatChrInt(
    feature_chrom = hits_chrs,
    feature = i,
    method = "sum",
    value = "into",
    filled = TRUE,
    missing = "rowmin_half",
    msLevel = 1,
    save_loc = "/test/",
    device = "pdf",
    feat_pairs = FALSE
  )
  print(tmp_hits_chrs$combined)
  stop_loop <- readline("Enter for next, break for stop: ")
  if (stop_loop == "break") {
    break
  }
}

# Feature pair graphs
# TODO
# go through all of this tomorrow
# currently at hit 319
# add feat1: FT02089, feat2: FT06587 to int_pairs
# add feat1: FT08181, feat2: FT09844
# int_pairs <- tibble::tibble()
for (i in 319:nrow(hits)) { # seq_len(nrow(hits))
  ft_pair_p <- plotFeatPairs(
    feature_chrom = hits_chrs,
    filt.match.row = hits[i, ],
    method = "sum",
    value = "into",
    filled = TRUE,
    missing = 0,
    msLevel = 1,
    save_pairs_loc = "/test/",
    device = "pdf"
  )
  print(ft_pair_p)
  idx <- which(df$rpair_num == hits[i, ][["name"]])
  message(
    "Feature pair: ", i, " / ", nrow(hits),
    "\nFeature X ± ", df[idx, ]$name1, " <=> ", "Feature Y ± ", df[idx, ]$name2,
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
}

print(int_pairs)

# FT02089 is apigenin
# FT08181 is apiin
# FT08191 is likely apiin too
