# Chemical similarity here
for (i in pred_peak_ids) {
  tmp_anno_chr <- plot_feat_chrom_int(
    feature_chrom = pred_chrs,
    feature = i,
    method = "sum",
    value = "into",
    filled = TRUE,
    missing = 0,
    ms_level = 1,
    save_loc = NULL,
    device = NULL,
    feat_pairs = FALSE
  )
  print(tmp_anno_chr$combined)
  stop_loop <- readline("Enter for next, break for end: ")
  if (stop_loop == "break") {
    break
  }
}

# FT06025
# FT08376

# REALLY INTERESTING!
to_look_at <- c(
  "FT03581", # 1
  "FT02413", # 2
  "FT08376",
  "FT09870",
  "FT06149",
  "FT04835"
)

# FIX the filtering first
test <- predicted_feats %>%
  dplyr::filter(feature %in% to_look_at) %>%
  dplyr::arrange(InChIKey) %>%
  dplyr::mutate(name = paste(feature, InChIKey)) %>%
  # only for now
  dplyr::distinct(InChIKey, .keep_all = TRUE) # %>%
  # only for now
 #  dplyr::filter(grepl("FT03581", feature))

for (i in test$SMILES) {
  mol.2d <- parse.smiles(i)[[1]]
  img <- view.image.2d(mol.2d)
  img.grob <- grid::rasterGrob(img, interpolate = TRUE)
  img.p <- wrap_elements(img.grob, clip = TRUE)
  print(img.p) 
    stop_loop <- readline("Enter for next, break for end: ")
    if (stop_loop == "break") {
      break
    }
}

biot.sims <- test$SMILES
names(biot.sims) <- test$name
apigenin.smiles <- "C1=CC(=CC=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)O"

biot.sim.final <- mol_similarity(
  query_smiles = apigenin.smiles,
  target_smiles = biot.sims,
  kekulise = TRUE,
  omit_nulls = TRUE,
  fingerprint = "circular",
  circular_type = "ECFP6",
  method = "tanimoto"
)

biot.sim.final %>%
  dplyr::left_join(
    x = .,
    y = test,
    by = c("feature" = "name")
  ) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = feature, # target_name
      y = sim,
    )
  ) +
  # because there are duplicates - just choose the best one
  ggplot2::geom_col(
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
