source("scripts/functions.R")

###############################################################################
# Setup output folder ----------------------------------------------------------
################################################################################
output_folders <- list.files(
  "/Volumes/bluecub/aglycone_release_100um_24h/output",
  full.names = TRUE
)

################################################################################
# Read all bio_sims ------------------------------------------------------------
################################################################################
bio_sim_paths <- file.path(
  output_folders,
  "tables",
  "biotransformer_similarities.csv"
)
bio_sim_paths <- bio_sim_paths[file.exists(bio_sim_paths)]

bio_sims <- purrr::map(
  .x = bio_sim_paths,
  .f = ~ {
    bio_sim <- readr::read_csv(
      file = .x,
      progress = FALSE,
      show_col_types = FALSE
    )

    sig_diffs <- readRDS(
       file.path(
         dirname(dirname(.x)),
        "snakemake_objects",
        "15_prep_annotation_biotransformation.rds"
      )
    )[["all_sig_diff"]]

    bio_sim %>%
      dplyr::filter(feature %in% sig_diffs)
  }
)
names(bio_sims) <- basename(dirname(dirname(bio_sim_paths)))
bio_sims2 <- purrr::keep(bio_sims, ~ nrow(.x) > 0)

all_bio_sims <- purrr::map_dfr(
  .x = bio_sims2,
  .f = ~ .x,
  .id = "experiment"
) %>%
  dplyr::relocate("feature", .after = "experiment")

################################################################################
# Read all anno_sims -----------------------------------------------------------
################################################################################
anno_sim_paths <- file.path(
  output_folders,
  "tables",
  "anno_similarities.csv"
)
anno_sim_paths <- anno_sim_paths[file.exists(anno_sim_paths)]

anno_sims <- purrr::map(
  .x = anno_sim_paths,
  .f = ~ {
    anno_sim <- readr::read_csv(
      file = .x,
      progress = FALSE,
      show_col_types = FALSE
    )

    sig_diffs <- readRDS(
      file.path(
        dirname(dirname(.x)),
        "snakemake_objects",
        "15_prep_annotation_biotransformation.rds"
      )
    )[["all_sig_diff"]]

    anno_sim %>%
      dplyr::filter(peak_id %in% sig_diffs)
  }
)
names(anno_sims) <- basename(dirname(dirname(anno_sim_paths)))
anno_sims2 <- purrr::keep(anno_sims, ~ nrow(.x) > 0)

all_anno_sims <- purrr::map_dfr(
  .x = anno_sims2,
  .f = ~ .x,
  .id = "experiment"
)

################################################################################
# Extract only interesting significant differences where G + Samp > Samp -------
################################################################################

# 1. define interesting comparisons
# 2. write out logic e.g. all sig in certain comps & p < 0.01
# 3. ascertain that i can extract the interesting ones by looking at chroms
# 4. plot heatmap somehow

extract_feats_peaks <- all_anno_sims %>%
  dplyr::filter(
    experiment == "afzelin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon"
  ) %>%
  dplyr::group_by(adduct) %>%
  dplyr::distinct(target_inchikey, .keep_all = TRUE) %>%
  dplyr::filter(sim > 0.3) %>%
  dplyr::pull(peak_id)

extract_bio_sims_peaks <- all_bio_sims %>%
  dplyr::filter(
    experiment == "afzelin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon"
  ) %>%
  dplyr::group_by(adduct) %>%
  dplyr::distinct(InChIKey, .keep_all = TRUE) %>%
  dplyr::filter(sim > 0.2) %>%
  dplyr::pull(feature)

both_tib_feats <- unique(c(extract_feats_peaks, extract_bio_sims_peaks))


extract_feats <- all_anno_sims %>%
  dplyr::filter(
    experiment == "afzelin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon"
  ) %>%
  dplyr::group_by(adduct) %>%
  dplyr::distinct(target_inchikey, .keep_all = TRUE) %>%
  dplyr::filter(sim > 0.3) %>%
  dplyr::filter(peak_id %in% both_tib_feats)

extract_bio_sims <- all_bio_sims %>%
  dplyr::filter(
    experiment == "afzelin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon"
  ) %>%
  dplyr::group_by(adduct) %>%
  dplyr::distinct(InChIKey, .keep_all = TRUE) %>%
  dplyr::filter(feature %in% both_tib_feats)


extract_feats$peak_id

extract_feats_int <- xcms::featureValues(
  object = test,
  method = "sum",
  value = "into",
  intensity = "into",
  filled = TRUE,
  missing = 0
) %>%
  tibble::as_tibble(., rownames = "feature") %>%
  dplyr::filter(feature %in% both_tib_feats) %>% # extract_feats$peak_id

extract_feats_int2 <- extract_feats_int %>%
  column_to_rownames(var = "feature") %>%
  t() %>%
  as.data.frame() %>%
  dplyr::mutate(
    dplyr::across(
      .cols = dplyr::everything(),
      # .fns = ~ as.vector(scale(.x, center = FALSE, scale = TRUE))
      # No scaling for now
      .fns = ~ .x
    )
  ) %>%
  t() %>%
  as.data.frame() %>%
  tibble::as_tibble(., rownames = "feature") %>%
  tidyr::pivot_longer(cols = dplyr::contains(".mzML")) %>%
  dplyr::left_join(
    x = .,
    y = dplyr::select(meta, sample, group, path),
    by = c("name" = "sample")
  )

p1 <- extract_feats %>%
  ggplot(
    aes(
      x = "", # feature
      y = peak_id,
      fill = sim
    )
  ) +
  geom_tile() +
  guides(x = guide_axis(angle = -45)) +
  scale_fill_gradient(
    limits = c(0, 1)
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(y = "feature")
p1

p2 <- extract_feats_int2 %>%
  # dplyr::filter(feature %in% temp_feats) %>%
  ggplot(
    aes(
      x = value,
      y = feature
    )
  ) +
  # ggplot2::stat_summary(
   #  geom = "point",
    # fun = "mean",
    # aes(color = group)
  # ) +
  # ggplot2::stat_summary(
  #   geom = "pointrange",
  #   fun.data = "mean_se",
  #   aes(color = group),
  #   position = position_jitter(height = 0.1)
  #   # width = 0.15
  # ) +
  geom_point(aes(color = group)) +
  guides(x = guide_axis(angle = -45)) +
  theme_bw() +
  theme(
    axis.title.x = element_blank()
  ) +
  labs(y = "feature")

extract_bio_sims_int <- xcms::featureValues(
  object = test2,
  method = "sum",
  value = "into",
  intensity = "into",
  filled = TRUE,
  missing = 0
) %>%
  tibble::as_tibble(., rownames = "feature") %>%
  dplyr::filter(feature %in% both_tib_feats) # extract_bio_sims$feature

extract_bio_sims_int2 <- extract_bio_sims_int %>%
  column_to_rownames(var = "feature") %>%
  t() %>%
  as.data.frame() %>%
  dplyr::mutate(
    dplyr::across(
      .cols = dplyr::everything(),
      # .fns = ~ as.vector(scale(.x, center = FALSE, scale = TRUE))
      # No scaling for now
      .fns = ~ .x
    )
  ) %>%
  t() %>%
  as.data.frame() %>%
  # feature is not unique here
  tibble::as_tibble(., rownames = "feature") %>%
  tidyr::pivot_longer(cols = dplyr::contains(".mzML")) %>%
  dplyr::left_join(
    x = .,
    y = dplyr::select(meta, sample, group, path),
    by = c("name" = "sample")
  )

p3 <- extract_bio_sims_int2 %>%
  # dplyr::filter(feature %in% temp_feats) %>%
  ggplot(
    aes(
      x = value,
      y = feature
    )
  ) +
  # ggplot2::stat_summary(
   #  geom = "point",
    # fun = "mean",
    # aes(color = group)
  # ) +
  # ggplot2::stat_summary(
  #   geom = "pointrange",
  #   fun.data = "mean_se",
  #   aes(color = group),
  #   position = position_jitter(height = 0.1)
  #   # width = 0.15
  # ) +
  geom_point(aes(color = group)) +
  guides(x = guide_axis(angle = -45)) +
  theme_bw() +
  theme(
    axis.title.x = element_blank()
  ) +
  labs(y = "feature")

(p1 + p2 +
  patchwork::plot_layout(
    guides = "collect",
    axes = "collect",
    widths = c(0.05, 0.95)
  )
)

# TODO
# Put p3 separately with a p1 that is biotransformer similarity



all_bio_sims %>%
  dplyr::filter(sim > 0.4)


################################################################################
# Heatmap prep -----------------------------------------------------------------
################################################################################
# TODO FIX THIS
yo <- anno_sims_all %>%
  dplyr::filter(
    experiment == "afzelin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon"
  ) %>%
  # Make a distribution instead
  dplyr::filter(sim > 0.3) %>%
  pull(peak_id) %>%
  unique()

pred_chrs_paths <- file.path(
  output_folders,
  "objects",
  "pred_chrs.rds"
)
pred_chrs_paths <- pred_chrs_paths[file.exists(pred_chrs_paths)]

"anno_chrs"
"pred_chrs"


# NEEDED FOR EACH PLOTTING
test <- readRDS(
  paste0(
    "/Volumes/bluecub/aglycone_release_100um_24h/output/",
    "afzelin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon",
    "/objects/anno_chrs.rds"
  )
)

test2 <- readRDS(
  paste0(
    "/Volumes/bluecub/aglycone_release_100um_24h/output/",
    "afzelin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon",
    "/objects/pred_chrs.rds"
  )
)

meta <- readr::read_csv(
  paste0(
    "/Volumes/bluecub/aglycone_release_100um_24h/output/",
    "afzelin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon",
    "/tables/metadata.csv"
  ),
  progress = FALSE,
  show_col_types = FALSE
) %>%
  dplyr::mutate(sample = basename(path)) %>%
  dplyr::relocate("sample", .before = "group")

full_limma <- readr::read_csv(
  paste0(
    "/Volumes/bluecub/aglycone_release_100um_24h/output/",
    "afzelin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon",
    "/tables/limma_norm_fill_imp.csv"
  )
)

groups_to_use <- unique(meta$group)
group_colors <- paste0(
  RColorBrewer::brewer.pal(
    n = length(groups_to_use), "Set1"
  )[seq_along(groups_to_use)]
)
group_colors <- setNames(group_colors, groups_to_use)

temp_int <- c(
  "FT22291",
  "FT12509",
  "FT11084",
  "FT10098",
  "FT11794",
  "FT07124"
)

for (i in temp_int) {
  tin <- plot_feat_chrom_int(
    feature_chrom = test,
    feature = i,
    method = "sum",
    value = "into",
    filled = TRUE,
    missing = 0,
    ms_level = 1,
    # This should be the same as the folder it came from
    save_loc = NULL,
    device = "pdf",
    feat_pairs = FALSE,
    overwrite = FALSE
  )
  print(tin$combined)
  readline("Enter for next: ")
}

tin <- plot_feat_chrom_int(
  feature_chrom = test,
  feature = "FT11084",
  method = "sum",
  value = "into",
  filled = TRUE,
  missing = 0,
  ms_level = 1,
  # This should be the same as the folder it came from
  save_loc = NULL,
  device = "pdf",
  feat_pairs = FALSE,
  overwrite = FALSE
)
print(tin$combined)


# TESTING

tester <- readRDS(
      file.path(
        dirname(dirname(bio_sim_paths[1])),
        "snakemake_objects",
        "15_prep_annotation_biotransformation.rds"
      )
    )


all_anno_sims %>%
  dplyr::filter(
    experiment == "afzelin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon"
  ) %>%
  dplyr::filter(feature %in% temp_int)

