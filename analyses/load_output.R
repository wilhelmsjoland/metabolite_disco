source("scripts/functions.R")
suppressPackageStartupMessages(
  {
    library(tidyverse)
    library(ComplexHeatmap)
    library(patchwork)
  }
)


###############################################################################
# Setup output folder ----------------------------------------------------------
################################################################################
output_folders <- list.files(
  "/Volumes/bluecub/aglycone_release_100um_24h/output",
  full.names = TRUE
)

# TODO
# TODO
# Mabye this should be dependent on distribution or different betweem
# the bio_sims and the anno_sims??
fold_change_min <- 2
# It means: the minimum of the substrate group means must be at least
# fold_change_min times larger than the maximum of the glucose group means.
# With fold_change_min <- 10: if the highest glucose group mean is 5,000,
# both substrate group means must be above 50,000 to pass.
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

    exp_dir <- dirname(dirname(.x))
    upset_int <- readRDS(
      file.path(exp_dir, "snakemake_objects", "13_upset.rds")
    )[["upset_intersect"]]

    set_names <- ComplexHeatmap::set_name(upset_int)
    # Contrasts where a substrate group is compared to a glucose group
    # Split contrast on "-" and check that exactly one side has ycfa_glucose
    afz_mask <- purrr::map_lgl(set_names, ~ {
      sides <- strsplit(.x, "-")[[1]]
      sum(grepl("ycfa_glucose", sides)) == 1
    })

    combs <- ComplexHeatmap::comb_name(upset_int, readable = FALSE)
    keep <- purrr::keep(combs, ~ {
      bits <- as.integer(strsplit(.x, "")[[1]])
      all(bits[afz_mask] == 1)
    })
    sig_diffs <- purrr::map(keep, ~ {
      ComplexHeatmap::extract_comb(upset_int, .x)
    }) %>%
      unlist(use.names = FALSE) %>%
      unique()

    # Filter: both substrate groups must have higher peak area than both glucose
    exp_meta <- readr::read_csv(
      file.path(exp_dir, "tables", "metadata.csv"),
      progress = FALSE, show_col_types = FALSE
    ) %>%
      dplyr::mutate(sample = basename(path))

    int_tib <- readr::read_csv(
      file.path(exp_dir, "tables", "norm_fill_imp_untransformed.csv"),
      progress = FALSE, show_col_types = FALSE
    ) %>%
      dplyr::filter(feature %in% sig_diffs) %>%
      dplyr::select(feature, dplyr::matches("\\.mzML$")) %>%
      tidyr::pivot_longer(-feature, names_to = "sample", values_to = "intensity") %>%
      dplyr::left_join(
        dplyr::select(exp_meta, sample, group),
        by = "sample"
      ) %>%
      dplyr::group_by(feature, group) %>%
      dplyr::summarise(mean_int = mean(intensity, na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = group, values_from = mean_int)

    glc_cols <- colnames(int_tib)[grepl("ycfa_glucose", colnames(int_tib))]
    sub_cols <- setdiff(colnames(int_tib), c("feature", glc_cols))

    sig_higher <- int_tib %>%
      dplyr::filter(
        pmin(!!!rlang::syms(sub_cols)) > fold_change_min * pmax(!!!rlang::syms(glc_cols))
      ) %>%
      dplyr::pull(feature)

    bio_sim %>%
      dplyr::filter(feature %in% sig_higher)
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

    exp_dir <- dirname(dirname(.x))
    upset_int <- readRDS(
      file.path(exp_dir, "snakemake_objects", "13_upset.rds")
    )[["upset_intersect"]]

    set_names <- ComplexHeatmap::set_name(upset_int)
    # Contrasts where a substrate group is compared to a glucose group
    # Split contrast on "-" and check that exactly one side has ycfa_glucose
    afz_mask <- purrr::map_lgl(set_names, ~ {
      sides <- strsplit(.x, "-")[[1]]
      sum(grepl("ycfa_glucose", sides)) == 1
    })

    combs <- ComplexHeatmap::comb_name(upset_int, readable = FALSE)
    keep <- purrr::keep(combs, ~ {
      bits <- as.integer(strsplit(.x, "")[[1]])
      all(bits[afz_mask] == 1)
    })
    sig_diffs <- purrr::map(keep, ~ {
      ComplexHeatmap::extract_comb(upset_int, .x)
    }) %>%
      unlist(use.names = FALSE) %>%
      unique()

    # Filter: both substrate groups must have higher peak area than both glucose
    exp_meta <- readr::read_csv(
      file.path(exp_dir, "tables", "metadata.csv"),
      progress = FALSE, show_col_types = FALSE
    ) %>%
      dplyr::mutate(sample = basename(path))

    int_tib <- readr::read_csv(
      file.path(exp_dir, "tables", "norm_fill_imp_untransformed.csv"),
      progress = FALSE, show_col_types = FALSE
    ) %>%
      dplyr::filter(feature %in% sig_diffs) %>%
      dplyr::select(feature, dplyr::matches("\\.mzML$")) %>%
      tidyr::pivot_longer(-feature, names_to = "sample", values_to = "intensity") %>%
      dplyr::left_join(
        dplyr::select(exp_meta, sample, group),
        by = "sample"
      ) %>%
      dplyr::group_by(feature, group) %>%
      dplyr::summarise(mean_int = mean(intensity, na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = group, values_from = mean_int)

    glc_cols <- colnames(int_tib)[grepl("ycfa_glucose", colnames(int_tib))]
    sub_cols <- setdiff(colnames(int_tib), c("feature", glc_cols))

    sig_higher <- int_tib %>%
      dplyr::filter(
        pmin(!!!rlang::syms(sub_cols)) > fold_change_min * pmax(!!!rlang::syms(glc_cols))
      ) %>%
      dplyr::pull(feature)

    anno_sim %>%
      dplyr::filter(peak_id %in% sig_higher)
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
# Load experiment-specific data ------------------------------------------------
################################################################################
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

################################################################################
# Extract only interesting significant differences where G + Samp > Samp -------
################################################################################

# 1. define interesting comparisons
# 2. write out logic e.g. all sig in certain comps & p < 0.01
# 3. ascertain that i can extract the interesting ones by looking at chroms
# 4. plot heatmap somehow

sim_filter <- 0.2

extract_feats_peaks <- all_anno_sims %>%
  dplyr::filter(
    experiment == "afzelin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon"
  ) %>%
  dplyr::group_by(adduct) %>%
  dplyr::distinct(target_inchikey, .keep_all = TRUE) %>%
  dplyr::filter(sim > sim_filter) %>%
  dplyr::pull(peak_id)

extract_bio_sims_peaks <- all_bio_sims %>%
  dplyr::filter(
    experiment == "afzelin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon"
  ) %>%
  dplyr::group_by(adduct) %>%
  dplyr::distinct(InChIKey, .keep_all = TRUE) %>%
  dplyr::filter(sim > sim_filter) %>%
  dplyr::pull(feature)

both_tib_feats <- unique(c(extract_feats_peaks, extract_bio_sims_peaks))


extract_feats <- all_anno_sims %>%
  dplyr::filter(
    experiment == "afzelin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon"
  ) %>%
  dplyr::group_by(adduct) %>%
  dplyr::distinct(target_inchikey, .keep_all = TRUE) %>%
  dplyr::filter(sim > sim_filter) %>%
  dplyr::filter(peak_id %in% both_tib_feats)

extract_bio_sims <- all_bio_sims %>%
  dplyr::filter(
    experiment == "afzelin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon"
  ) %>%
  dplyr::group_by(adduct) %>%
  dplyr::distinct(InChIKey, .keep_all = TRUE) %>%
  dplyr::filter(sim > sim_filter) %>%
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
  dplyr::filter(feature %in% both_tib_feats) # extract_feats$peak_id

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

extract_bio_sims_int <- xcms::featureValues(
  object = test2,
  method = "sum",
  value = "into",
  intensity = "into",
  filled = TRUE,
  missing = 0
) %>%
  tibble::as_tibble(., rownames = "feature") %>%
  dplyr::filter(feature %in% both_tib_feats)

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

shared_features <- purrr::reduce(
  .x = list(
    unique(extract_feats$peak_id),
    unique(extract_feats_int2$feature),
    unique(extract_bio_sims$feature),
    unique(extract_bio_sims_int2$feature)
  ),
  .f = ~ union(.x, .y)
)

p1 <- extract_feats %>%
  dplyr::group_by(peak_id) %>%
  dplyr::mutate(hit = dplyr::row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(peak_id = factor(peak_id, levels = shared_features)) %>%
  ggplot(
    aes(
      x = hit,
      y = peak_id,
      fill = sim
    )
  ) +
  geom_tile(color = "black") +
  scale_y_discrete(drop = FALSE) +
  # scale_fill_viridis_c(option = "turbo", limits = c(0, 1)) +
  # ggplot2::scale_fill_gradient(
  #   low = "#e4f1e1",
  #   high = "#0d585f",
  #   limits = c(0, 1)
  # ) +
  # scale_fill_gradient(
  #   # low = "#a1d99b", high = "#00441b",
  #   # low =  "#b4d9cc", high = "#0d585f",
  #   low = "#e4f1e1", high = "#023c3f",
  #   limits = c(0.2, 1) # ,
  #   # oob = scales::squish
  # ) +
  scale_fill_gradientn(
    colours = c("#e4f1e1", "#4a9a8e", "#023c3f"),
    values = scales::rescale(c(0.2, 0.4, 1)),
    limits = c(0.2, 1)
  ) +
  scale_x_continuous(expand = expansion(c(0, 0))) +
  theme_bw() +
  theme(
    axis.ticks.x = element_blank()
  ) +
  labs(
    x = "n predictions",
    y = "feature",
    title = stringr::str_wrap("Massbank annotated features", 20)
  )

p2 <- extract_feats_int2 %>%
  dplyr::filter(feature %in% extract_feats$peak_id) %>%
  dplyr::mutate(feature = factor(feature, levels = shared_features)) %>%
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
  scale_y_discrete(drop = FALSE) +
  scale_x_continuous(
    # transform = "sqrt"
    transform = scales::pseudo_log_trans(sigma = 1e5)
    # transform = scales::transform_boxcox(p = 0.4)
  ) +
  guides(x = guide_axis(angle = -45)) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(y = "feature")

p3 <- extract_bio_sims %>%
  dplyr::group_by(feature) %>%
  dplyr::mutate(hit = dplyr::row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(feature = factor(feature, levels = shared_features)) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = hit,
      y = feature,
      fill = sim
    )
  ) +
  ggplot2::geom_tile(color = "black") +
  ggplot2::scale_y_discrete(drop = FALSE) +
  # scale_fill_viridis_c(option = "turbo", limits = c(0, 1)) +
  # ggplot2::scale_fill_gradient(
  #   low = "#e4f1e1",
  #   high = "#0d585f",
  #   limits = c(0, 1)
  # ) +
  # scale_fill_gradient(
  #   #low =  "#b4d9cc", high = "#0d585f",
  #   low = "#e4f1e1", high = "#023c3f",
  #   limits = c(0.2, 1)
  #   # oob = scales::squish
  # ) +
  scale_fill_gradientn(
    colours = c("#e4f1e1", "#4a9a8e", "#023c3f"),
    values = scales::rescale(c(0.2, 0.4, 1)),
    limits = c(0.2, 1)
  ) +
  scale_x_continuous(expand = expansion(c(0, 0))) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.ticks.x = ggplot2::element_blank()
  ) +
  ggplot2::labs(
    x = "n predictions",
    y = "feature",
    title = stringr::str_wrap("Biotransformer 3.0 predictions", 20)
  )

p4 <- extract_bio_sims_int2 %>%
  dplyr::filter(feature %in% extract_bio_sims$feature) %>%
  dplyr::mutate(feature = factor(feature, levels = shared_features)) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = value,
      y = feature
    )
  ) +
  ggplot2::geom_point(aes(color = group)) +
  ggplot2::scale_y_discrete(drop = FALSE) +
  scale_x_continuous(
    # transform = "sqrt"
    transform = scales::pseudo_log_trans(sigma = 1e5)
    # transform = scales::transform_boxcox(p = 0.4)
  ) +
  ggplot2::guides(x = ggplot2::guide_axis(angle = -45)) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank()
  ) +
  ggplot2::labs(
    y = "feature"
  )

p1 + p2 + p3 + p4 +
  patchwork::plot_layout(
    guides = "collect",
    axes = "collect",
    axis_titles = "collect",
    widths = c(0.15, 0.35, 0.15, 0.35)
  )


# TESTING
# temp_int <- c(
#   "FT22291",
#   "FT12509",
#   "FT11084",
#   "FT10098",
#   "FT11794",
#   "FT07124"
# )

# for (i in temp_int) {
#   tin <- plot_feat_chrom_int(
#     feature_chrom = test,
#     feature = i,
#     method = "sum",
#     value = "into",
#     filled = TRUE,
#     missing = 0,
#     ms_level = 1,
#     # This should be the same as the folder it came from
#     save_loc = NULL,
#     device = "pdf",
#     feat_pairs = FALSE,
#     overwrite = FALSE
#   )
#   print(tin$combined)
#   readline("Enter for next: ")
# }

"FT16082"
