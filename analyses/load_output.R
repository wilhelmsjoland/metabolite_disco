suppressPackageStartupMessages(
  {
    library(tidyverse)
    library(ComplexHeatmap)
    library(patchwork)
    library(xcms)
    library(Spectra)
    library(rstatix)
    library(Rdisop)
    library(MetaboCoreUtils)
    library(MsCoreUtils)
    library(Spectra)
  }
)
source("scripts/functions.R")
source("analyses/analyses_functions.R")

###############################################################################
# Setup output folder ----------------------------------------------------------
################################################################################
output_path <- "/Volumes/bluecub/aglycone_release_100um_24h/output/experiment"
exp_to_use <- "afzelin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon"
output_folders <- list.files(
  output_path,
  full.names = TRUE
)

fold_change_min <- 1
# It means: the minimum of the substrate group means must be at least
# fold_change_min times larger than the maximum of the glucose group means.
# With fold_change_min <- 10: if the highest glucose group mean is 5,000,
# both substrate group means must be above 50,000 to pass.
cores <- 4
BiocParallel::register(BiocParallel::MulticoreParam(workers = cores))

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
    tryCatch(
      {
        bio_sim <- readr::read_csv(
          file = .x,
          progress = FALSE,
          show_col_types = FALSE
        )
        sig_higher <- extract_sig_higher(
          exp_dir = dirname(dirname(.x)),
          fold_change_min = fold_change_min,
          top_pct = NULL
        )
        bio_sim %>%
          dplyr::filter(feature %in% sig_higher)
      },
      error = function(e) {
        cli::cli_alert_warning(
          "Skipping {.path {basename(dirname(dirname(.x)))}}: {e$message}"
        )
        tibble::tibble()
      }
    )
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
    tryCatch(
      {
        anno_sim <- readr::read_csv(
          file = .x,
          progress = FALSE,
          show_col_types = FALSE
        )
        sig_higher <- extract_sig_higher(
          exp_dir = dirname(dirname(.x)),
          fold_change_min = fold_change_min,
          top_pct = NULL
        )
        anno_sim %>%
          dplyr::filter(peak_id %in% sig_higher)
      },
      error = function(e) {
        cli::cli_alert_warning(
          "Skipping {.path {basename(dirname(dirname(.x)))}}: {e$message}"
        )
        tibble::tibble()
      }
    )
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
# Read all significant features with intensities -------------------------------
################################################################################
int_paths <- file.path(
  output_folders,
  "tables",
  "norm_fill_imp_untransformed.csv"
)
int_paths <- int_paths[file.exists(int_paths)]

xchr9_datas <- purrr::map(
  .x = int_paths,
  .f = ~ {
    tryCatch(
      {
        intensities <- readr::read_csv(
          file = .x,
          progress = FALSE,
          show_col_types = FALSE
        )
        sig_higher <- extract_sig_higher(
          exp_dir = dirname(dirname(.x)),
          fold_change_min = 5, # 70
          top_pct = 0.05
          # fold_change_min
        )
        intensities %>%
          dplyr::filter(feature %in% sig_higher) %>%
          tidyr::pivot_longer(cols = dplyr::contains(".mzML"))
      },
      error = function(e) {
        cli::cli_alert_warning(
          "Skipping {.path {basename(dirname(dirname(.x)))}}: {e$message}"
        )
        tibble::tibble()
      }
    )
  }
)
names(xchr9_datas) <- basename(dirname(dirname(int_paths)))
xchr9_datas2 <- purrr::keep(xchr9_datas, ~ nrow(.x) > 0)

all_xchr9_data <- purrr::map_dfr(
  .x = xchr9_datas2,
  .f = ~ .x,
  .id = "experiment"
) %>%
  dplyr::relocate("feature", .after = "experiment")

all_xchr9_data %>%
  dplyr::filter(experiment == paste0(
    exp_to_use
  )
)

################################################################################
# Load experiment-specific data ------------------------------------------------
################################################################################
anno_chrs <- readRDS(
  paste0(
    output_path, "/",
    exp_to_use,
    "/objects/anno_chrs.rds"
  )
)

pred_chrs <- readRDS(
  paste0(
    output_path, "/",
    exp_to_use,
    "/objects/pred_chrs.rds"
  )
)
meta <- readr::read_csv(
  paste0(
    output_path, "/",
    exp_to_use,
    "/tables/metadata.csv"
  ),
  progress = FALSE,
  show_col_types = FALSE
) %>%
  dplyr::mutate(sample = basename(path)) %>%
  dplyr::relocate("sample", .before = "group")

full_limma <- readr::read_csv(
  paste0(
    output_path, "/",
    exp_to_use,
    "/tables/limma_norm_fill_imp.csv"
  ),
  progress = FALSE,
  show_col_types = FALSE
)

xchr9 <- readRDS(
  paste0(
    output_path, "/",
    exp_to_use,
    "/objects/xchr9.rds"
  )
)
sp <- spectra(xchr9)
sp@backend@spectraData$dataStorage <- gsub(
  r"(V:\aglycone_release_100um_24h\data\mzml_files\)",
  "/Volumes/bluecub/aglycone_release_100um_24h/data/experiment/mzml/",
  sp@backend@spectraData$dataStorage,
  fixed = TRUE
)
xchr9@spectra <- sp
xchr9@spectra <- Spectra::setBackend(
  spectra(xchr9),
  MsBackendMemory()
)

# After loading xchr9, sync sampleData with local meta
sd <- MsExperiment::sampleData(xchr9)
sd$path <- meta$path
rownames(sd) <- meta$sample
sd$spectraOrigin <- gsub(
  r"(V:\aglycone_release_100um_24h\data\mzml_files\)",
  "/Volumes/bluecub/aglycone_release_100um_24h/data/experiment/mzml/",
  sd$spectraOrigin,
  fixed = TRUE
)
MsExperiment::sampleData(xchr9) <- sd
meta$path <- gsub(
  "V:/aglycone_release_100um_24h",
  "/Volumes/bluecub/aglycone_release_100um_24h",
  meta$path,
  fixed = TRUE
)

################################################################################
# Extract only interesting significant differences where G + Samp > Samp -------
################################################################################
# 1. define interesting comparisons
# 2. write out logic e.g. all sig in certain comps & p < 0.01
# 3. ascertain that i can extract the interesting ones by looking at chroms
# 4. plot heatmap somehow

sim_filter <- 0.2

extract_feats <- all_anno_sims %>%
  dplyr::filter(
    experiment == paste0(
      exp_to_use
    )
  ) %>%
  dplyr::group_by(adduct) %>%
  dplyr::distinct(target_inchikey, .keep_all = TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::filter(sim > sim_filter)

extract_bio_sims <- all_bio_sims %>%
  dplyr::filter(
    experiment == paste0(
      exp_to_use
    )
  ) %>%
  dplyr::group_by(adduct) %>%
  dplyr::distinct(InChIKey, .keep_all = TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::filter(sim > sim_filter)

xchr9_all_ints <- all_xchr9_data %>%
  dplyr::filter(experiment == paste0(
    exp_to_use
  )
) %>%
  dplyr::left_join(
    x = .,
    y = dplyr::select(meta, sample, group, path),
    by = c("name" = "sample")
  ) %>%
  dplyr::relocate(group, .after = "feature") %>%
  dplyr::mutate(
    title = paste0(
      feature, " ",
      round(mzmed, 2), " ",
      round(rtmed, 2)
    )
  )

shared_features <- unique(
  c(
    extract_feats$peak_id,
    extract_bio_sims$feature,
    xchr9_all_ints$feature
  )
)

feat_map <- xcms::featureDefinitions(xchr9) %>%
  tibble::as_tibble(., rownames = "feature") %>%
  dplyr::select(feature, mzmed, rtmed) %>%
  dplyr::arrange(mzmed) %>%
  dplyr::mutate(
    title = paste0(
      feature, "_",
      round(mzmed, 2), "_",
      round(rtmed, 2)
    )
  ) %>%
  dplyr::filter(feature %in% shared_features) %>%
  dplyr::pull(title, name = feature)
# names = feature IDs, values = display titles, sorted by mzmed

feature_levels <- names(feat_map)

extract_feats_int <- xcms::featureValues(
  object = anno_chrs,
  method = "sum",
  value = "into",
  intensity = "into",
  filled = TRUE,
  missing = 0
) %>%
  tibble::as_tibble(., rownames = "feature") %>%
  dplyr::filter(feature %in% shared_features) # extract_feats$peak_id

extract_feats_long <- extract_feats_int %>%
  tibble::column_to_rownames(var = "feature") %>%
  t() %>%
  as.data.frame() %>%
  dplyr::mutate(
    dplyr::across(
      .cols = dplyr::everything(),
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
  object = pred_chrs,
  method = "sum",
  value = "into",
  intensity = "into",
  filled = TRUE,
  missing = 0
) %>%
  tibble::as_tibble(., rownames = "feature") %>%
  dplyr::filter(feature %in% shared_features)

extract_bio_sims_long <- extract_bio_sims_int %>%
  tibble::column_to_rownames(var = "feature") %>%
  t() %>%
  as.data.frame() %>%
  dplyr::mutate(
    dplyr::across(
      .cols = dplyr::everything(),
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

################################################################################
# Plotting of the interesting features -----------------------------------------
################################################################################
p1 <- extract_feats %>%
  dplyr::group_by(peak_id) %>%
  dplyr::mutate(hit = dplyr::row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(peak_id = factor(peak_id, levels = feature_levels)) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = hit,
      y = peak_id,
      fill = sim
    )
  ) +
  ggplot2::geom_tile(color = "black") +
  ggplot2::scale_y_discrete(drop = FALSE, labels = feat_map) +
  ggplot2::scale_fill_gradientn(
    colours = c("#e4f1e1", "#4a9a8e", "#023c3f"),
    values = scales::rescale(c(0.2, 0.4, 1)),
    limits = c(0.2, 1)
  ) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(c(0, 0))) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.ticks.x = ggplot2::element_blank()) +
  ggplot2::labs(
    x = "n predictions",
    y = "feature",
    # title = stringr::str_wrap("Massbank annotated features", 20),
    title = "Massbank annotated features"
  )

p2 <- extract_feats_long %>%
  dplyr::filter(feature %in% extract_feats$peak_id) %>%
  dplyr::mutate(feature = factor(feature, levels = feature_levels)) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = value,
      y = feature
    )
  ) +
  ggplot2::geom_point(aes(color = group)) +
  ggplot2::scale_y_discrete(drop = FALSE, labels = feat_map) +
  ggplot2::scale_x_continuous(
    # transform = "sqrt"
    transform = scales::pseudo_log_trans(sigma = 1e5)
    # transform = scales::transform_boxcox(p = 0.4)
  ) +
  ggplot2::guides(x = ggplot2::guide_axis(angle = -45)) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank()
  ) +
  ggplot2::labs(y = "feature")

p3 <- extract_bio_sims %>%
  dplyr::group_by(feature) %>%
  dplyr::mutate(hit = dplyr::row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(feature = factor(feature, levels = feature_levels)) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = hit,
      y = feature,
      fill = sim
    )
  ) +
  ggplot2::geom_tile(color = "black") +
  ggplot2::scale_y_discrete(drop = FALSE, labels = feat_map) +
  ggplot2::scale_fill_gradientn(
    colours = c("#e4f1e1", "#4a9a8e", "#023c3f"),
    values = scales::rescale(c(0.2, 0.4, 1)),
    limits = c(0.2, 1)
  ) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(c(0, 0))) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.ticks.x = ggplot2::element_blank()
  ) +
  ggplot2::labs(
    x = "n predictions",
    y = "feature",
    # title = stringr::str_wrap("Biotransformer 3.0 predictions", 20),
    title = "Biotransformer 3.0 predictions"
  )

p4 <- extract_bio_sims_long %>%
  dplyr::filter(feature %in% extract_bio_sims$feature) %>%
  dplyr::mutate(feature = factor(feature, levels = feature_levels)) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = value,
      y = feature # feature
    )
  ) +
  ggplot2::geom_point(aes(color = group)) +
  ggplot2::scale_y_discrete(drop = FALSE, labels = feat_map) +
  ggplot2::scale_x_continuous(
    # transform = "sqrt"
    transform = scales::pseudo_log_trans(sigma = 1e5)
    # transform = scales::transform_boxcox(p = 0.4)
  ) +
  ggplot2::guides(x = ggplot2::guide_axis(angle = -45)) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank()
  ) +
  ggplot2::labs(y = "feature")

p5 <- xchr9_all_ints %>%
  dplyr::mutate(feature = factor(feature, levels = feature_levels)) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = value,
      y = feature, # feature
      color = group
    )
  ) +
  ggplot2::geom_point() +
  ggplot2::theme_bw() +
  ggplot2::scale_y_discrete(
    drop = FALSE,
    labels = feat_map
  ) +
  ggplot2::scale_x_continuous(
    # transform = "sqrt"
    transform = scales::pseudo_log_trans(sigma = 1e5)
    # transform = scales::transform_boxcox(p = 0.4)
  ) +
  ggplot2::guides(x = ggplot2::guide_axis(angle = -45)) +
  ggplot2::labs(
    title = "Largest delta peak area"
  )

final_p <- p1 + p2 + p3 + p4 + p5 +
  patchwork::plot_layout(
    guides = "collect",
    axes = "collect",
    axis_titles = "collect",
    widths = c(0.1, 0.3, 0.1, 0.3, 0.3)
  ) &
  ggplot2::theme(
    axis.title.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(size = 7)
  )
final_p
################################################################################
# Extracting chromatograms -----------------------------------------------------
################################################################################
xchr9_int_chrs <- xcms::featureChromatograms(
  BPPARAM = bpparam(),
  chunkSize = cores,
  object = xchr9,
  expandRt = 0,
  expandMz = 0,
  aggregationFun = "sum",
  filled = TRUE,
  features = unique(xchr9_all_ints$feature),
  missing = 0,
  msLevel = 1L,
  return.type = "XChromatograms"
)
colnames(xchr9_int_chrs) <- sub(".*[/\\\\]", "", colnames(xchr9_int_chrs))

for (i in unique(xchr9_all_ints$feature)) {
  xchr9_int_chr_result <- plot_feature(
    feature_chrom = xchr9_int_chrs,
    feature = i,
    meta = meta,
    limma_results = full_limma,
    method = "sum",
    value = "into",
    filled = TRUE,
    missing = 0,
    ms_level = 1L,
    # This should be the same as the folder it came from
    save_loc = NULL,
    device = "pdf",
    overwrite = FALSE
  )
  print(xchr9_int_chr_result$full)
  # readline("Enter for next: ")
}

plot_feature(
  feature_chrom = xchr9_int_chrs,
  feature = "FT18047",
  meta = meta,
  limma_results = full_limma,
  method = "sum",
  value = "into",
  filled = TRUE,
  missing = 0,
  ms_level = 1L,
  # This should be the same as the folder it came from
  save_loc = NULL,
  device = "pdf",
  overwrite = FALSE
)$full

# TODO
# export mgf for sirius
# Run through all the standards
# and select the peaks that match the aglycone in the end
# Find all the standards that have been used
# and find their retention time in addition to their mass
# Import all gly_agly & match based on mass?
# Use an excel table for this?


all_anno_sims %>%
  dplyr::filter(
    experiment == paste0(
      exp_to_use
    )
  ) %>%
  dplyr::filter(feature %in% unique(xchr9_all_ints$feature))

all_bio_sims %>%
  dplyr::filter(
    experiment == paste0(
      exp_to_use
    )
  ) %>%
  dplyr::filter(feature %in% unique(xchr9_all_ints$feature))

# subset_matched_diffs <- readRDS(
#   paste0(
#     output_path, "/",
#     "/afzelin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon",
#     "/objects/subset_matched_diffs.rds"
#   )
# )

