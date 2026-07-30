library(tidyverse)
library(RSQLite)
library(xcms)
source("scripts/functions.R")

exp_path <- paste0(
  "/Volumes/bluecub/aglycone_release_100um_24h/output/",
  "afzelin_b_thetaiotaomicron_vpi_5482"
)

subset_diffs <- readRDS(
  file.path(
    exp_path,
    "objects",
    "subset_matched_diffs.rds"
  )
)

test <- readRDS(file.path(exp_path, "snakemake_objects", "10_limma.rds"))
xchr9 <- readRDS(file.path(exp_path, "objects", "xchr9.rds"))
setup <- readRDS(file.path(exp_path, "snakemake_objects", "01_setup.rds"))
meta <- tibble::as_tibble(setup$meta, rownames = "sample")
group_colors <- setup$group_colors
limma <- readRDS(file.path(exp_path, "snakemake_objects", "10_limma.rds"))
full_limma <- limma$full_limma
pot_glycosides <- readr::read_csv(
  file.path(exp_path, "tables", "gly_agly.csv"),
  progress = FALSE,
  show_col_types = FALSE
) %>%
  dplyr::pull(feature) %>%
  unique()

# filter to good peaks
samp_data <- MsExperiment::sampleData(xchr9) %>%
  tibble::as_tibble(rownames = "sample") %>%
  dplyr::mutate(
    samp_idx = dplyr::row_number(),
    .before = "sample"
  ) %>%
  dplyr::select(samp_idx, "sample_name" = "sample", group)

xchr9_peaks <- xcms::chromPeaks(xchr9) %>%
  tibble::as_tibble(rownames = "peak") %>%
  dplyr::left_join(
    x = .,
    y = samp_data,
    by = c("sample" = "samp_idx")
  )

xchr9_feats <- xcms::featureDefinitions(xchr9) %>%
  tibble::as_tibble(rownames = "feature")

xchr9_feat_filtered <- xchr9_feats %>%
  dplyr::filter(
    purrr::map_lgl(
      .x = peakidx,
      .f = ~ {
        pk <- xchr9_peaks[.x, ] %>%
          dplyr::filter(
            group %in% c(
              "b_thetaiotaomicron_vpi_5482_afzelin"
            )
          )

        # For now it's all that have to be above the cutoff
        nrow(pk) > 0 && all(pk$beta_cor > 0.9)
      }
    )
  )

# do the correlations
samp_ints <- test$intensities$norm_fill_imp$log2_scale %>%
  dplyr::select(feature, contains(".mzML")) %>%
  tidyr::pivot_longer(cols = contains(".mzML")) %>%
  tidyr::pivot_wider(
    names_from = "feature",
    values_from = "value"
  ) %>%
  dplyr::rename("sample" = "name") %>%
  # Filtering to only corrs in the afzelin groups
  dplyr::left_join(
    x = .,
    y = meta,
    by = c("sample" = "sample")
  ) %>%
  dplyr::filter(stringr::str_detect(group, "afzelin")) %>%
  # Center each feature within its own group before correlating, so a
  # control-vs-treatment mean shift shared by two features doesn't produce
  # a spurious correlation between them
  dplyr::group_by(group) %>%
  dplyr::mutate(
    dplyr::across(dplyr::where(is.numeric), ~ .x - mean(.x, na.rm = TRUE))
  ) %>%
  dplyr::ungroup()

tin <- subset_diffs %>%
  dplyr::filter(
    dplyr::if_any(
      .cols = c("feat1", "feat2"),
      .fn = ~ . %in% xchr9_feat_filtered$feature
    )
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    cor = cor(samp_ints[[feat1]], samp_ints[[feat2]])
  ) %>%
  dplyr::ungroup()


raw <- readr::read_tsv("mass_shift/rpairs.tsv")

tin %>%
  dplyr::filter(
    dplyr::if_all(
      .cols = c("adduct1", "adduct2"),
      .fns = ~ .x == "[M-H]-"
    )
  ) %>%
  dplyr::filter(feat1 %in% c("FT07647")) %>%
  dplyr::filter(
    ppm_mz < 10 & ppm_diff < 15, cor > 0.8
  ) %>%
  dplyr::arrange(desc(abs(cor)), ppm_diff)

raw %>%
  slice_head(n = 4245) %>% tail(1)

tester <- tin %>%
  dplyr::filter(
    dplyr::if_all(
      .cols = c("adduct1", "adduct2"),
      .fns = ~ .x == "[M-H]-"
    )
  ) %>%
  # dplyr::filter(feat1 %in% c("FT07647")) %>%
  dplyr::filter(ppm_mz < 10) %>%
  dplyr::arrange(desc(abs(cor)), ppm_diff) %>%
  dplyr::slice_head(n = 1)

tester


# plot these pairs
feats_use <- c(
  tester$feat1,
  tester$feat2
  # pot_glycos
)

feats_chrs <- xcms::featureChromatograms(
  object = xchr9,
  expandRt = 0,
  expandMz = 0,
  aggregationFun = "sum",
  filled = TRUE,
  features = feats_use,
  missing = 0,
  return.type = "XChromatograms"
)

feat_p <- plot_feat_chrom_int(
  feature_chrom = feats_chrs,
  feature = feats_use[1],
  method = "sum",
  value = "into",
  filled = TRUE,
  missing = 0,
  ms_level = 1,
  save_loc = NULL,
  device = NULL,
  feat_pairs = FALSE,
  overwrite = FALSE
)

feat_p$combined

# TODO 
# Fix those this is a proper functions or script
# Rank hits by if they have more mass/larger than the glycoside/aglycone