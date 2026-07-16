suppressPackageStartupMessages(
  {
    library(tidyverse)
    library(xcms)
    library(Rdisop)
    library(MetaboCoreUtils)
    library(MsCoreUtils)
    library(Spectra)
    library(MsBackendMgf)
    library(BiocParallel)
  }
)
source("scripts/functions.R")
source("analyses/analyses_functions.R")

# TODO
# - Move extract_sig_higher to the functions file instead of analyses_functions
# - Integrate this into the pipeline

###############################################################################
# Setup output folder ----------------------------------------------------------
################################################################################
input_path <- file.path(
  "/Volumes/bluecub/aglycone_release_100um_24h/output",
  "afzelin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon"
)
mgf_export_path <- file.path(input_path, "sirius")
sim_filter <- 0.2
fold_change_min <- 1
# It means: the minimum of the substrate group means must be at least
# fold_change_min times larger than the maximum of the glucose group means.
# With fold_change_min <- 10: if the highest glucose group mean is 5,000,
# both substrate group means must be above 50,000 to pass.

################################################################################
# Read bio_sim -----------------------------------------------------------------
################################################################################
bio_sim_path <- file.path(
  file.path(input_path),
  "tables",
  "biotransformer_similarities.csv"
)
bio_sim_path <- bio_sim_path[file.exists(bio_sim_path)]

bio_sim <- tryCatch(
  {
    bio_sim <- readr::read_csv(
      file = bio_sim_path,
      progress = FALSE,
      show_col_types = FALSE
    )
    sig_higher <- extract_sig_higher(
      exp_dir = dirname(dirname(bio_sim_path)),
      fold_change_min = fold_change_min,
      top_pct = NULL
    )
    bio_sim %>%
      dplyr::filter(feature %in% sig_higher)
  },
  error = function(e) {
    cli::cli_alert_warning(
      "Skipping {.path {basename(dirname(dirname(bio_sim_path)))}}: {e$message}"
    )
    tibble::tibble()
  }
)

if (nrow(bio_sim) == 0) {
  cli::cli_alert_warning(
    "No significant features found in {.path {bio_sim_path}}"
  )
}

################################################################################
# Read anno_sim ----------------------------------------------------------------
################################################################################
anno_sim_path <- file.path(
  file.path(input_path),
  "tables",
  "anno_similarities.csv"
)
anno_sim_path <- anno_sim_path[file.exists(anno_sim_path)]

anno_sim <- tryCatch(
  {
    anno_sim <- readr::read_csv(
      file = anno_sim_path,
      progress = FALSE,
      show_col_types = FALSE
    )
    sig_higher <- extract_sig_higher(
      exp_dir = dirname(dirname(anno_sim_path)),
      fold_change_min = fold_change_min,
      top_pct = NULL
    )
    anno_sim %>%
      dplyr::filter(peak_id %in% sig_higher)
  },
  error = function(e) {
    cli::cli_alert_warning(
      paste0(
        "Skipping",
        " {.path {basename(dirname(dirname(anno_sim_path)))}}: {e$message}"
      )

    )
    tibble::tibble()
  }
)

if (nrow(anno_sim) == 0) {
  cli::cli_alert_warning(
    "No significant features found in {.path {bio_sim_path}}"
  )
}

################################################################################
# Read significant features with intensities -----------------------------------
################################################################################
int_path <- file.path(
  file.path(input_path),
  "tables",
  "norm_fill_imp_untransformed.csv"
)
int_path <- int_path[file.exists(int_path)]

xchr9_data <- tryCatch(
  {
    intensities <- readr::read_csv(
      file = int_path,
      progress = FALSE,
      show_col_types = FALSE
    )
    sig_higher <- extract_sig_higher(
      exp_dir = dirname(dirname(int_path)),
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
      "Skipping {.path {basename(dirname(dirname(int_path)))}}: {e$message}"
    )
    tibble::tibble()
  }
)

if (nrow(bio_sim) == 0) {
  cli::cli_alert_warning(
    "No significant features found in {.path {int_path}}"
  )
}

################################################################################
# Import the rest --------------------------------------------------------------
################################################################################
anno_chrs <- readRDS(
  file = file.path(input_path, "objects", "anno_chrs.rds")
)

pred_chrs <- readRDS(
  file = file.path(input_path, "objects", "pred_chrs.rds")
)

meta <- readr::read_csv(
  file = file.path(input_path, "tables", "metadata.csv"),
  progress = FALSE,
  show_col_types = FALSE
) %>%
  dplyr::mutate(sample = basename(path)) %>%
  dplyr::relocate("sample", .before = "group")

full_limma <- readr::read_csv(
  file = file.path(input_path, "tables", "limma_norm_fill_imp.csv"),
  progress = FALSE,
  show_col_types = FALSE
)

xchr9 <- readRDS(
  file = file.path(input_path, "objects", "xchr9.rds")
)

################################################################################
# Extract only interesting significant differences where G + Samp > Samp -------
################################################################################
# 1. define interesting comparisons
# 2. write out logic e.g. all sig in certain comps & p < 0.01
# 3. ascertain that i can extract the interesting ones by looking at chroms
# 4. plot heatmap somehow

extract_feats <- anno_sim %>%
  dplyr::group_by(adduct) %>%
  dplyr::distinct(target_inchikey, .keep_all = TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::filter(sim > sim_filter)

extract_bio_sims <- bio_sim %>%
  dplyr::group_by(adduct) %>%
  dplyr::distinct(InChIKey, .keep_all = TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::filter(sim > sim_filter)

xchr9_all_ints <- xchr9_data %>%
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

################################################################################
# Export mgfs ------------------------------------------------------------------
################################################################################

# Extract MS1 spectra at feature apexes
feat_spectra <- xcms::featureSpectra(
  xchr9,
  msLevel = 1L,
  method = "closest_rt", # Makes sense for MS1 but not for MS2
  BPPARAM = BiocParallel::SerialParam(),
  features = feature_levels
)

feat_spectra$precursorCharge <- -1L

feat_export <- Spectra::selectSpectraVariables(
  feat_spectra,
  c(
    "scanIndex",
    "mz",
    "intensity",
    "rtime",
    "precursorMz",
    "precursorCharge",
    "msLevel"
  )
)

map <- c(
  feature_id = "TITLE",
  Spectra::spectraVariableMapping(MsBackendMgf::MsBackendMgf())
)

dir.create(
  mgf_export_path,
  recursive = TRUE,
  showWarnings = FALSE
)

# Export into the output folder
Spectra::export(
  object = feat_export,
  BPPARAM = BiocParallel::SerialParam(),
  backend = MsBackendMgf::MsBackendMgf(),
  mapping = map,
  file = file.path(mgf_export_path, "ms1_features.mgf")
)
