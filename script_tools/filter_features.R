suppressPackageStartupMessages(
  {
    library(dplyr)
    library(xcms)
    library(RSQLite)
    library(optparse)
  }
)
################################################################################
# Optparse arguments -----------------------------------------------------------
################################################################################
option_list <- list(
  optparse::make_option(
    c("-i", "--input"),
    type = "character",
    default = NULL,
    help = "Path to metabolite_disco output directory",
    metavar = "character"
  ),
  optparse::make_option(
    c("-g", "--grouping"),
    type = "character",
    default = NULL,
    help = "Group(s) where good peak shapes are required (samples)",
    metavar = "character"
  ),
  optparse::make_option(
    c("-s", "--similarity_filter"),
    type = "double",
    default = 0.1,
    help = "Minimum similarity to keep [default %default]",
    metavar = "double"
  ),
  optparse::make_option(
    c("-f", "--lfc"),
    type = "double",
    default = 3,
    help = "Minimum fold change difference between groups [default %default]",
    metavar = "double"
  ),
  optparse::make_option(
    c("-q", "--qval"),
    type = "double",
    default = 0.1,
    help = "Minimum adjusted p value between groups [default %default]",
    metavar = "double"
  ),
  optparse::make_option(
    c("-b", "--beta_cor"),
    type = "double",
    default = 0.6,
    help = "Minimum peak shape correlation to bell curve [default %default]",
    metavar = "double"
  ),
  optparse::make_option(
    c("-o", "--output"),
    type = "character",
    default = NULL,
    help = "Path to output [default: <input>/report/retained_features.csv]",
    metavar = "character"
  )
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

if (is.null(opt$input) && !interactive()) {
  optparse::print_help(opt_parser)
  stop("--input must be supplied", call. = TRUE)
}

################################################################################
# Setup ------------------------------------------------------------------------
################################################################################
# Only for testing interactively
# input_path <- file.path(
#   "/Volumes/bluecub/aglycone_release_100um_24h/output",
#   "afzelin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon"
# )
# output_file <- file.path(input_path, "report", "retained_features.csv")

sim_filter <- opt$similarity_filter
beta_cor <- opt$beta_cor
lfc <- opt$lfc
qval <- opt$qval
grouping <- unlist(strsplit(x = opt$grouping, split = ","))

input_path <- opt$input
output_file <- if (is.null(opt$output)) {
  file.path(input_path, "report", "retained_features.csv")
} else {
  opt$output
}

xchr9 <- readRDS(file.path(input_path, "objects", "xchr9.rds"))
setup <- readRDS(file.path(input_path, "snakemake_objects", "01_setup.rds"))
meta <- tibble::as_tibble(setup$meta, rownames = "sample")
group_colors <- setup$group_colors
full_limma <- readr::read_csv(
  file = file.path(input_path, "tables", "limma_norm_fill_imp.csv"),
  progress = FALSE,
  show_col_types = FALSE
)

################################################################################
# Wrangling --------------------------------------------------------------------
################################################################################
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
            group %in% grouping
          )

        # For now it's all that have to be above the cutoff
        nrow(pk) > 0 && all(pk$beta_cor > beta_cor)
      }
    )
  )

good_features <- unique(xchr9_feat_filtered$feature)

# next - define significant features with lfc threshold
good_limma_features <- full_limma %>%
  dplyr::filter(adj.P.Val < qval) %>%
  dplyr::filter(abs(logFC) > lfc) %>%
  dplyr::pull(feature)

# next bio sim
bio_sim_path <- file.path(
  file.path(input_path),
  "tables",
  "biotransformer_similarities.csv"
)
bio_sim_path <- bio_sim_path[file.exists(bio_sim_path)]

bio_sim <- readr::read_csv(
  file = bio_sim_path,
  progress = FALSE,
  show_col_types = FALSE
)

# next anno sim
anno_sim_path <- file.path(
  file.path(input_path),
  "tables",
  "anno_similarities.csv"
)
anno_sim_path <- anno_sim_path[file.exists(anno_sim_path)]

anno_sim <- readr::read_csv(
  file = anno_sim_path,
  progress = FALSE,
  show_col_types = FALSE
)

extract_bio_sims <- bio_sim %>%
  dplyr::group_by(adduct) %>%
  dplyr::distinct(InChIKey, .keep_all = TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::filter(feature %in% good_features) %>%
  dplyr::filter(feature %in% good_limma_features) %>%
  dplyr::filter(sim > sim_filter)

extract_feats <- anno_sim %>%
  dplyr::group_by(adduct) %>%
  dplyr::distinct(target_inchikey, .keep_all = TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::filter(peak_id %in% good_features) %>%
  dplyr::filter(peak_id %in% good_limma_features) %>%
  dplyr::filter(sim > sim_filter)

shared_features <- unique(
  c(
    extract_feats$peak_id,
    extract_bio_sims$feature # ,
    # xchr9_all_ints$feature
  )
)

retained_features <- shared_features

dir.create(
  path = file.path(input_path, "report"),
  recursive = TRUE,
  showWarnings = FALSE
)

readr::write_csv(
  x = tibble::tibble(feature = retained_features),
  file = output_file
)