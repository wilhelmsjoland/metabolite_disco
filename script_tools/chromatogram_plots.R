suppressPackageStartupMessages(
  {
    library(readr)
    library(dplyr)
    library(xcms)
    library(RSQLite)
    library(svglite)
    library(arrow)
    library(optparse)
  }
)
source("scripts/functions.R")

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
    c("-f", "--features"),
    type = "character",
    default = NULL,
    help = "Features to produce chromatograms from",
    metavar = "character"
  ),
  optparse::make_option(
    c("-o", "--output"),
    type = "character",
    default = NULL,
    help = "Path to output [default: <input>/report/features.parquet]",
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

input_path <- opt$input
features <- readr::read_csv(
  file = opt$features,
  show_col_types = FALSE,
  progress = FALSE
)
output_file <- if (is.null(opt$output)) {
  file.path(input_path, "report", "features.parquet")
} else {
  opt$output
}
xchr9 <- readRDS(file.path(input_path, "objects", "xchr9.rds"))
setup <- readRDS(file.path(input_path, "snakemake_objects", "01_setup.rds"))
meta <- tibble::as_tibble(setup$meta, rownames = "sample")
group_colors <- setup$group_colors
limma <- readRDS(file.path(input_path, "snakemake_objects", "10_limma.rds"))
full_limma <- limma$full_limma

output_path_msg <- file.path(
  basename(dirname(output_file)),
  basename(output_file)
)

################################################################################
# Produce chromatograms and .svgs ----------------------------------------------
################################################################################
cli::cli_progress_step(
  msg = c(
    "Producing chromatogram plots for: ",
    "{.path {basename(basename(input_path))}}"
  ),
  msg_done = paste0(
    "Saved plots to: {.path {output_file}}"
  )
)

feature_levels <- dplyr::pull(features)

feature_chrs <- xcms::featureChromatograms(
  object = xchr9,
  expandRt = 0,
  expandMz = 0,
  aggregationFun = "sum",
  filled = TRUE,
  features = feature_levels,
  missing = 0,
  return.type = "XChromatograms"
)

svgs <- feature_levels %>%
  purrr::map(
    .x = .,
    .f = ~ {
      tmp_p <- plot_feat_chrom_int(
        feature_chrom = feature_chrs,
        feature = .x,
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

      cli::cli_alert(paste0("Saving feature: ", .x))
      svg_string <- svglite::svgstring(standalone = FALSE)
      print(tmp_p$combined)
      dev.off()

      svg_stored <- setNames(.x, svg_string())
      svg_tibble <- tibble::tibble(
        "feature" = svg_stored,
        "chromatogram" = names(svg_stored)
      )
    }
  ) %>%
  dplyr::bind_rows()

dir.create(
  path = file.path(input_path, "report"),
  recursive = TRUE,
  showWarnings = FALSE
)

arrow::write_parquet(
  x = svgs,
  sink = output_file
)

cli::cli_progress_done()