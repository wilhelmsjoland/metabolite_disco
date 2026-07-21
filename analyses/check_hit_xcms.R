library(tidyverse)
library(xcms)
library(RSQLite)
library(svglite)
library(arrow)
source("scripts/functions.R")

xchr9 <- readRDS(
  file = paste0(
    "/Volumes/bluecub/aglycone_release_100um_24h/output/",
    "afzelin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon",
    "/objects/",
    "xchr9.rds"

  )
)

setup <- readRDS(
  file = paste0(
    "/Volumes/bluecub/aglycone_release_100um_24h/output/",
    "afzelin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon",
    "/snakemake_objects/",
    "01_setup.rds"
  )
)
meta <- tibble::as_tibble(setup$meta, rownames = "sample")
group_colors <- setup$group_colors

limma_data <- readRDS(
  file = paste0(
    "/Volumes/bluecub/aglycone_release_100um_24h/output/",
    "afzelin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon",
    "/snakemake_objects/",
    "10_limma.rds"
  )
)
full_limma <- limma_data$full_limma

feats_to_use <- feature_levels

feature_chrs <- xcms::featureChromatograms(
  object = xchr9,
  expandRt = 0,
  expandMz = 0,
  aggregationFun = "sum",
  filled = TRUE,
  features = feats_to_use,
  missing = 0,
  return.type = "XChromatograms"
)

svgs <- feats_to_use %>%
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

      print(paste0("Running: ", .x))
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

arrow::write_parquet(
  x = svgs,
  sink = "/Users/wilhelm/Desktop/features.parquet"
)
