library(tidyverse)
library(pacman)
library(SQLite)
library(xcms)
library(MsFeatures)

# Add filtering of rt < 70
# add peak-grouping to the data
# add filtering per peak-group to only the dominant peak in the data
# Also add analyze_mass_shifts.R
# Try making my own pickaxe

use_path <- "/Volumes/bluecub/aglycone_release_100um_24h/output/afzelin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon"

xchr9 <- readRDS(file.path(use_path, "objects", "xchr9.rds"))
limma <- readRDS(file.path(use_path, "snakemake_objects", "10_limma.rds"))


xchr10 <- MsFeatures::groupFeatures(
  xchr9,
  param = SimilarRtimeParam(diffRt = 10, groupFun = groupConsecutive)
)

xchr11 <- MsFeatures::groupFeatures(
  xchr10,
  param = AbundanceSimilarityParam(
    threshold = 0.9, # 0.9
    transform = log2,
  )
)

plotFeatureGroups(
  xchr11,
  pch = 21,
  lwd = 2,
  col = "#00000040",
  bg = "#00000020"
)

# xchr12 <- MsFeatures::groupFeatures(
#   xchr11,
#   param = EicSimilarityParam(
#     threshold = 0.9,
#     n = 3,
#     BPPARAM = BiocParallel::MulticoreParam(10)
#   )
# )

as_tibble(featureDefinitions(xchr12), rownames = "feature") %>%
  dplyr::group_by(feature_group) %>%
  dplyr::slice_max(order_by = npeaks, n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::filter(rtmed > 60)
