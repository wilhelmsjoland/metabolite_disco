library(tidyverse)
library(arrow)

exp_dir <- "/Volumes/bluecub/aglycone_release_100um_24h/output"
res_out <- "/Volumes/bluecub/aglycone_release_100um_24h/combined_results"

exp_names <- list.files(
  exp_dir,
  recursive = FALSE,
  full.names = FALSE
)

all_hits <- exp_names %>%
  purrr::map(
    .x = .,
    .f = ~ {
      bio <- file.path(
        exp_dir,
        .x,
        "report",
        "biotransformer_similarities.parquet"
      )

      anno <- file.path(
        exp_dir,
        .x,
        "report",
        "annotation_similarities.parquet"
      )

      anno_wrangle <- arrow::read_parquet(file = anno) %>%
        dplyr::mutate(
          experiment = .x,
          source = "anno",
          sim_mol = NA
        ) %>%
        dplyr::select(
          "experiment",
          "source",
          "sim_mol",
          "feature" = peak_id,
          "feature_pred" = feature,
          "adduct",
          "mass" = "exact_mass",
          "mz",
          "rtime",
          "sim",
          "formula" = "target_formula",
          "smiles" = "target_smiles",
          "inchikey" = "target_inchikey",
          "ppm_error"
        )

      bio_wrangle <- arrow::read_parquet(file = bio) %>%
        dplyr::mutate(
          experiment = .x,
          source = "bio",
          feature_pred = paste0(feature, ".", met_id),
          ppm_error = abs((mzmed - mz) / mz * 1e6)
        ) %>%
        dplyr::select(
          "experiment",
          "source",
          "sim_mol",
          "feature",
          "feature_pred" = feature_pred,
          "adduct",
          "mass",
          "mz" = "mzmed",
          "rtime" = "rtmed",
          "sim",
          "formula" = `Molecular formula`,
          "smiles" = "SMILES",
          "inchikey",
          "ppm_error"
        )

      dplyr::bind_rows(
        ... = list(
          anno_wrangle,
          bio_wrangle
        )
      )
    }
  ) %>%
  dplyr::bind_rows()

arrow::write_parquet(
  x = all_hits,
  sink = file.path(
    res_out,
    "aglycone_release_100um_24h_hits.parquet"
  )
)
