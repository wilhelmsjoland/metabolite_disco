library(tidyverse)
library(readr)

meta_path <-  file.path(
  "/Volumes/bluecub/data/260805_saligenin_1week",
  "original_metadata.tsv"
)
dir.create(file.path(dirname(meta_path), "metadata"), FALSE, TRUE)
dir.create(file.path(dirname(meta_path), "config"), FALSE, TRUE)

meta <- readr::read_tsv(meta_path)

meta %>%
  dplyr::filter(ATTRIBUTE_SampleType == "Sample") %>%
  split(.$ATTRIBUTE_Tissue) %>%
  purrr::map(
    .x = .,
    .f = ~ dplyr::select(
      .x,
      "sample" = "filename",
      "group" = "ATTRIBUTE_Group"
    ) %>%
    dplyr::mutate(sample = sub("\\.mzML$", "-MS1.mzML", sample))
  ) %>%
  purrr::walk2(
    .x = .,
    .y = names(.),
    .f = ~ {
      # Save the metadata CSV for this experiment
      readr::write_csv(
        .x,
        file.path(
          dirname(meta_path),
          "metadata",
          paste0(.y, ".csv")
        )
      )

      # # Create the config and write as YAML
      config <- list(
        output = file.path(
          "/cfs/klemming/projects/supr/sjoland_naiss/data/260805_saligenin_1week/output/260805_saligenin_1week_ms1",
          .y
        ),
        data_path = file.path(
          "/cfs/klemming/projects/supr/sjoland_naiss/data/260805_saligenin_1week/2026-08-08-MS1-mzML"
        ),
        meta_file = file.path(
          "/cfs/klemming/projects/supr/sjoland_naiss/data/260805_saligenin_1week",
          "metadata",
          paste0(.y, ".csv")
        ),
        internal_standard = FALSE,
        is_adduct = FALSE,
        ppm_global = 25,
        sn_threshold = 10,
        mzdiff = 0.01,
        bw_first_grouping = 3,
        bw_second_grouping = 0.5,
        peak_anchor_sd = 2,
        min_fraction_align = 0.9,
        extra_peaks = 0,
        span = 0.6,
        missingness = 50,
        gap_filling = "norm_fill_imp",
        qvalue = 0.1,
        polarity = "negative",
        mass_shift_path = "mass_shift/mass_shifts.csv",
        rpairs_path = "mass_shift/rpairs.tsv",
        metabolite_search = "C7H8O2",
        ppm_match = 5,
        all_vs_all = FALSE,
        smiles = "saligenin,C1=CC=C(C(=C1)CO)O",
        annotate_path = file.path(
          "/cfs/klemming/projects/supr/sjoland_naiss/database/pubchem",
          "260730_pubchem.db"
        ),
        seed = 123,
        cores = 6
      )

      yaml::write_yaml(
        config,
        file.path(
          dirname(meta_path),
          "config",
          paste0(.y, ".yaml")
        )
      )
    }
  )
