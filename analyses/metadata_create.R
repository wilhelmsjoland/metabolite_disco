# ==============================================================================
# Create configs and metadata files --------------------------------------------
# ==============================================================================
config_path <- "V:/aglycone_release_100um_24h/configs"
metadata_path <- "V:/aglycone_release_100um_24h/metadata"
dir.create("V:/aglycone_release_100um_24h/configs", FALSE, TRUE)
dir.create("V:/aglycone_release_100um_24h/metadata", FALSE, TRUE)

experiments %>%
  split(.$experiment) %>%
  purrr::walk2(
    .x = .,
    .y = names(.),
    .f = ~ {
      # Save the metadata CSV for this experiment
      readr::write_csv(.x, file.path(metadata_path, paste0(.y, ".csv")))

      mol_info <- glycone_pairs_metadata %>%
        dplyr::filter(
          glycoside == unique(
            .x$unclean_condition[.x$condition != "ycfa_glucose"]
          )
        )

      # # Create the config and write as YAML
      config <- list(
        output = paste0("V:/aglycone_release_100um_24h/output/", .y),
        data_path = "V:/aglycone_release_100um_24h/data/mzml_files",
        meta_file = file.path(metadata_path, paste0(.y, ".csv")),
        internal_standard = "C7H8O2",
        is_adduct = "[M-H]-",
        ppm_global = 25,
        sn_threshold = 10,
        mzdiff = 0.01,
        beta_cor_threshold = 0.8,
        beta_snr_threshold = 3,
        bw_first_grouping = 3,
        peak_anchor_sd = 2,
        min_fraction_align = 0.9,
        extra_peaks = 0,
        span = 0.6,
        bw_second_grouping = 0.5,
        missingness = 50,
        gap_filling = "norm_fill_imp",
        qvalue = 0.05,
        polarity = "negative",
        biotransf_file = paste0(
          "V:/aglycone_release_100um_24h/data/biotransformations",
          "/biotransformations.csv"
        ),
        rpairs_path = paste0(
          "C:/Users/wilhelm/Documents/MEGA/01_juniper/01_arbete/01_projekt/",
          "03_psm/scripts/search_compounds/output/rpairs.tsv"
        ),
        glycoside = mol_info$glycoside_form,
        aglycone = mol_info$aglycone_form,
        ppm_match = 5,
        all_vs_all = FALSE,
        smiles = mol_info$aglycone_SMILES,
        biot_dir = paste0(
          "C:/Users/wilhelm/Documents/MEGA/01_juniper/01_arbete/01_projekt/",
          "03_psm/biotransformer3.0jar"
        ),
        seed = 123,
        cores = 6
      )

      yaml::write_yaml(config, file.path(config_path, paste0(.y, ".yaml")))
    }
  )

length(list.files(config_path)) == max(experiments$experiment_id)
length(list.files(metadata_path)) == max(experiments$experiment_id)
