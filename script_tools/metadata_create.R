# Run
source("script_tools/metadata_prep.R")
# Need to set to TRUE for this to run
run_create_metadata <- TRUE

# ==============================================================================
# Create configs and metadata files --------------------------------------------
# ==============================================================================
# "input/experiment/configs"
config_path <- "input/aglycone_release_100um_24h/config"
# "input/experiment/metadata"
metadata_path <- "input/aglycone_release_100um_24h/metadata"
dir.create(config_path, FALSE, TRUE)
dir.create(metadata_path, FALSE, TRUE)

# Add bad samples to this list
bad_samples <- readxl::read_xlsx(
  path = file.path(exp_path, "data", "experiment", "info", "bad_samples.xlsx"),
  sheet = "bad_samples"
)

# Bad samples
# experiments %>%
#   dplyr::filter(sample %in% bad_samples$sample)

if (isTRUE(run_create_metadata)) {
  cli::cli_alert_info("Creating new configs and metadata")
  experiments %>%
    dplyr::filter(!sample %in% bad_samples$sample) %>%
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
              .x$condition[.x$condition != "ycfa_glucose"]
            )
          )

        # # Create the config and write as YAML
        config <- list(
          output = paste0(
            # "/Volumes/bluecub/aglycone_release_100um_24h/output/",
            "/cfs/klemming/projects/supr/sjoland_naiss/project/",
            "metabolite_disco/output/",
            .y
          ),
          data_path = paste0(
            # "/Volumes/bluecub/aglycone_release_100um_24h/data/experiment/mzml"
            "/cfs/klemming/projects/supr/sjoland_naiss",
            "/data/aglycone_release_100um_24h/mzml"
          ),
          meta_file = file.path(metadata_path, paste0(.y, ".csv")),
          internal_standard = "C7H8O2",
          is_adduct = "[M-H]-",
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
          qvalue = 0.05,
          polarity = "negative",
          mass_shift_path = "mass_shift/mass_shifts.csv",
          rpairs_path = "mass_shift/rpairs.tsv",
          metabolite_search = paste0(
            mol_info$glycoside_form,
            ";",
            mol_info$aglycone_form
          ),
          ppm_match = 5,
          all_vs_all = FALSE,
          smiles = paste0(
            mol_info$aglycone, ",",
            mol_info$aglycone_SMILES,
            ";",
            mol_info$glycoside, ",",
            mol_info$glycoside_SMILES
          ),
          annotate_path = paste0(
            "/cfs/klemming/projects/supr/sjoland_naiss/database/pubchem/",
            "260730_pubchem.db"
          ),
          seed = 123,
          cores = 6
        )

        yaml::write_yaml(config, file.path(config_path, paste0(.y, ".yaml")))
      }
    )
}

length(list.files(config_path)) == max(experiments$experiment_id)
length(list.files(metadata_path)) == max(experiments$experiment_id)