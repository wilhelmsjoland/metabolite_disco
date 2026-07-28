# ==============================================================================
# Mock snakemake object for interactive R sessions
# ==============================================================================
source("analyses/analyses_functions.R")
library(tools)

pipeline_scripts <- list.files(
  path = "scripts",
  pattern = "^[0-9]+_.*\\.R$",
  full.names = TRUE
)

pipeline_scripts[19]

for (script in pipeline_scripts) {
  rule <- gsub("[0-9]+_", "", basename(script))
  rule <- tools::file_path_sans_ext(rule)
  snakemake <- mock_snakemake(
    rule = rule,
    config_file = paste0(
      # "V:/aglycone_release_100um_24h/configs/",
      # "avicularin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon.yaml"
      # "/Users/wilhelm/Documents/from_ssd/aglycone_release_100um_24h/configs/",
      # "afzelin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon.yaml"
      # "example_config.yaml"
      "/Users/wilhelm/Proton/01_juniper/01_arbete/01_projekt/03_psm/input/",
      "experiment/",
      "afzelin_test.yaml"

    )
  )
  source(script)
}