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

for (script in pipeline_scripts[18]) {
  rule <- gsub("[0-9]+_", "", basename(script))
  rule <- tools::file_path_sans_ext(rule)
  snakemake <- mock_snakemake(
    rule = rule,
    config_file = paste0(
      "/Users/wilhelm/Proton/01_juniper/01_arbete/01_projekt/03_psm/input/",
      "experiment/",
      "afzelin_test.yaml"

    )
  )
  source(script)
}
