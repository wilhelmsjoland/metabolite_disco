source("scripts/functions.R")
source("analyses/analyses_functions.R")
source("analyses/met_analysis_args.R")

library(cli)
library(tidyverse)
library(optparse)

# create folders
dir.create(file.path(opt$output), FALSE, TRUE)
dir.create(file.path(opt$output, "graphs", "features"), FALSE, TRUE)
#

snake_files <- list.files(
  file.path(
    opt$input,
    "snakemake_objects"
  ),
  full.names = TRUE
)

purrr::walk(
  .x = snake_files,
  .f = ~ {
    read_snakes(.x)
    cli::cli_alert_success("Imported {.val {.x}}")
  }
)

# ==============================================================================
# Create list of interesting comparisons ---------------------------------------
# ==============================================================================
# TODO
# Make a list out of this -> of interesting comparisons
int_upset_comp <- c(
  "bu_mutant_apiin-bu_mutant_control",
  "bu_mutant_apiin-bu_wt_control",
  "bu_mutant_control-bu_wt_apiin",
  "bu_wt_apiin-bu_wt_control"
)

int_upset_comp2 <- c(
  "bu_mutant_apiin-bu_mutant_control",
  "bu_mutant_apiin-bu_wt_control",
  "bu_mutant_control-bu_wt_apiin",
  "bu_wt_apiin-bu_wt_control",
  "bu_mutant_apiin-bu_wt_apiin"
)

int_upset_comps <- list(
  int_upset_comp,
  int_upset_comp2
)

# ==============================================================================
# Run analysis scripts ---------------------------------------------------------
# ==============================================================================

source("analyses/01_generate_comparisons.R")
# source("analysis/02_produce_chromatograms.R")
# TODO
# fix this part
# source("analysis/03_plotting_features.R")
