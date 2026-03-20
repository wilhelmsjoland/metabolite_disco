source("scripts/functions.R")
source("analyses/analyses_functions.R")
source("analyses/met_analysis_args.R")

library(cli)
library(tidyverse)
library(optparse)

# create folders
dir.create(file.path(opt$output), FALSE, TRUE)
dir.create(file.path(opt$output, "graphs", "features"), FALSE, TRUE)

# wrangle snake files
snake_files <- list.files(
  file.path(
    opt$input,
    "snakemake_objects"
  ),
  full.names = TRUE
)
output_folder_names <- basename(dirname(dirname(snake_files)))
output_folder_file <- basename(snake_files)
output_names <- paste0(output_folder_names, "_", output_folder_file)
snake_files <- setNames(snake_files, output_names)
# End of wrangle snake files

# Here - Extract the interesting ones
# with molecular similarity
rds_to_select <- snake_files[grepl("19_molecular_similarity.rds", snake_files)]
prep <- snake_files[grepl("15_prep_annotation_biotransformation.rds", snake_files)]
# read_snakes(rds_to_select)
# read_snakes(prep)
# assign(
#   x = names(rds_to_select),
#   value = anno_sims_final,
#   envir = .GlobalEnv
# )
# rm(anno_sims_final)
# rm(apiin_bu_25_ppm_19_molecular_similarity.rds)

anno_sims_final %>%
  dplyr::filter(peak_id %in% all_sig_diff)

chem_pred_feats %>%
  dplyr::filter(feature %in% all_sig_diff)

# End of extract

# Import all snake files
purrr::walk(
  .x = snake_files,
  .f = ~ {
    read_snakes(.x)
    cli::cli_alert_success("Imported {.val {.x}}")
  }
)
# End of import all snake files

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
# Probably not needed now
# source("analysis/02_produce_chromatograms.R")
# TODO
# fix this part
# source("analysis/03_plotting_features.R")
