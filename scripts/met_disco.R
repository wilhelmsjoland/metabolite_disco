# ==============================================================================
# TODO Add optparse specifications here -> make it a func ----------------------
# ==============================================================================
source("scripts/01_setup.R")
source("scripts/02_bpc.R")
# TODO
# Isn't the last step of creating a new chromatogram unnecessary?
# Can't I just use on of the chromatograms I used before?
# 1. full
# 2. narrow
# 3. wide
# 4. call peaks on wide
# 5. take peakwidths from wide
# ==============================================================================
source("scripts/03_internal_standard.R")

# TODO
# Anything under has not been formatted or fixed yet.
# TODO
# Potentially fix so that parameters can be chosen as well
# TODO
# Extract some peaks during peka inspection and check quality of peak picking
# TODO
# potentially remove or move the inspecting of called peaks?
source("scripts/04_peak_calling.R")
source("scripts/05_alignment.R")
source("scripts/06_correspondence.R")
source("scripts/07_gap_filling.R")
# TODO
# Add control flow for using filtering or not
source("scripts/08_filter_features.R")
source("scripts/09_scaling.R")
source("scripts/10_pca.R")
source("scripts/11_limma.R")
source("scripts/12_volcano.R")
source("scripts/13_upset.R")
source("scripts/14_intersecting_features.R")
source("scripts/15_prep_annotation_biotransformation.R")
source("scripts/16_mz_predictions.R")
source("scripts/17_annotation.R")
source("scripts/18_biotransformer.R")
source("scripts/19_molecular_similarity.R")
source("scripts/20_produce_chromatograms.R")
source("scripts/21_plotting_features.R")
