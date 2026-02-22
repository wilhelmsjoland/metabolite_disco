source("scripts/01_setup.R")
source("scripts/02_bpc.R")
source("scripts/03_internal_standard.R")
source("scripts/04_peak_calling.R")
source("scripts/05_alignment.R")

# TODO 03_internal_standard.R
# Isn't the last step of creating a new chromatogram unnecessary?
# Can't I just use on of the chromatograms I used before?
# 1. full
# 2. narrow
# 3. wide
# 4. call peaks on wide
# 5. take peakwidths from wide
# TODO
# Potentially fix so that parameters can be chosen as well
# TODO
# Change the first peak grouping to more broad
# so I don't get a million anchor peaks?
# TODO
# Could add these as names for the list
# xchr5@processHistory[[2]]@type
# TODO
# Could potentially move the BPC stuff from 02 to 05 instead?
# Was in the 05_script
# - From Sattely paper
# - Retention time correction was performed using the obiwarp method, with a
# - step size of m/z 0.5. Peak alignment was performed with bandwidth
# - of 3 seconds and minimum fraction (minfrac) of samples
# - necessary for a valid group of 0.5.
# TODO
# Add control flow for using filtering in 08_filter_features or not

# TODO
# Anything under has not been formatted or fixed yet.
source("scripts/06_correspondence.R")
source("scripts/07_gap_filling.R")
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
