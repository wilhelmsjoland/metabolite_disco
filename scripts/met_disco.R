source("scripts/01_setup.R")
source("scripts/02_bpc.R")
source("scripts/03_internal_standard.R")
source("scripts/04_peak_calling.R")
source("scripts/05_alignment.R")
source("scripts/06_correspondence.R")
source("scripts/07_gap_filling.R")
source("scripts/08_filter_features.R")
source("scripts/09_scaling.R")
source("scripts/10_limma.R")
source("scripts/11_pca.R")
source("scripts/12_volcano.R")
source("scripts/13_upset.R")
source("scripts/14_intersecting_features.R")
source("scripts/15_prep_annotation_biotransformation.R")
source("scripts/16_mz_predictions.R")
source("scripts/17_annotation.R")
source("scripts/18_biotransformer.R")
source("scripts/19_molecular_similarity.R")
end_log()

# TODO
# FILTERING XCHR9 too hard, losing apigenin!!!!!!!!!
# beta_cor isnt everything, mabye make ppm smaller?
# this is likely because the peak is wide as hell???
# Mabye fix the model they use OR
# use beta_snr OR beta_cor as cutoff - either one has to be GOOD?
# Another thing to scale down the biotransformations is finding
# a good way to classify what peaks are good ->
# mabye use the beta_snr and beta_cor more intelligently
# Should probably do both ppms at 5 and not one at 10 and one at 5
# Better to start looking at interesting features and annotate them!
# Look at what seth marked me in, in longs PR
# TODO
# THE UPSET_COMP COMPARISONS WORKS THE BEST HONESTLY
# FOR EVERY EXPERIMENT -> DEFINE THE INTERESTING COMPS!!!!!!!!!!!
# Run the similarity with several methods
# and combine them -> they all are important anyway
# ECPF4 & EPCF6

# Implement MsFeatures - grouping of features!!!!
# https://bioconductor.org/packages/3.22/bioc/manuals/MsFeatures/man/MsFeatures.pdf
# https://bioconductor.org/packages/3.22/bioc/vignettes/MsFeatures/inst/doc/MsFeatures.html
# https://sneumann.github.io/xcms/articles/LC-MS-feature-grouping.html
# Try CAMERA AGAIN!!!!!!!!!!!!!!!!!!

# TODO 03_internal_standard.R
# Isn't the last step of creating a new chromatogram unnecessary?
# Can't I just use on of the chromatograms I used before?
# 1. full
# 2. narrow
# 3. wide
# 4. call peaks on wide
# 5. take peakwidths from wide
# TODO - genreal
# Potentially fix so that parameters can be chosen as well
# TODO - 05
# Change the first peak grouping to more broad
# so I don't get a million anchor peaks?
# TODO - general
# Could add these as names for the list
# xchr5@processHistory[[2]]@type
# TODO - 05 -> 02
# Could potentially move the BPC stuff from 02 to 05 instead?
# TODO - 08
# Add control flow for using filtering in 08_filter_features or not
# TODO 10_limma
# Could make p.adjust.method a choice in optparse
# perhaps change the ordering of the left_joins in this script.
# x = .x, y = res_def, not opposite
# TODO 15_prep_annotation... .R
# search_compound is sourced from within the file.... -> FIX THIS
# TODO 16_mz_predictions.R
# Look at m/z predictions and distributions of them
# Filter the m/z predictions to the most reasonable..........
# The fill matched_diffs needs filtering before saving
# 3 Gigs is way too much for a table
