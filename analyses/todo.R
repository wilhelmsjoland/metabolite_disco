################################################################################
# High priority
################################################################################
# TODO
# Fix under "14. plotting interesting feature"
# Filter harder!! and make this a selection in the original pipeline call
# Update the annotation database to the prepackaged full one!
# Later fix and build my own one


# Create a a table of all genes that are different in all comparisons
# Add more features to the biotransformation dataset
# DO it through the RHEA thing -> check GPT suggestion
# Make a function of the standard checking
# My implementation of biotransformation checks can
# right now only handle if they are 1 negative/positive charge
# Check all the chromatogram objects and see if
# they include filled peaks with the parameter "filled = TRUE"
# Filter what I put into the feature.chrs object
# Mabye based on beta gaussian etc etc or something?
# Fix so the path isn't hardcoded in check_saved() and
# that it follows any path input
# Change the plotting in the standard_disco to ggplots ->
# take the idea from 251230_sqarc.R
# Do this also for the plotting functions -> may be faster
# Add the content from 251230_sqarc.R into the standard_disco.R and met_disco.R
# Fix the glycoside script so I can always use that as well!
# Programatically create the biotransformation dataset from some database
# Look through the annotation file and look what molecular formulas overlap
# and potentially make sense depending on the adduct that is added or removed?
# Check so delta m/z is correct for glycosides and aglycones
# Find a physical harddrive on amazon and send to Gabriel!!

# add the observed ppm diff to the biotransformations file
# In the output I could also add what other transformations it could be
# Using the Rdisop decomposemass and get the formula
# Compare the patterns in mz1_forms vs mz2_forms to determine
# what formulas even can match and make sense!
# Filter significant ones based on the delta mass dataframe

################################################################################
# Low priority
################################################################################

# Läs om glykosider osv
# Mejla Trevor/Nicole Levesque om slack osv....
# Move the database storage from a temporary location
# Add a setting where I specify the order of the groups
# Calculate the min and max peakwidth for several peaks
# and not just for the IS as it is now
# - Look into how I should determine molecular fingerprinting
# - Molecular similiarity - Tanimoto

# TODO read- CentWaveParam()
# Create an estimation of the noise in each sample
# and use as input for the Peak calling
# Look into different ways of clustering m/zs and if they are related