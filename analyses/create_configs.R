dir.create("configs", FALSE, TRUE)

# TODO
# Now create the config file for every experiment.....

yam <- yaml::read_yaml("config.yaml")

test <- experiments %>%
  dplyr::filter(experiment == 1)

# make a species columns as the name of the experiment
# make a glycoside columns as the name of the experiment
# glycoside_species

# TODO
# the metadata files need to be split and saved somewhere

list(
  output = paste0("output", "column"),
  data_path = exp_path,
  meta_file = "",
  internal_standard = "", # Check with long
  is_adduct = "", # Check this first with standard_disco.R
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
  biotransf_file = "biotransformations.csv",
  rpairs_path = "scripts/search_compounds/output/rpairs.tsv",
  glycoside = "", # fetch from glycoside name
  aglycone = "", # fetch from glycoside name
  ppm_match = 5,
  all_vs_all = "false",
  smiles =  "", # fetch from glycoside name
  biot_dir = "biotransformer3.0jar",
  seed = 123,
  cores = 10
)