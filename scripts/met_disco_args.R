option_list <- list(
  optparse::make_option(
    c("--wd"),
    type = "character",
    default = paste0(
      "C:/Users/wilhelm/Documents/MEGA/",
      "01_juniper/01_arbete/01_projekt/03_psm"
    ),
    help = "Working directory [default: %default]"
  ),
  optparse::make_option(
    c("--biot_dir"),
    type = "character",
    default = paste0(
      "C:/Users/wilhelm/Documents/MEGA/",
      "01_juniper/01_arbete/01_projekt/03_psm/biotransformer3.0jar"
    ),
    help = "Biotransformer directory [default: %default]"
  ),
  optparse::make_option(
    c("--smiles"),
    type = "character",
    default = "C1=CC(=CC=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)O",
    help = "SMILES for prediction of metabolites [default: %default]"
  ),
  optparse::make_option(
    c("-d", "--data_path"),
    type = "character",
    default = "input_data",
    help = "Input data directory [default: %default]"
  ),
  optparse::make_option(
    c("-m", "--meta_file"),
    type = "character",
    default = "input_data/metadata_filtered.csv",
    help = "Metadata CSV [default: %default]"
  ),
  optparse::make_option(
    c("-b", "--biotransf_file"),
    type = "character",
    default = "biotransformations.csv",
    help = "Biotransformations.csv file path [default: %default]"
  ),
  optparse::make_option(
    c("-o", "--output"),
    type = "character",
    default = "apiin_bu",
    help = "Results folder [default: %default]"
  ),
  optparse::make_option(
    c("-f", "--internal_standard"),
    type = "character",
    default = "C7H8O2",
    help = "Chemical formula of the internal standard [default: %default]"
  ),
  optparse::make_option(
    c("-p", "--polarity"),
    type = "character",
    default = "negative",
    help = "Polarity for adduct formation [default: %default]"
  ),
  optparse::make_option(
    c("--glycoside"),
    type = "character",
    default = "C26H28O14",
    help = "Chemical formula of the glycoside [default: %default]"
  ),
  optparse::make_option(
    c("--aglycone"),
    type = "character",
    default = "C15H10O5",
    help = "Chemical formula of the aglycone [default: %default]"
  ),
  optparse::make_option(
    c("-a", "--is_adduct"),
    type = "character",
    "default" = "[M-H]-",
    help = "Adduct of the internal standard [default: %default]"
  ),
  optparse::make_option(
    c("--ppm_global"),
    type = "double",
    default = 15,
    help = "Global ppm tolerance [default: %default]"
  ),
  optparse::make_option(
    c("--ppm_match"),
    type = "double",
    default = 5, # 10
    help = "ppm tolerance for matching [default: %default]"
  ),
  optparse::make_option(
    c("--bw_first_grouping"),
    type = "double",
    default = 3,
    help = "Bandwidth for first grouping [default: %default]"),
  optparse::make_option(
    c("--bw_second_grouping"),
    type = "double",
    default = 0.5,
    help = "Bandwidth for second grouping [default: %default]"),
  optparse::make_option(
    c("-x", "--seed"),
    type = "integer",
    default = 123,
    help = "Random seed [default: %default]"
  ),
  optparse::make_option(
    c("-c", "--cores"),
    type = "integer",
    default = parallel::detectCores() - 1,
    help = "Workers, will default to max available - 1 [default: %default]."
  ),
  optparse::make_option(
    c("-l", "--missingness"),
    type = "integer",
    default = 50, # 100 # 70 # 50 is okay since it's per group with f
    help = paste0(
      "Filter features based on proportion missingness per group (0 - 100)",
      " [default: %default]."
    )
  ),
  optparse::make_option(
    c("-q", "--qvalue"),
    type = "double",
    default = 0.05,
    help = "P value to use for linear models and plots [default: %default]."
  ),
  optparse::make_option(
    c("--sn_threshold"),
    type = "double",
    default = 10,
    help = "Signal: Noise threshold [default: %default]."
  ),
  optparse::make_option(
    c("--beta_cor_threshold"),
    type = "double",
    default = 0.8, # 0.3
    help = "Beta correlation threshold [default: %default]."
  ),
  optparse::make_option(
    c("--beta_snr_threshold"),
    type = "double",
    default = 3, # 6
    help = "Beta snr threshold [default: %default]."
  ),
  optparse::make_option(
    c("--gap_filling"),
    type = "character",
    default = "norm_fill_imp",
    help = paste0(
      "Can choose between 'norm', 'norm_fill', 'norm_fill_imp'",
      " [default: %default]."
    )
  ),
  optparse::make_option(
    c("--rpairs_path"),
    type = "character",
    default = "scripts/search_compounds/output/rpairs.tsv",
    help = "Path of rpairs.tsv [default: %default]."
  )
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

setwd(opt$wd)
set.seed(opt$seed)
register_parallel(opt$cores)

# ==============================================================================
# Parse validity of options in optparse list -----------------------------------
# ==============================================================================
valid_gap_fill <- c("norm", "norm_fill", "norm_fill_imp")
if (!opt$gap_filling %in% valid_gap_fill) {
  cli::cli_abort(
    paste0(
      "--gap_filling must be one of: ",
      paste(valid_gap_fill, collapse = ", ")
    )
  )
}
