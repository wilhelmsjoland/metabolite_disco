option_list <- list(
  optparse::make_option(
    c("--wd"),
    type = "character",
    default = "C:/Users/wilhelm/Documents/MEGA/01_juniper/01_arbete/01_projekt/03_psm",
    help = "Working directory [default: %default]"
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
    help = "Biotransformations file (relative to data_path or absolute) [default: %default]"
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
    c("--glycoside_ppm"),
    type = "double",
    default = 2000,
    help = "ppm tolerance for glycoside/aglycone search [default: %default]"
  ),
  optparse::make_option(
    c("-a", "--adduct"),
    type = "character",
    "default" = "[M-H]-",
    help = "Chemical formula of the internal standard [default: %default]"
  ),
  optparse::make_option(
    c("--ppm_global"),
    type = "double",
    default = 15,
    help = "Global ppm tolerance [default: %default]"
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
    default = parallel::detectCores() - 2,
    help = "Workers, will default to max available - 2 [default: %default]."
  ),
  optparse::make_option(
    c("-l", "--missingness"),
    type = "integer",
    default = 100,
    help = paste0(
      "Filter features based on missingness from 0 - 100 percent",
      " [default: %default]."
    )
  ),
  optparse::make_option(
    c("-p", "--pvalue"),
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
    default = 0.3,
    help = "Beta correlation threshold [default: %default]."
  ),
  optparse::make_option(
    c("--beta_snr_threshold"),
    type = "double",
    default = 6,
    help = "Beta snr threshold [default: %default]."
  )
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

# ---- deterministic working dir ----
setwd(opt$wd)
set.seed(opt$seed)

# ---- map options to variable names ----

data_path <- opt$data_path
meta_file <- opt$meta_file
biotransf_file <- opt$biotransf_file
res_folder <- opt$output
internal_standard <- opt$internal_standard
glycoside_form <- opt$glycoside
aglycone_form <- opt$aglycone
glycoside_ppm <- opt$glycoside_ppm
adduct <- opt$adduct
polarity <- opt$polarity
ppm_global <- opt$ppm_global
bw_first_grouping <- opt$bw_first_grouping
bw_second_grouping <- opt$bw_second_grouping
workers <- opt$cores
missing_threshold <- opt$missingness
p_value_global <- opt$pvalue
sn_threshold <- opt$sn_threshold
beta_cor_threshold <- opt$beta_cor_threshold
beta_snr_threshold <- opt$beta_snr_threshold

register_parallel(workers)
