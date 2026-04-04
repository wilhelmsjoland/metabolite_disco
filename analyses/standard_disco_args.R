option_list <- list(
  optparse::make_option(
    c("--wd"),
    type = "character",
    default = paste0(
      "/Users/wilhelm/MEGA/01_juniper/01_arbete/01_projekt/03_psm"
    ),
    help = "Working directory [default: %default]"
  ),
  optparse::make_option(
    c("-d", "--data_path"),
    type = "character",
    default = paste0(
      "/Volumes/bluecub/aglycone_release_100um_24h/",
      "data/experiment/mzml"
    ),
    help = "Input data directory [default: %default]"
  ),
  optparse::make_option(
    c("-m", "--meta_file"),
    type = "character",
    default = paste0(
      "/Volumes/bluecub/aglycone_release_100um_24h/standard_output/",
      "aglycone_metadata.csv"
    ),
    help = "Metadata CSV [default: %default]"
  ),
  optparse::make_option(
    c("-o", "--output"),
    type = "character",
    default = paste0(
      "/Volumes/bluecub/aglycone_release_100um_24h/standard_output/aglycone"
    ),
    help = "Results folder [default: %default]"
  ),
  optparse::make_option(
    c("--polarity"),
    type = "character",
    default = "negative",
    help = "Polarity [default: %default]"
  ),
  optparse::make_option(
    c("-f", "--is_formula"),
    type = "character",
    "default" = "C7H8O2",
    help = "Chemical formula of the internal standard [default: %default]"
  ),
  optparse::make_option(
    c("-a", "--is_adduct"),
    type = "character",
    "default" = "[M-H]-",
    help = "Major adduct of the internal standard [default: %default]"
  ),
  optparse::make_option(
    c("--ppm_global"),
    type = "double",
    default = 25,
    help = "Global ppm tolerance [default: %default]"
  ),
  optparse::make_option(
    c("-x", "--seed"),
    type = "integer",
    default = 123,
    help = "Random seed [default: %default]"
  ),
  optparse::make_option(
    c("-c", "--cores"),
    type = "integer",
    default = 4,
    help = "Workers, will default to max available - 2 [default: %default]."
  )
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

setwd(opt$wd)
set.seed(opt$seed)