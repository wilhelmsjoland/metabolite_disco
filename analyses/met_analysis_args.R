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
    c("--input"),
    type = "character",
    default = "output/apiin_bu_25_ppm",
    help = "Input data directory [default: %default]"
  ),
  optparse::make_option(
    c("--output"),
    type = "character",
    default = "output/apiin_bu_25_ppm/analysis",
    help = "Output data directory [default: %default]"
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
    default = 1, # parallel::detectCores() - 1,
    help = "Workers, will default to max available - 1 [default: %default]."
  )
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

setwd(opt$wd)
set.seed(opt$seed)
register_parallel(opt$cores)