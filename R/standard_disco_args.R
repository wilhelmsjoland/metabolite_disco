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
        default = "./standards/apiin_stds",
        help = "Input data directory [default: %default]"
    ),
    optparse::make_option(
        c("-m","--meta_file"), 
        type = "character", 
        default = "./standards/apiin_stds/apiin_stds_meta.csv",
        help = "Metadata CSV [default: %default]"
    ),
    optparse::make_option(
        c("-o","--output"), 
        type = "character", 
        default = "./standard_peaks/apiin_standards",
        help = "Results folder [default: %default]"
    ),
    optparse::make_option(
        c("-f", "--standard"),
        type = "character",
        "default" = "C26H28O14",
        help = "Chemical formula of the internal standard [default: %default]"
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
        c("--expand_rt"), 
        type = "double", 
        default = 8,
        help = "Expansion of retention time of chromatograms in seconds [default: %default]"
    ),
    optparse::make_option(
        c("--expand_mz"), 
        type = "double", 
        default = 0.05,
        help = "Expansion of m/z in chromatograms [default: %default]"
    ),
    optparse::make_option(
        c("-x","--seed"), 
        type = "integer",
        default = 123,
        help = "Random seed [default: %default]"
    ),
    optparse::make_option(
        c("-c","--cores"), 
        type = "integer", 
        default = parallel::detectCores() - 2,
        help = "Workers, will default to max available - 2 [default: %default]."
    )
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

setwd(opt$wd)
set.seed(opt$seed)
data.path <- opt$data_path
meta.file <- opt$meta_file
stds.output.path <- opt$output
glycoside.form <- opt$standard
std.ppm <- opt$ppm_global
std.adduct <- opt$adduct
expand.rt <- opt$expand_rt
expand.mz <- opt$expand_mz
workers <- opt$cores

registerParallel(workers)