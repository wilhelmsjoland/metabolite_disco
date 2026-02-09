# TODO FIX the function for this
# it already worked separetely
run_biotransformer <- function() {
}

old_wd <- getwd()
new_wd <- file.path(old_wd, "predict_metabolism", "biotransformer3.0jar")

output_name <- "apigenin2.csv"
biot_output_loc <- file.path(old_wd, res.folder, "tables")

cmd <- paste(
  "java -jar BioTransformer3.0_20230525.jar",
  "-k pred",
  "-b superbio",
  # "-isdf", biot_output_loc,
  "-ismi C1=CC(=CC=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)O",
  "-ocsv", paste0(biot_output_loc, "/", output_name),
  "-a"
)

setwd(new_wd)
biot_output <- system(cmd, intern = TRUE)
setwd(old_wd)

log_file <- file.path(biot_output_loc, "biotransformer.log")
writeLines(biot_output, log_file)
