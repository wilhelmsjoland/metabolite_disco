output_folders <- list.files(
  "/Volumes/bluecub/aglycone_release_100um_24h/output",
  full.names = TRUE
)

test <- file.path(
  output_folders[1],
  "snakemake_objects",
  "19_molecular_similarity.rds"
)

test2 <- readRDS(test)

test2
