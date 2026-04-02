library(tidyverse)
library(readxl)
library(writexl)
source("scripts/functions.R")

glycoside_stds_path <- paste0(
  "/Volumes/bluecub/aglycone_release_100um_24h/data/standards/"
)

glycoside_stds <- readxl::read_xlsx(
  path = paste0(
    glycoside_stds_path,
    "032626_G&AG_STD_metadata_ws.xlsx"
  ),
  sheet = "glycoside_file info"
) %>%
  dplyr::mutate(
    file_name = gsub(
      ".d",
      ".mzML",
      `File name`,
      fixed = TRUE
    )
  )



# So many missing files so I can't do this
# Long needs to send me the files but actually non-empty ones now
glycosides_stds_files <- list.files(file.path(glycoside_stds_path, "mzml"))
glycoside_stds
glycoside_stds[glycoside_stds$file_name %in% glycosides_stds_files,]
missing_files <- glycoside_stds[!glycoside_stds$file_name %in% glycosides_stds_files,]

write_xlsx(
  x = missing_files,
  path = "missing_files.xlsx"
)


aglycone_stds <- readxl::read_xlsx(
  path = paste0(
    "/Volumes/bluecub/aglycone_release_100um_24h/data/standards/",
    "032626_G&AG_STD_metadata_ws.xlsx"
  ),
  sheet = "aglycone_file info"
) %>%
  dplyr::mutate(
    file_name = gsub(
      ".d",
      ".mzML",
      `Data File`,
      fixed = TRUE
    )
  )

aglycone_stds_files <- list.files(
  file.path(
    paste0(
      "/Volumes/bluecub/aglycone_release_100um_24h/data/experiment/",
      "mzml"
    )
  )
)

aglycone_metadata <- aglycone_stds[
  aglycone_stds$file_name %in% aglycone_stds_files,
] %>%
  dplyr::select(
    "sample" = "file_name",
    "group" = "Sample",
    "unit",
    "bacteria"
  ) %>%
  dplyr::filter(unit > 50)

readr::write_csv(
  x = aglycone_metadata,
  file = paste0(
    "/Volumes/bluecub/aglycone_release_100um_24h/standards/",
    "aglycone_metadata.csv"
  )
)

# aglycones are in the experiment path

stds <- import_mzml(
  data_path =   file.path(
    paste0(
      "/Volumes/bluecub/aglycone_release_100um_24h/data/experiment/",
      "mzml"
    )
  ),
  meta_file = paste0(
    "/Volumes/bluecub/aglycone_release_100um_24h/standards/",
    "aglycone_metadata.csv"
  )
)

stds_exp <- MsExperiment::readMsExperiment(
  spectraFiles = stds$path,
  sampleData = stds
)

stds_exp