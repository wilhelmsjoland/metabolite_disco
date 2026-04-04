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

# Det ska vara 48 filer totalt

# So many missing files so I can't do this
# Long needs to send me the files but actually non-empty ones now
glycosides_stds_files <- list.files(file.path(glycoside_stds_path, "mzml"))


glycosides_stds_files_df <- tibble::tibble(
  actual_file = glycosides_stds_files
) %>%
  dplyr::mutate(
    file_name = stringr::str_remove(actual_file, "-r\\d+(?=\\.mzML$)")
  )

glycoside_stds_matched <- glycoside_stds %>%
  dplyr::left_join(y = glycosides_stds_files_df, by = "file_name")

"sample,group,unit,bacteria,formula" 

glycoside_stds_matched %>%
  dplyr::select(
    "long_spec_d" = `File name`,
    "long_spec" = "file_name",
    "sample" = "actual_file",
    "group" = `Glycoside group`,
    "bacteria" = "Strain"
  )

writexl::write_xlsx(
  x = glycoside_stds_matched,
  path = "glycoside_std_metadata.xlsx"
)

# Seems like the files are still missing but I have enough now at least to
# find the standard locations


glycoside


glycoside_stds %>% print(n = 20)
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
    "aglycone_metadata_automated.csv"
  )
)

# aglycones are in the experiment path

stds <- import_mzml(
  data_path = file.path(
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
