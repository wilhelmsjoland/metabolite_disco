library(tidyverse)
library(readxl)
library(CompoundDb)
library(MetaboAnnotation)

db_path <- paste0(
  "/Users/wilhelm/Proton/01_juniper/01_arbete/01_projekt/03_psm/molekyler/",
  "information/plant_small_molecules_db"
)
db_out <- "/Volumes/bluecub/databases/psm"

db_info <- readxl::read_xlsx(
  path = file.path(db_path, "psm_standards.xlsx"),
  sheet = "standards",
  col_types = c(
    c(
      rep("text", 7),
      "numeric",
      rep("text", 2),
      rep("numeric", 2)
    )
  )
)

cmps <- db_info %>%
  tidyr::drop_na(cid)

metadata <- CompoundDb::make_metadata(
  source = "psm_standards",
  url = "",
  source_version = "1.0",
  source_date = as.character(Sys.Date()),
  organism = NA_character_
)

db_file <- CompoundDb::createCompDb(
  x = cmps,
  metadata = metadata,
  path = db_out
)

cdb <- CompoundDb::CompDb(db_file)
