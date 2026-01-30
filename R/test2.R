library(tidyverse)

data <- read_tsv(file.path(
  "parse_kegg",
  "search_compounds",
  "output",
  "rpairs.tsv"
))

data |>
  dplyr::select(entry, rpair, comp1, comp2, name1, name2, delta_mass) |>
  # dplyr::distinct(delta_mass, .keep_all = TRUE) %>%
  dplyr::mutate(number = row_number()) |>
  dplyr::relocate(number, .before = "entry") |>
  dplyr::mutate(
    dplyr::across(
      .cols = all_of(c("name1", "name2")),
      .fns = ~ stringr::str_trim(stringr::str_extract(., "^[^;]+"))
    )
  ) |>
  dplyr::filter(
    dplyr::if_any(
      .cols = all_of(c("name1", "name2")),
      .fns = ~ stringr::str_detect(., regex("apige", TRUE))
    )
  ) |>
  dplyr::distinct(delta_mass, .keep_all = TRUE)
