library(tidyverse)

file_names <- list.files(
  path = "/Volumes/bluecub/aglycone_release_100um_24h/output",
  full.names = TRUE,
  recursive = TRUE,
  pattern = "*.snakemake.log"
)

file_names %>%
  split(x = ., f = dirname(dirname(.))) %>%
  purrr::map(
    .x = .,
    .f = ~ {
      lines <- .x %>%
        purrr::map(readr::read_lines) %>%
        unlist()

      dplyr::if_else(
        any(stringr::str_detect(lines, fixed("(100%) done"))),
        "Finished",
        "Failed"
      )
    }
  ) %>%
  setNames(basename(names(.))) %>%
  t() %>%
  t() %>%
  tibble::as_tibble(rownames = "experiment") %>%
  tidyr::unnest(V1) %>%
  dplyr::filter(V1 == "Failed") %>%
  print(n = Inf)

# pred_biot in mz_predictions seems to be the problem mostly
# one problem in hclust in bpc