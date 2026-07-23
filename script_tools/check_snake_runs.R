library(tidyverse)

file_names <- list.files(
  path = "/Volumes/bluecub/aglycone_release_100um_24h/output",
  full.names = TRUE,
  recursive = TRUE,
  pattern = "*.pipeline.log"
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
        any(stringr::str_detect(lines, fixed("Saved molecular similarity of biotransformer predicted features to:"))),
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

# A tibble: 3 × 2
#   experiment                                                                           V1    
#   <chr>                                                                                <chr> 
# 1 clitorin_b_thetaiotaomicron_vpi_5482                                                 Failed
# 2 quercetin_7_o_b_d_glucopyranoside_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon Failed
# 3 rutin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon                             Failed


file_names %>%
  readr::read_lines() %>%
  stringr::str_subset(., fixed("Saved molecular similarity of biotransformer"))

file_names %>% length
