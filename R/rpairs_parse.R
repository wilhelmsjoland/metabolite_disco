library(tidyverse)

data <- read_tsv("parse_kegg/search_compounds/output/rpairs.tsv", show_col_types = FALSE)

df <- data %>%
    dplyr::select(entry, rpair, comp1, comp2, name1, name2, delta_mass) %>%
    # dplyr::distinct(delta_mass, .keep_all = TRUE) %>%
    dplyr::mutate(number = row_number()) %>%
    dplyr::relocate(number, .before = "entry") %>%
    dplyr::mutate(
        across(
            .cols = all_of(c("name1", "name2")),
            .fns = ~ stringr::str_trim(stringr::str_extract(., "^[^;]+"))
        )
    )

df %>%
    print(n = "all")

nrow(df) * 17
