
suppressWarnings(
    suppressPackageStartupMessages({
        library(tidyverse)
        library(writexl)
        library(readxl)
    })
)

setwd(file.path("C:/Users/wilhelm/Documents/MEGA/01_juniper/01_arbete/01_projekt/03_psm"))
setwd(file.path("C:/Users/wilhelm/Documents/MEGA/01_juniper/01_arbete/01_projekt/03_psm", "parse_kegg"))

rpair <- read_tsv(
    file.path(getwd(), "kegg_reaction_classes.tsv"),
    show_col_types = FALSE
)

rpair.filt <- rpair %>%
    dplyr::mutate(rpair = strsplit(rpair, "\\s", perl = TRUE)) %>%
    tidyr::unnest(cols = "rpair") %>%
    dplyr::mutate(
        comp1 = stringr::str_split_i(rpair, "_", 1),
        comp2 = stringr::str_split_i(rpair, "_", 2)
    ) %>%
    dplyr::select(c(
        "entry",
        "rpair",
        "comp1",
        "comp2",
        "definition",
        "reaction",
        "enzyme",
        "pathway",
        "orthology"
    ))

compounds <- unique(c(rpair.filt$comp1, rpair.filt$comp2))
compounds.split <- split(compounds, ceiling(seq_along(compounds) / 10))
comp.tib <- tibble::tibble()
for (i in seq_along(compounds.split)) {
    
    message("Getting compounds: ", paste0(compounds.split[[i]], collapse = ", "))
    data <- KEGGREST::keggGet(dbentries = compounds.split[[i]])
    
    tmp.tib <- dplyr::bind_rows(
        lapply(
            X = data,
            FUN = function(x) {
                
                tib <- tibble::tibble(
                    entry = x$ENTRY,
                    name = paste0(x$NAME, collapse = " "),
                    formula = x$FORMULA,
                    dblinks = list(x$DBLINKS)
                )
                return(tib)
            }
        )
    )
    comp.tib <- bind_rows(comp.tib, tmp.tib)
}

saveRDS(object = final.tib, file = file.path(getwd(), "compounds.rds"))

final.tib <- comp.tib %>%
    dplyr::filter(!grepl("[nR()]", formula))

removed <- comp.tib %>%
    dplyr::filter(grepl("[nR()]", formula))
message("Removed ", nrow(removed), " rows with 'n', 'R', '(', or ')'")

reac.pair <- rpair.filt %>%
    dplyr::left_join(
        x = .,
        y = final.tib %>%
            dplyr::select(entry, formula) %>%
            dplyr::rename("formula_1" = "formula"),
        by = c("comp1" = "entry")
    ) %>%
    dplyr::left_join(
        x = .,
        y = final.tib %>%
            dplyr::select(entry, formula) %>%
            dplyr::rename("formula_2" = "formula"),
        by = c("comp2" = "entry")
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
        monomass_1 = Rdisop::getMass(Rdisop::getMolecule(formula_1)),
        monomass_2 = Rdisop::getMass(Rdisop::getMolecule(formula_2))
    ) %>%
    dplyr::mutate(
        delta_monomass = abs(monomass_1 - monomass_2)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(c(
        "entry",
        "rpair",
        "comp1",
        "comp2",
        "formula_1",
        "formula_2",
        "monomass_1",
        "monomass_2",
        "delta_monomass",
        "definition",
        "reaction",
        "enzyme",
        "pathway",
        "orthology"
    ))

readr::write_tsv(
    x = reac.pair,
    file = file.path(getwd(), "rpairs.tsv")
)

message("Done.")