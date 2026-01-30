library(Rdisop)
library(tidyverse)

tmp.all <- c()

# TODO
# Map this instead of nested loops
# Prefilter based on maximum masses
# Fix so that C0's mass is 0

for (i in 0:5) {
    for (j in 0:5) {
        for (k in 0:5) {
            for (l in 0:5) {
                tmp <- paste0("C", i, "H", j, "N", k, "O", l)
                tmp.all <- c(tmp.all, tmp)
            }
        }
    }
}


tibble::tibble(
    "name" = tmp.all,
    "chem_formula" = paste0("± ", tmp.all),
    "allowed_n" = 1,
    "multiplier" = 1,
    ) %>%
    # Speed this up instead of doing rowwise
    rowwise() %>%
    dplyr::mutate(delta_mass = Rdisop::getMass(Rdisop::getMolecule(name))) %>%
    ungroup()

bio.transf

Rdisop::getMass(Rdisop::getMolecule("C0H0N4O11"))


# TODO
# Define some reasonable rules for this and what to use
try <- expand.grid(
    C = 0:10,
    H = 5:10,
    N = 0:8, 
    O = 0:8, 
    P = 0:10,
    S = 0:10
    ) %>%
    dplyr::mutate(chem_form = paste0(
        "C", C, 
        "H", H, 
        "N", N, 
        "O", O, 
        "P", P, 
        "S", S
    )) %>%
    tibble::as_tibble(.)

try2 <- try %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
        chem_form = gsub("([A-Z][a-z]?)0(?![0-9])", "", chem_form, perl = TRUE)
    ) %>%
    dplyr::mutate(
        chem_form2 = if (chem_form == "") {
            NA
        } else {
            Rdisop::getMass(Rdisop::getMolecule(chem_form))
        }
    )

try2
system.time(
getMass(getMolecule("C0H0N0O0S0"))
)

