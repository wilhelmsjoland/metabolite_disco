library(tidyverse)

data <- read_tsv(
    file = "testing/chemical_data.tsv",
    col_types = paste0(strrep("c", 3), strrep("n", 4), "c", "l")
    )

data.filt <- data %>%
    dplyr::select(compound_id, formula)

rhea <- read_tsv(
    file = "testing/filt_rhea.tsv",
    col_types = paste0(strrep("c", 3))
    ) %>%
    dplyr::rename(
        "reaction" = `Reaction identifier`,
        "equation" = Equation,
        "chebi_id" = `ChEBI identifier`
        )

try <- rhea %>%
    dplyr::mutate(chebi_id = gsub("CHEBI:", "", chebi_id)) %>%
    dplyr::mutate(chebi_id = strsplit(chebi_id, ";")) %>%
    tidyr::unnest(cols = chebi_id) %>%
    dplyr::left_join(
        x = .,
        y = data.filt,
        by = c(chebi_id = "compound_id")
    )

rem.reac <- try %>%
    dplyr::filter(grepl("R", formula)) %>%
    dplyr::pull(reaction)

try2 <- try %>%
    dplyr::filter(!reaction %in% rem.reac) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(mono_mass = getMass(getMolecule(formula))) %>%
    dplyr::group_by(reaction, equation) %>%
    dplyr::summarize(
        chebi_id = list(chebi_id),
        formula = list(formula),
        mono_mass = list(mono_mass),
        .groups = "keep"
    )

test <- try2 %>%
    dplyr::filter(reaction %in% c("RHEA:83987")) %>%
    tidyr::unnest(cols = c(mono_mass))

# need to know which ones are on which side of the equation to do the 
# next steps here
# basically:

test[1,]$mono_mass + test[2,]$mono_mass
test[3,]$mono_mass + test[4,]$mono_mass

# Chem A + Chem B = Chem C + Chem D
# The total masses should be the same on both sides
# However Chem A + Chem B - Chem D = Chem C
test[1,]$mono_mass + test[2,]$mono_mass - test[4,]$mono_mass
test[3,]$mono_mass

# what is what (confusing since online rhea db shows different)
test[1,]$mono_mass # secologanin
test[2,]$mono_mass # serotonin
test[3,]$mono_mass # seco+sero
test[4,]$mono_mass # h20

# secologanin (1) + serotonin (2) = seco+sero (3) + h20 (4)
# leading to:
# Chem B - Chem D = Chem C - Chem A
test[2,]$mono_mass - test[4,]$mono_mass
test[3,]$mono_mass - test[1,]$mono_mass
# serotonin - h20 = seco+sero - secologanin

# Example: https://www.rhea-db.org/rhea/83987
# So the difference between serotonin and whats added to a molecule is:
getMass(getMolecule("C10H13N2O")) - getMass(getMolecule("H2O"))
getMass(getMolecule("C27H35N2O10")) - getMass(getMolecule("C17H24O10"))

# This means 159.0922 is added when serotonin is added to a molecule like this

# What if I have a bunch of reaction participants?
getFormula(decomposeMass(159.0922))

# try once
try3 <- try2[which(try2$reaction == "RHEA:83987"), ]
abs(try3$mono_mass[[1]][3] - try3$mono_mass[[1]][1])
abs(try3$mono_mass[[1]][2] - try3$mono_mass[[1]][4])

try2[4,]

try2[4,]$chebi_id

abs(try2[4,]$mono_mass[[1]][3] - try2[4,]$mono_mass[[1]][1])
abs(try2[4,]$mono_mass[[1]][2] - try2[4,]$mono_mass[[1]][4])

try2[4,]$mono_mass[[1]][1] - try2[4,]$mono_mass[[1]][5]

try2[4,]$mono_mass[[1]][3] + try2[4,]$mono_mass[[1]][4] - try2[4,]$mono_mass[[1]][2]
