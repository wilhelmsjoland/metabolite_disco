library(RMariaDB)
library(AnnotationHub)

# met.to.check <- "C39H44N4O12"
# 
# tibble::as_tibble(featureDefinitions(xchr9), rownames = "feature") %>%
#     dplyr::filter(
#         between(
#             mzmed,
#             getMass(getMolecule(met.to.check)) - ppm_to_num(15),
#             getMass(getMolecule(met.to.check)) + ppm_to_num(15)
#         )
#     )

# TRYING with annotationhub
# from https://bioconductor.org/packages/release/bioc/vignettes/AHMassBank/inst/doc/creating-MassBank-CompDbs.html

ah <- AnnotationHub()
query(ah, "MassBank")
qr <- query(ah, c("MassBank", "2024.11"))
cdb <- qr[[1]]

# Creating a CompDb from MassBank
# con <- dbConnect(MariaDB(), host = "localhost", user = <username>,
#                  pass = <password>, dbname = "MassBank")
# source(system.file("scripts", "massbank_to_compdb.R", package = "CompoundDb"))
# massbank_to_compdb(con)

anno.hmdb %>%
    dplyr::arrange(mz) %>% 
    dplyr::filter(adduct == "[M-H]-") %>%
    view()

anno.mass %>%
    dplyr::arrange(mz) %>% 
    dplyr::filter(adduct == "[M-H]-") %>%
    view()

anno %>%
    dplyr::arrange(mz) %>% 
    dplyr::filter(adduct == "[M-H]-") %>%
    view()

anno %>%
    dplyr::filter(adduct == "[M-H]-") %>%
    view()
