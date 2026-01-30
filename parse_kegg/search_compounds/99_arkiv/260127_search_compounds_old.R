args <- commandArgs(trailingOnly = TRUE)
cmd.call <- commandArgs()

input.idx <- match(TRUE, args %in% c("-i", "--input"))
outdir.idx <- match(TRUE, args %in% c("-o", "--output"))

input <- args[input.idx + 1]
outdir <- args[outdir.idx + 1]

run.flag <- (length(args) == 4) &&
            (!is.na(input.idx) && !is.na(outdir.idx)) &&
            (input.idx < length(args) && outdir.idx < length(args)) &&
            (!input %in% c("-i", "--input", "-o", "--output")) &&
            (!outdir %in% c("-i", "--input", "-o", "--output"))

if (isFALSE(run.flag)) {
    stop(
        "\nUsage: \n\'Rscript -i <input> -o <output>' or \n'Rscript --input <input> --output <output>'\n",
        "\nYou need to supply exactly two arguments:\n\t", 
        "input with either '-i', or '--input'\n\t",
        "output with either '-o' or '--output'\n\n\t"
        )
}

suppressWarnings(
    suppressPackageStartupMessages({
        library(dplyr)
        library(KEGGREST)
    })
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
time.run <- format(Sys.time(), usetz = TRUE)

message("Downloading database info...")
for (i in c("compound", "reaction", "rclass")) {
    writeLines(
        text = KEGGREST::keggInfo(i),
        con = file.path(outdir, paste0("db_", i, "_info.txt"))
    )
}

data <- read.delim(
    file.path(input), 
    header = FALSE, 
    col.names = "query",
    comment.char = "#"
    )

message("Finding compounds related to search terms...")
n.iters <- length(data$query)
iter.chars <- nchar(n.iters)
sprintf.fmt <- paste0("%0", iter.chars, "d")
compounds <- c()
for (i in seq_along(data$query)) {
    
    num.fmt <- sprintf(sprintf.fmt, i)
    message("[", num.fmt, "] ", "Fetching compounds related to: ", paste0(data$query[[i]], collapse = ", "))
    
    compound <- KEGGREST::keggFind("compound", data$query[[i]])
    compounds <- c(compounds, compound)
}

compounds.sort <- sort(unique(names(compounds)))

message("Finding mapping between compounds -> reactions -> reaction classes...")
compound.map <- KEGGREST::keggLink("reaction", compounds.sort)
reaction.map  <- KEGGREST::keggLink("rclass", compound.map)
reaction.map.sort <- sort(unique(reaction.map))

ten.rclass.list <- split(reaction.map.sort, ceiling(seq_along(reaction.map.sort) / 10))
ten.rclass.plus <- lapply(
    X = ten.rclass.list, 
    FUN = function(x) {
        gsub("rc:", "", paste0(x, collapse = "+"))
    }
)

message("Downloading reaction classes...")
n.iters <- length(ten.rclass.plus)
iter.chars <- nchar(n.iters)
sprintf.fmt <- paste0("%0", iter.chars, "d")
rclasses <- list()
for (i in seq_along(ten.rclass.plus)) {
    
    num.fmt <- sprintf(sprintf.fmt, i)
    message("[", num.fmt, "] ", "Downloading: ", paste0(ten.rclass.plus[[i]], collapse = ", "))
    
    rclass <- KEGGREST::keggGet(ten.rclass.plus[[i]])
    rclasses <- c(rclasses, rclass)
}

message("Parsing reaction classes...")
rclass.tib <- dplyr::bind_rows(
    lapply(
        X = seq_along(rclasses),
        FUN = function (x) {
            tibble(
                entry = rclasses[[x]]$ENTRY,
                definition = list(rclasses[[x]]$DEFINITION),
                rpair = list(unlist(strsplit(rclasses[[x]]$RPAIR, "\\s+"))),
                reaction = list(unlist(strsplit(rclasses[[x]]$REACTION, "\\s+"))),
                enzyme = list(rclasses[[x]]$ENZYME),
                pathway = list(rclasses[[x]]$PATHWAY),
                orthology = list(rclasses[[x]]$ORTHOLOGY),
                rmodule = list(rclasses[[x]]$RMODULE)
            )
        }
    )
)

message("Filtering reaction pairs to compounds derived from search terms...")
filtered.rpairs <- rclass.tib %>%
    dplyr::select(c("entry", "rpair")) %>%
    tidyr::unnest(rpair) %>%
    dplyr::mutate(
        comp1 = stringr::str_split_i(rpair, "_", 1),
        comp2 = stringr::str_split_i(rpair, "_", 2)
    ) %>%
    dplyr::filter(
        dplyr::if_any(
            .cols = c("comp1", "comp2"),
            .fns = ~ .x %in% gsub("cpd:", "", compounds.sort)
        )
    )

# Why is this shorter?
# Probably because at least one of the compounds may not be in a specific reaction
# -> I need reactions where the molecule of interest in not necessarily in the name
# of both compounds in the reaction
filt.comps <- sort(unique(c(filtered.rpairs$comp1, filtered.rpairs$comp2)))
filt.comps.orig <- sort(unique(gsub("cpd:", "", compounds.sort)))

ten.comp.list <- split(filt.comps, ceiling(seq_along(filt.comps) / 10))
ten.comp.plus <- lapply(
    X = ten.comp.list, 
    FUN = function(x) {
        paste0(x, collapse = "+")
    }
)

message("Downloading compunds part of reaction pairs based on original search terms...")
n.iters <- length(ten.comp.plus)
iter.chars <- nchar(n.iters)
sprintf.fmt <- paste0("%0", iter.chars, "d")
compounds.flat <- list()
for (i in seq_along(ten.comp.plus)) {
    
    num.fmt <- sprintf(sprintf.fmt, i)
    message("[", num.fmt, "] ", "Downloading: ", paste0(ten.comp.plus[[i]], collapse = ", "))
    
    tmp.compounds <- KEGGREST::keggGet(ten.comp.plus[[i]])
    compounds.flat <- c(compounds.flat, tmp.compounds)
}

message("Parsing compunds...")
compounds.tib <- dplyr::bind_rows(
    lapply(
        X = seq_along(compounds.flat),
        FUN = function (x) {
            tibble::tibble(
                entry = compounds.flat[[x]]$ENTRY,
                name = paste0(gsub(";", "", compounds.flat[[x]]$NAME), collapse = "; "),
                formula = compounds.flat[[x]]$FORMULA,
                exact_mass = compounds.flat[[x]]$EXACT_MASS,
                mol_weight = compounds.flat[[x]]$MOL_WEIGHT,
                reaction = list(unlist(strsplit(compounds.flat[[x]]$REACTION, "\\s+"))),
                pathway = list(compounds.flat[[x]]$PATHWAY),
                enzyme = list(compounds.flat[[x]]$ENZYME),
                brite = list(compounds.flat[[x]]$BRITE),
                dblinks = paste0(compounds.flat[[x]]$DBLINK, collapse = "; "),
                atom = list(compounds.flat[[x]]$ATOM),
                bond = list(compounds.flat[[x]]$BOND)
            )
        }
    )
)

message("Combining reaction pairs tied to compounds derived from search terms...")
final.pairs <- filtered.rpairs %>%
    dplyr::left_join(
        x = .,
        y = compounds.tib %>%
            dplyr::select(entry, formula, name, dblinks) %>%
            dplyr::rename(
                "formula1" = "formula",
                "name1" = "name",
                "dblink1" = "dblinks"
                ), 
        by = c("comp1" = "entry")
    ) %>%
    dplyr::left_join(
        x = .,
        y = compounds.tib %>%
            dplyr::select(entry, formula, name, dblinks) %>%
            dplyr::rename(
                "formula2" = "formula",
                "name2" = "name",
                "dblink2" = "dblinks"
                ),
        by = c("comp2" = "entry")
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
        mass1 = Rdisop::getMass(Rdisop::getMolecule(formula1)),
        mass2 = Rdisop::getMass(Rdisop::getMolecule(formula2)),
        delta_mass = abs(mass1 - mass2)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(entry, rpair, comp1, comp2, formula1, formula2, mass1, mass2, delta_mass, name1, name2, dblink1, dblink2)

message("Writing results to .tsv...")
write.table(
    x = final.pairs,
    file = file.path(outdir, "rpairs.tsv"),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
)

output.info <- paste0(
    "Call: ", paste0(cmd.call, collapse = " "), "\n",
    "search_compounds.sh run on: ", time.run, "\n",
    file.path(outdir, "db_compound.info.txt"), ": information on compound database used\n",
    file.path(outdir, "db_reaction.info.txt"), ": information on reaction database used\n",
    file.path(outdir, "db_rclass.info.txt"), ": information on reaction class database used\n",
    file.path(outdir, "rpairs.tsv"), ": dataframe with reaction pairs matched to search terms\n",
    "Search terms used from ", input, ": \n\t", paste0(data$query, collapse = "\n\t")
)

writeLines(text = output.info, con = file.path(outdir, "output_info.txt"))

message("Done.")
