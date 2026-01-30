extract_brite <- function(txt) {
  # If txt contains literal backslash escapes,
  # turn them into real tabs/newlines
  txt <- gsub("\\\\t", "\t", txt)
  txt <- gsub("\\\\n", "\n", txt)
  lines <- strsplit(txt, "\n", fixed = TRUE)[[1]]
  # Match: leading level letter (B/C/D/...), then KEGG compound ID, then name
  pat <- "^\\s*[A-Z]\\s+(C\\d{5})\\s+(.+?)\\s*$"

  m <- regexec(pat, lines, perl = TRUE)
  hits <- regmatches(lines, m)
  hits <- hits[vapply(hits, length, integer(1)) == 3]  # keep only matches
  df <- tibble::tibble(
    compound = vapply(hits, `[`, "", 2),
    name = vapply(hits, `[`, "", 3))
  unique(df)
}

rate_limit <- local({
  last_time <- Sys.time() - 10
  function(max_per_sec = 2) {
    min_gap <- 1 / max_per_sec
    now <- Sys.time()
    gap <- as.numeric(difftime(now, last_time, units = "secs"))
    if (gap < min_gap) Sys.sleep(min_gap - gap)
    last_time <<- Sys.time()
  }
})

args <- commandArgs(trailingOnly = TRUE)
cmd.call <- commandArgs()

input.idx <- match(TRUE, args %in% c("-i", "--input"))
outdir.idx <- match(TRUE, args %in% c("-o", "--output"))
append.idx <- match(TRUE, args %in% c("-a", "--append"))

input <- args[input.idx + 1]
outdir <- args[outdir.idx + 1]
append <- args[append.idx + 1]

input  <- if (!is.na(input.idx)  && input.idx  < length(args)) args[input.idx  + 1] else NA
outdir <- if (!is.na(outdir.idx) && outdir.idx < length(args)) args[outdir.idx + 1] else NA
append <- if (!is.na(append.idx) && append.idx < length(args)) args[append.idx + 1] else NA

run.flag <- (length(args) == 4 || length(args) == 6) &&
    (!is.na(input.idx) && !is.na(outdir.idx)) &&
    (input.idx < length(args) && outdir.idx < length(args)) &&
    (is.na(append.idx) || append.idx < length(args)) &&
    (!input %in% c("-i", "--input", "-o", "--output", "-a", "--append")) &&
    (!outdir %in% c("-i", "--input", "-o", "--output", "-a", "--append")) &&
    (is.na(append.idx) || !append %in% c("-i", "--input", "-o", "--output", "-a", "--append"))

if (isFALSE(run.flag)) {
    stop(
        "\nUsage: 'Rscript -i <input> -o <output> -a <brite_categories>'",
        "\nYou need to supply either two or three arguments:\n\t", 
        "Required: input with either '-i', or '--input'\n\t",
        "Required: output with either '-o' or '--output'\n\t",
        "Optional: append with either '-a' or '--append'\n\t",
        "Example: 'Rscript search_compounds.R -i comps.tsv -o out -a br:br08021,br:br08005'"
    )
}

suppressWarnings(
    suppressPackageStartupMessages({
        library(dplyr)
        library(KEGGREST)
        library(Rdisop)
    })
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
time.run <- format(Sys.time(), usetz = TRUE)

message("Downloading database info...")
for (i in c("compound", "reaction", "rclass", "brite")) {
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
    rate_limit(2)
    num.fmt <- sprintf(sprintf.fmt, i)
    message("[", num.fmt, "] ", "Fetching compounds related to: ", paste0(data$query[[i]], collapse = ", "))
    
    compound <- KEGGREST::keggFind("compound", data$query[[i]])
    compounds <- c(compounds, compound)
}

compounds.sort <- sort(unique(names(compounds)))

if (length(args) == 6) {
    
    brites <- strsplit(args[6], ",")[[1]]
    
    message("Downloading supplied brite categories...")
    x.all <- c()
    for (i in brites) {
        rate_limit(2)
        # Fix the formatting here
        x <- KEGGREST::keggGet(i)[[1]]
        x.all <- c(x.all, x)
    }
    
    message("Parsing supplied brite categories")
    brite.compounds <- dplyr::bind_rows(
        lapply(
            X = x.all,
            FUN = function (x) {
                extract_brite(x)
            })
    )
    
    message("Appending compounds to search term compounds...")
    compounds.sort <- sort(unique(c(compounds.sort, paste0("cpd:", brite.compounds$compound))))
}

message("Finding mapping between compounds -> reactions -> reaction classes...")

ten.compound.list <- split(compounds.sort, ceiling(seq_along(compounds.sort) / 10))
ten.compound.plus <- lapply(
    X = ten.compound.list, 
    FUN = function(x) {
        gsub("rc:", "", paste0(x, collapse = "+"))
    }
)

n.iters <- length(ten.compound.plus)
iter.chars <- nchar(n.iters)
sprintf.fmt <- paste0("%0", iter.chars, "d")
compound.maps <- c()
for (i in seq_along(ten.compound.plus)) {
    rate_limit(2)
    num.fmt <- sprintf(sprintf.fmt, i)
    message("[", num.fmt, "] ", "Mapping: ", paste0(ten.compound.plus[[i]], collapse = ", "))
    
    compound.map <- KEGGREST::keggLink("reaction", ten.compound.plus[[i]])
    compound.maps <- c(compound.maps, compound.map)
}

sort.compound.maps <- sort(unique(compound.maps))

ten.reaction.list <- split(sort.compound.maps, ceiling(seq_along(sort.compound.maps) / 10))
ten.reaction.plus <- lapply(
    X = ten.reaction.list, 
    FUN = function(x) {
        gsub("rc:", "", paste0(x, collapse = "+"))
    }
)

n.iters <- length(ten.reaction.plus)
iter.chars <- nchar(n.iters)
sprintf.fmt <- paste0("%0", iter.chars, "d")
reaction.maps <- c()
for (i in seq_along(ten.reaction.plus)) {
    rate_limit(2)
    num.fmt <- sprintf(sprintf.fmt, i)
    message("[", num.fmt, "] ", "Mapping: ", paste0(ten.reaction.plus[[i]], collapse = ", "))
    
    reaction.map <- KEGGREST::keggLink("rclass", ten.reaction.plus[[i]])
    reaction.maps <- c(reaction.maps, reaction.map)
}

reaction.map.sort <- sort(unique(reaction.maps))

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
    rate_limit(2)
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
    rate_limit(2)
    num.fmt <- sprintf(sprintf.fmt, i)
    message("[", num.fmt, "] ", "Downloading: ", paste0(ten.comp.plus[[i]], collapse = ", "))
    
    tmp.compounds <- KEGGREST::keggGet(ten.comp.plus[[i]])
    compounds.flat <- c(compounds.flat, tmp.compounds)
}

message("Parsing compounds...")
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
                reaction = if (!is.null(compounds.flat[[x]]$REACTION)) {
                    list(unlist(strsplit(compounds.flat[[x]]$REACTION, "\\s+")))
                } else {
                    list(NULL)
                },
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

compounds.tib.filt <- compounds.tib %>%
    dplyr::filter(!grepl("[nR()]", formula))

compounds.tib.removed <- compounds.tib %>%
    dplyr::filter(grepl("[nR()]", formula))

message("Removed ", nrow(compounds.tib.removed), " rows with `n`, `R`, `(`, or `)`")

message("Combining reaction pairs tied to compounds derived from search terms...")
final.pairs <- filtered.rpairs %>%
    dplyr::filter(
        dplyr::if_any(
            .cols = all_of(c("comp1", "comp2")),
            .fns = ~ !.x %in% unname(compounds.tib.removed$entry)
        )
    ) %>%
    dplyr::left_join(
        x = .,
        y = compounds.tib.filt %>%
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
        y = compounds.tib.filt %>%
            dplyr::select(entry, formula, name, dblinks) %>%
            dplyr::rename(
                "formula2" = "formula",
                "name2" = "name",
                "dblink2" = "dblinks"
            ),
        by = c("comp2" = "entry")
    ) %>%
    tidyr::drop_na(formula1, formula2) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
        mass1 = Rdisop::getMass(Rdisop::getMolecule(formula1)),
        mass2 = Rdisop::getMass(Rdisop::getMolecule(formula2)),
        delta_mass = abs(mass1 - mass2)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(entry, rpair, comp1, comp2, formula1, formula2, mass1, mass2, delta_mass, name1, name2, dblink1, dblink2)

message("Writing results to rpairs.tsv...")
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
    file.path(outdir, "db_brite.info.txt"), ": information on brite database used\n",
    file.path(outdir, "rpairs.tsv"), ": dataframe with reaction pairs matched to search terms\n",
    "Search terms used from ", input, ": \n\t", paste0(data$query, collapse = "\n\t")
)

writeLines(text = output.info, con = file.path(outdir, "output_info.txt"))

message("Done.")
