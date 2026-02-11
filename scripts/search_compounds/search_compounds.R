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
    name = vapply(hits, `[`, "", 3)
  )
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
cmd_call <- commandArgs()

input_idx <- match(TRUE, args %in% c("-i", "--input"))
outdir_idx <- match(TRUE, args %in% c("-o", "--output"))
append_idx <- match(TRUE, args %in% c("-a", "--append"))

input <- args[input_idx + 1]
outdir <- args[outdir_idx + 1]
append <- args[append_idx + 1]

input  <- if (!is.na(input_idx)  && input_idx  < length(args)) {
  args[input_idx  + 1]
} else {
  NA
}
outdir <- if (!is.na(outdir_idx) && outdir_idx < length(args)) {
  args[outdir_idx + 1]
} else {
  NA
}
append <- if (!is.na(append_idx) && append_idx < length(args)) {
  args[append_idx + 1]
} else {
  NA
}

usable_flags <- c("-i", "--input", "-o", "--output", "-a", "--append")

run_flag <- (length(args) == 4 || length(args) == 6) &&
  (!is.na(input_idx) && !is.na(outdir_idx)) &&
  (input_idx < length(args) && outdir_idx < length(args)) &&
  (is.na(append_idx) || append_idx < length(args)) &&
  (!input %in% usable_flags) &&
  (!outdir %in% usable_flags) &&
  (is.na(append_idx) || !append %in% usable_flags)

if (isFALSE(run_flag)) {
  stop(
    "\nUsage: 'Rscript -i <input> -o <output> -a <brite_categories>'",
    "\nYou need to supply either two or three arguments:\n\t",
    "Required: input with either '-i', or '--input'\n\t",
    "Required: output with either '-o' or '--output'\n\t",
    "Optional: append with either '-a' or '--append'\n\t",
    "Example: 'Rscript search_compounds.R -i",
    "comps.tsv -o out -a br:br08021,br:br08005'"
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
time_run <- format(Sys.time(), usetz = TRUE)

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
n_iters <- length(data$query)
iter_chars <- nchar(n_iters)
sprintf_fmt <- paste0("%0", iter_chars, "d")
compounds <- c()
for (i in seq_along(data$query)) {
  rate_limit(2)
  num_fmt <- sprintf(sprintf_fmt, i)
  message(
    "[", num_fmt, "] ", "Fetching compounds related to: ",
    paste0(data$query[[i]], collapse = ", ")
  )

  compound <- KEGGREST::keggFind("compound", data$query[[i]])
  compounds <- c(compounds, compound)
}

compounds_sort <- sort(unique(names(compounds)))

if (length(args) == 6) {
  brites <- strsplit(args[6], ",")[[1]]

  message("Downloading supplied brite categories...")
  x_all <- c()
  for (i in brites) {
    rate_limit(2)
    # Fix the formatting here
    x <- KEGGREST::keggGet(i)[[1]]
    x_all <- c(x_all, x)
  }

  message("Parsing supplied brite categories")
  brite_compounds <- dplyr::bind_rows(
    lapply(
      X = x_all,
      FUN = function(x) {
        extract_brite(x)
      }
    )
  )

  message("Appending compounds to search term compounds...")
  compounds_sort <- sort(
    unique(
      c(
        compounds_sort,
        paste0("cpd:", brite_compounds$compound)
      )
    )
  )
}

message("Finding mapping between compounds -> reactions -> reaction classes...")

ten_compound_list <- split(
  compounds_sort,
  ceiling(seq_along(compounds_sort) / 10)
)
ten_compound_plus <- lapply(
  X = ten_compound_list,
  FUN = function(x) {
    gsub("rc:", "", paste0(x, collapse = "+"))
  }
)

n_iters <- length(ten_compound_plus)
iter_chars <- nchar(n_iters)
sprintf_fmt <- paste0("%0", iter_chars, "d")
compound_maps <- c()
for (i in seq_along(ten_compound_plus)) {
  rate_limit(2)
  num_fmt <- sprintf(sprintf_fmt, i)
  message(
    "[", num_fmt, "] ", "Mapping: ",
    paste0(ten_compound_plus[[i]], collapse = ", ")
  )

  compound_map <- KEGGREST::keggLink("reaction", ten_compound_plus[[i]])
  compound_maps <- c(compound_maps, compound_map)
}

sort_compound_maps <- sort(unique(compound_maps))

ten_reaction_list <- split(
  sort_compound_maps,
  ceiling(seq_along(sort_compound_maps) / 10)
)
ten_reaction_plus <- lapply(
  X = ten_reaction_list,
  FUN = function(x) {
    gsub("rc:", "", paste0(x, collapse = "+"))
  }
)

n_iters <- length(ten_reaction_plus)
iter_chars <- nchar(n_iters)
sprintf_fmt <- paste0("%0", iter_chars, "d")
reaction_maps <- c()
for (i in seq_along(ten_reaction_plus)) {
  rate_limit(2)
  num_fmt <- sprintf(sprintf_fmt, i)
  message(
    "[", num_fmt, "] ", "Mapping: ",
    paste0(ten_reaction_plus[[i]], collapse = ", ")
  )

  reaction_map <- KEGGREST::keggLink("rclass", ten_reaction_plus[[i]])
  reaction_maps <- c(reaction_maps, reaction_map)
}

reaction_map_sort <- sort(unique(reaction_maps))

ten_rclass_list <- split(
  reaction_map_sort,
  ceiling(seq_along(reaction_map_sort) / 10)
)

ten_rclass_plus <- lapply(
  X = ten_rclass_list,
  FUN = function(x) {
    gsub("rc:", "", paste0(x, collapse = "+"))
  }
)

message("Downloading reaction classes...")
n_iters <- length(ten_rclass_plus)
iter_chars <- nchar(n_iters)
sprintf_fmt <- paste0("%0", iter_chars, "d")
rclasses <- list()
for (i in seq_along(ten_rclass_plus)) {
  rate_limit(2)
  num_fmt <- sprintf(sprintf_fmt, i)
  message(
    "[", num_fmt, "] ", "Downloading: ",
    paste0(ten_rclass_plus[[i]], collapse = ", ")
  )

  rclass <- KEGGREST::keggGet(ten_rclass_plus[[i]])
  rclasses <- c(rclasses, rclass)
}

message("Parsing reaction classes...")
rclass_tib <- dplyr::bind_rows(
  lapply(
    X = seq_along(rclasses),
    FUN = function(x) {
      tibble::tibble(
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
filtered_rpairs <- rclass_tib %>%
  dplyr::select(c("entry", "rpair")) %>%
  tidyr::unnest(rpair) %>%
  dplyr::mutate(
    comp1 = stringr::str_split_i(rpair, "_", 1),
    comp2 = stringr::str_split_i(rpair, "_", 2)
  ) %>%
  dplyr::filter(
    dplyr::if_any(
      .cols = c("comp1", "comp2"),
      .fns = ~ .x %in% gsub("cpd:", "", compounds_sort)
    )
  )

# Why is this shorter?
# Probably because at least one of the compounds may not be in a
# specific reaction -> I need reactions where the molecule of interest in
# not necessarily in the name of both compounds in the reaction
filt_comps <- sort(unique(c(filtered_rpairs$comp1, filtered_rpairs$comp2)))
filt_comps_orig <- sort(unique(gsub("cpd:", "", compounds_sort)))

ten_comp_list <- split(filt_comps, ceiling(seq_along(filt_comps) / 10))
ten_comp_plus <- lapply(
  X = ten_comp_list,
  FUN = function(x) {
    paste0(x, collapse = "+")
  }
)

message(
  "Downloading compunds part of reaction pairs, ",
  "based on original search terms..."
)
n_iters <- length(ten_comp_plus)
iter_chars <- nchar(n_iters)
sprintf_fmt <- paste0("%0", iter_chars, "d")
compounds_flat <- list()
for (i in seq_along(ten_comp_plus)) {
  rate_limit(2)
  num_fmt <- sprintf(sprintf_fmt, i)
  message(
    "[", num_fmt, "] ", "Downloading: ",
    paste0(ten_comp_plus[[i]], collapse = ", ")
  )

  tmp_compounds <- KEGGREST::keggGet(ten_comp_plus[[i]])
  compounds_flat <- c(compounds_flat, tmp_compounds)
}

message("Parsing compounds...")
compounds_tib <- dplyr::bind_rows(
  lapply(
    X = seq_along(compounds_flat),
    FUN = function(x) {
      tibble::tibble(
        entry = compounds_flat[[x]]$ENTRY,
        name = paste0(gsub(";", "", compounds_flat[[x]]$NAME), collapse = "; "),
        formula = compounds_flat[[x]]$FORMULA,
        exact_mass = compounds_flat[[x]]$EXACT_MASS,
        mol_weight = compounds_flat[[x]]$MOL_WEIGHT,
        reaction = if (!is.null(compounds_flat[[x]]$REACTION)) {
          list(unlist(strsplit(compounds_flat[[x]]$REACTION, "\\s+")))
        } else {
          list(NULL)
        },
        pathway = list(compounds_flat[[x]]$PATHWAY),
        enzyme = list(compounds_flat[[x]]$ENZYME),
        brite = list(compounds_flat[[x]]$BRITE),
        dblinks = paste0(compounds_flat[[x]]$DBLINK, collapse = "; "),
        atom = list(compounds_flat[[x]]$ATOM),
        bond = list(compounds_flat[[x]]$BOND)
      )
    }
  )
)

compounds_tib_filt <- compounds_tib %>%
  dplyr::filter(!grepl("[nR()]", formula))

compounds_tib_removed <- compounds_tib %>%
  dplyr::filter(grepl("[nR()]", formula))

message(
  "Removed ", nrow(compounds_tib_removed),
  " rows with `n`, `R`, `(`, or `)`"
)

message(
  "Combining reaction pairs tied to compounds ",
  "derived from search terms..."
)
final_pairs <- filtered_rpairs %>%
  dplyr::filter(
    dplyr::if_any(
      .cols = all_of(c("comp1", "comp2")),
      .fns = ~ !.x %in% unname(compounds_tib_removed$entry)
    )
  ) %>%
  dplyr::left_join(
    x = .,
    y = compounds_tib_filt %>%
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
    y = compounds_tib_filt %>%
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
  dplyr::select(
    entry, rpair, comp1, comp2, formula1, formula2,
    mass1, mass2, delta_mass, name1, name2, dblink1, dblink2
  )

message("Writing results to rpairs.tsv...")
write.table(
  x = final_pairs,
  file = file.path(outdir, "rpairs.tsv"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

output_info <- paste0(
  "Call: ", paste0(cmd_call, collapse = " "), "\n",
  "search_compounds.sh run on: ", time_run, "\n",
  file.path(outdir, "db_compound.info.txt"),
  ": information on compound database used\n",
  file.path(outdir, "db_reaction.info.txt"),
  ": information on reaction database used\n",
  file.path(outdir, "db_rclass.info.txt"),
  ": information on reaction class database used\n",
  file.path(outdir, "db_brite.info.txt"),
  ": information on brite database used\n",
  file.path(outdir, "rpairs.tsv"),
  ": dataframe with reaction pairs matched to search terms\n",
  "Search terms used from ", input, ": \n\t",
  paste0(data$query, collapse = "\n\t")
)

writeLines(text = output_info, con = file.path(outdir, "output_info.txt"))

message("Done.")
