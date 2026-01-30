mkdir -p db_info
mkdir -p output

# # download database info used
# echo "Downloading database info"
# curl -s "https://rest.kegg.jp/info/rclass" -o db_info/kegg_rclass.txt
# curl -s "https://rest.kegg.jp/info/reaction" -o db_info/kegg_reaction.txt
# curl -s "https://rest.kegg.jp/info/compound" -o db_info/kegg_compound.txt
# 
# echo "Finding reactions matched to search terms..."
# N_ITERS=${#SEARCH_TERMS[@]}      # like length(compounds.split)
# ITER_CHARS=${#N_ITERS}             # like nchar(n.iters)
# readarray -t SEARCH_TERMS < search_terms_compounds.txt
# printf "entry\tdescription\tquery\n" > output/search_compounds.tsv
# for IDX in "${!SEARCH_TERMS[@]}"; do
#         NUM=$((IDX+1))
#         NUM_FMT=$(printf "%0*d" "${ITER_CHARS}" "${NUM}")
#         SEARCH_TERM="${SEARCH_TERMS[IDX]}"
#         echo "[${NUM_FMT}] Downloading hits for: ${SEARCH_TERM}"
#         curl -s "https://rest.kegg.jp/find/compound/${SEARCH_TERM}" |
#         awk -v term="${SEARCH_TERM}" 'BEGIN{FS=OFS="\t"} {print $0, term}' >> output/search_compounds.tsv
# done

# # Now take compounds and find reactions that are involved in this
# echo "Finding reactions matched to compounds..."
# COMPOUNDS=$(awk 'BEGIN{FS = "\t"; ORS = "\n"} NR > 1 {print $1}' output/search_compounds.tsv | xargs -n 10 | sed 's/ /+/g')
# readarray -t COMPOUND_ARR < <(printf '%s\n' "$COMPOUNDS")
# : > output/compound_mapping.txt
# N_ITERS=${#COMPOUND_ARR[@]}
# ITER_CHARS=${#N_ITERS}
# for IDX in "${!COMPOUND_ARR[@]}"; do
#         NUM=$((IDX+1))
#         NUM_FMT=$(printf "%0*d" "${ITER_CHARS}" "${NUM}")
# 	  COMPOUND="${COMPOUND_ARR[IDX]}"
#         echo "[${NUM_FMT}] Fetching mapping for: ${COMPOUND}"
#         curl -s "https://rest.kegg.jp/link/reaction/${COMPOUND}" >> output/compound_mapping.txt
# done
 
# awk 'BEGIN{FS = OFS = "\t"} {gsub(/rc:/, "", $2); print $2}' output/compound_mapping.txt | sort -u > output/reactions.tsv

# # Download kegg reaction based on reactions.tsv
# echo "Downloading reactions..."
# REACTIONS=$(awk 'BEGIN{FS = "\t"; ORS = "\n"} NR > 1 {print $1}' output/reactions.tsv | xargs -n 10 | sed 's/ /+/g')
# readarray -t REACTION_ARR < <(printf '%s\n' "$REACTIONS")
# : > output/kegg_reactions.txt
# N_ITERS=${#REACTION_ARR[@]}
# ITER_CHARS=${#N_ITERS}
# for IDX in "${!REACTION_ARR[@]}"; do
#         NUM=$((IDX+1))
#         NUM_FMT=$(printf "%0*d" "${ITER_CHARS}" "${NUM}")
#         REACTION="${REACTION_ARR[IDX]}"
#         echo "[${NUM_FMT}] Downloading: ${REACTION}"
#         curl -s "https://rest.kegg.jp/get/${REACTION}" >> output/kegg_reactions.txt
# done
 
# # Parse the kegg reaction into .tsv
# echo "Parsing reactions into tabular format"
# {
#   echo -e "entry\tname\tdefinition\tequation\tcomment\trclass\tenzyme\tpathway\tbrite\tdblinks\torthology\tmodule"
#   awk -v OFS='\t' '
#     function trim(s){ sub(/^[ \t]+/, "", s); sub(/[ \t]+$/, "", s); return s }
#     function norm(s){ gsub(/[ \t]+/, " ", s); return trim(s) }
#     function add(cur, val, sep){
#       val = norm(val)
#       if(val=="") return cur
#       return (cur=="" ? val : cur sep val)
#     }
# 
#     BEGIN {
#       key=""
#       entry=name=def=eq=comment=rclass=enzyme=pathway=brite=dblinks=orth=module=""
#     }
# 
#     $0=="///" {
#       print entry, name, def, eq, comment, rclass, enzyme, pathway, brite, dblinks, orth, module
#       key=""
#       entry=name=def=eq=comment=rclass=enzyme=pathway=brite=dblinks=orth=module=""
#       next
#     }
# 
#     {
#       k = trim(substr($0, 1, 12))
#       v = substr($0, 13)
# 
#       if (k != "") key = k   # new section; else continuation line keeps previous key
# 
#       if (key=="ENTRY") {
#         v = norm(v)
#         split(v, a, " ")
#         entry = a[1]
#       }
#       else if (key=="NAME")       name    = add(name,    v, " | ")
#       else if (key=="DEFINITION") def     = add(def,     v, " ")
#       else if (key=="EQUATION")   eq      = add(eq,      v, " ")
#       else if (key=="COMMENT")    comment = add(comment, v, " ")
#       else if (key=="RCLASS")     rclass  = add(rclass,  v, " | ")
#       else if (key=="ENZYME")     enzyme  = add(enzyme,  v, " ")
#       else if (key=="PATHWAY")    pathway = add(pathway, v, " | ")
#       else if (key=="BRITE")      brite   = add(brite,   v, " > ")
#       else if (key=="DBLINKS")    dblinks = add(dblinks, v, " | ")
#       else if (key=="ORTHOLOGY")  orth    = add(orth,    v, " | ")
#       else if (key=="MODULE")     module  = add(module,  v, " | ")
#     }
#   ' output/kegg_reactions.txt
# } > output/kegg_reactions.tsv

# echo "Finding compounds within reactions based on search terms for compounds"
# "/mnt/c/Program Files/R/R-4.5.2/bin/x64/Rscript.exe" -e '
# 	data <- read.delim(file.path("output", "kegg_reactions.tsv"))
# 	compounds <- unlist(regmatches(data$equation, gregexpr("\\bC\\d+\\b", data$equation, perl = TRUE)))
# 	writeLines(compounds, file.path("output", "compounds_to_filter.tsv"))
# '

# TODO FROM HERE

# TODO
# Here I need to find reaction classes from reactions
#

# TODO
# Here I need to download reaction classes mapped
#

# TODO
# Here I need to take the rclasses downloaded kegg information and fetch only pairs that are in my compound list


# # Final parse to create reaction_pairs
# echo "Parsing reaction classes to produce reaction pairs"
# "/mnt/c/Program Files/R/R-4.5.2/bin/x64/Rscript.exe" -e '
# 
# if (!file.exists(file.path(getwd(), "output", "rpairs.tsv"))) {
#     suppressWarnings(
#         suppressPackageStartupMessages({
#             library(tidyverse)
#             library(writexl)
#             library(readxl)
#         })
#     )
# 
#     rpair <- read_tsv(
#         file.path(getwd(), "output", "kegg_reaction_classes.tsv"),
#         show_col_types = FALSE
#     )
# 
#     rpair.filt <- rpair %>%
#         dplyr::mutate(rpair = strsplit(rpair, "\\s", perl = TRUE)) %>%
#         tidyr::unnest(cols = "rpair") %>%
#         dplyr::mutate(
#             comp1 = stringr::str_split_i(rpair, "_", 1),
#             comp2 = stringr::str_split_i(rpair, "_", 2)
#         ) %>%
#         dplyr::select(c(
#             "entry",
#             "rpair",
#             "comp1",
#             "comp2",
#             "definition",
#             "reaction",
#             "enzyme",
#             "pathway",
#             "orthology"
#         ))
# 
#         if (file.exists(file.path(getwd(), "output", "compounds.rds"))) {
#                 message("Reading compounds.rds from file...")
#                 comp.tib <- readRDS(file = file.path(getwd(), "output", "compounds.rds"))
#         } else {
#                 message("Fetching KEGG compounds in reaction pairs...")
#                 compounds <- unique(c(rpair.filt$comp1, rpair.filt$comp2))
#                 compounds.split <- split(compounds, ceiling(seq_along(compounds) / 10))
#                 comp.tib <- tibble::tibble()
# 
#                 n.iters <- length(compounds.split)
#                 iter.chars <- nchar(n.iters)
#                 sprintf.fmt <- paste0("%0", iter.chars, "d")
#                 for (i in seq_along(compounds.split)) {
#                         num.fmt <- sprintf(sprintf.fmt, i)
#                         message("[", num.fmt, "] ", "Fetching compounds: ", paste0(compounds.split[[i]], collapse = ", "))
#                         data <- KEGGREST::keggGet(dbentries = compounds.split[[i]])
# 
#                         tmp.tib <- dplyr::bind_rows(
#                         lapply(
#                                 X = data,
#                                 FUN = function(x) {
#                                         tib <- tibble::tibble(
#                                         entry = x$ENTRY,
#                                         name = paste0(x$NAME, collapse = " "),
#                                         formula = x$FORMULA,
#                                         dblinks = list(x$DBLINKS)
#                                 )
#                                 return(tib)
#                                         }
#                                 )
#                         )
#                 comp.tib <- bind_rows(comp.tib, tmp.tib)
#                 }
#         saveRDS(object = comp.tib, file = file.path(getwd(), "output", "compounds.rds"))
#         }
# 
#         final.tib <- comp.tib %>%
#                 dplyr::filter(!grepl("[nR()]", formula))
# 
#         removed <- comp.tib %>%
#                 dplyr::filter(grepl("[nR()]", formula))
# 
#     message("Removed ", nrow(removed), " rows with `n`, `R`, `(`, or `)`")
# 
#     reac.pair <- rpair.filt %>%
#     dplyr::left_join(
#         x = .,
#         y = final.tib %>%
#             dplyr::select(entry, formula, name) %>%
#             dplyr::rename(
#                 "formula1" = "formula",
#                 "comp_name1" = "name"
#                 ),
#         by = c("comp1" = "entry")
#     ) %>%
#     dplyr::left_join(
#         x = .,
#         y = final.tib %>%
#             dplyr::select(entry, formula, name) %>%
#             dplyr::rename(
#                 "formula2" = "formula",
#                 "comp_name2" = "name"
#                 ),
#         by = c("comp2" = "entry")
#     ) %>%
#     tidyr::drop_na(formula1, formula2) %>%
#     dplyr::rowwise() %>%
#     dplyr::mutate(
#         monomass1 = Rdisop::getMass(Rdisop::getMolecule(formula1)),
#         monomass2 = Rdisop::getMass(Rdisop::getMolecule(formula2))
#     ) %>%
#     dplyr::mutate(
#         delta_monomass = abs(monomass1 - monomass2)
#     ) %>%
#     dplyr::ungroup() %>%
#     dplyr::select(c(
#         "entry",
#         "rpair",
#         "comp1",
#         "comp2",
#         "comp_name1",
#         "comp_name2",
#         "formula1",
#         "formula2",
#         "monomass1",
#         "monomass2",
#         "delta_monomass",
#         "definition",
#         "reaction",
#         "enzyme",
#         "pathway",
#         "orthology"
#     ))
# 
#     message("Saving reaction pairs to .tsv...")
#     readr::write_tsv(
#         x = reac.pair,
#         file = file.path(getwd(), "output", "rpairs.tsv")
#     )
# 
#     message("Done.")
# 
#     }        else {
#         message("rpairs.tsv already exists.")
#     }
# '
