data <- readr::read_tsv(
  "parse_kegg/search_compounds/output/rpairs.tsv",
  show_col_types = FALSE
)

df <- data %>%
  # dplyr::select(entry, rpair, comp1, comp2, name1, name2, delta_mass) %>%
  # dplyr::distinct(delta_mass, .keep_all = TRUE) %>%
  dplyr::mutate(rpair_num = paste0("RP", dplyr::row_number())) %>%
  dplyr::relocate(rpair_num, .before = "entry") %>%
  dplyr::mutate(
    across(
      .cols = all_of(c("name1", "name2")),
      .fns = ~ stringr::str_trim(stringr::str_extract(., "^[^;]+"))
    )
  ) %>%
  dplyr::mutate(delta_formula = dplyr::if_else(
    mass1 >= mass2,
    # These don't work for all formulas
    MetaboCoreUtils::subtractElements(formula1, formula2),
    MetaboCoreUtils::subtractElements(formula2, formula1)
  ))

biotransf.append.redundant <- df %>%
  dplyr::mutate(allowed_n = 1) %>%
  tidyr::uncount(
    data = ., 
    weights = allowed_n,
    .id = "multiplier",
    .remove = FALSE
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(delta_formula = {
    if (!is.na(delta_formula) && delta_mass != 0) {
      multChem(delta_formula, multiplier)
    } else {
      NA
    }
  }) %>%
  dplyr::mutate(delta_formula = paste0("± ", delta_formula)) %>%
  dplyr::ungroup() %>%
  dplyr::select(rpair_num, delta_formula, allowed_n, multiplier, delta_mass)

biotransf.append <- biotransf.append.redundant %>%
  dplyr::group_by(delta_formula, delta_mass) %>%
  dplyr::summarize(
    rpair_nums = paste(sort(unique(rpair_num)), collapse = ", "),
    n_rpairs   = dplyr::n_distinct(rpair_num),
    allowed_n  = dplyr::first(allowed_n),
    multiplier = dplyr::first(multiplier),
    .groups = "drop"
  ) %>%
  dplyr::select(
    name = rpair_nums,
    chem_formula = delta_formula,
    allowed_n,
    multiplier,
    delta_mass
  )
