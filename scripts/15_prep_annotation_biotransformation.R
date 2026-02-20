# ==============================================================================
# Prepare biotransformation checking -------------------------------------------
# ==============================================================================
message("Expanding possible adducts...")
xchr9_defs <- xcms::featureDefinitions(xchr9) %>%
  tibble::as_tibble(., rownames = "feature")

xchr9_mzs <- xchr9_defs$mzmed
names(xchr9_mzs) <- xchr9_defs$feature

possible_adducts <- MetaboCoreUtils::mz2mass(
  xchr9_mzs,
  adduct = adducts(polarity = opt$polarity)
) %>%
  tibble::as_tibble(., rownames = "feature") %>%
  tidyr::pivot_longer(
    cols = 2:ncol(.),
    names_to = "adduct",
    values_to = "mass"
  ) %>%
  dplyr::left_join(
    x = .,
    y = xchr9_defs,
    by = "feature"
  )

if (nrow(xchr9_defs) * 17 == nrow(possible_adducts)) {
  message("Adducts expanded")
} else {
  warning("Adducts did not correctly expanded")
}

message("Importing biotransformation file...")
bio_transf <- import_biotransform_meta(
  file = paste0(opt$data_path, "/", opt$biotransf_file)
)

if (!exists("biotransf_append")) {
  data <- readr::read_tsv(
    "scripts/search_compounds/output/rpairs.tsv",
    show_col_types = FALSE
  ) %>%
  dplyr::mutate(rpair_num = paste0("RP", dplyr::row_number()))

  df <- data %>%
    tidyr::drop_na(formula1, formula2) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      dplyr::across(
        .cols = dplyr::all_of(c("formula1", "formula2")),
        .fns = ~ MetaboCoreUtils::standardizeFormula(.)
      )
    ) %>%
    dplyr::mutate(
      delta_formula = {
        # Try subtracting both directions
        diff1 <- MetaboCoreUtils::subtractElements(formula1, formula2)
        diff2 <- MetaboCoreUtils::subtractElements(formula2, formula1)
        # Use whichever worked
        if (!is.na(diff1)) {
          diff1
        } else if (!is.na(diff2)) {
          diff2
        } else {
          paste0("\u00B1 ", formula1, " <=> ", formula2)
        }
      }
    ) %>%
    dplyr::mutate(delta_mass = Rdisop::getMonoisotopic(delta_formula)) %>%
    dplyr::ungroup() %>%
    dplyr::relocate(rpair_num, .before = "entry") %>%
    dplyr::mutate(
      across(
        .cols = tidyselect::all_of(c("name1", "name2")),
        .fns = ~ stringr::str_trim(stringr::str_extract(., "^[^;]+"))
      )
    ) %>%
    dplyr::mutate(allowed_n = 1) %>%
    tidyr::uncount(
      data = .,
      weights = allowed_n,
      .id = "multiplier",
      .remove = FALSE
    ) %>%
    dplyr::select(rpair_num, delta_formula, allowed_n, multiplier, delta_mass)

  biotransf_append <- df %>%
    dplyr::group_by(delta_mass) %>%
    dplyr::summarize(
      delta_formula = paste(sort(unique(delta_formula)), collapse = ", "),
      rpair_nums = paste(sort(unique(rpair_num)), collapse = ", "),
      n_rpairs   = dplyr::n_distinct(rpair_num),
      allowed_n  = dplyr::first(allowed_n),
      multiplier = dplyr::first(multiplier),
      .groups = "drop"
    ) %>%
    dplyr::select(
      name = rpair_nums,
      delta_formula,
      allowed_n,
      multiplier,
      delta_mass
    ) %>%
    dplyr::filter(!grepl("<=>", delta_formula)) %>%
    # nothing changed? or nothing has changed at least according
    # since the delta formulas probably only had 1n or 1c or similar and were
    # subtracted
    tidyr::drop_na(delta_mass)
} else {
  message("'biotransf_append' already exists.")
}

bio_transf2 <- dplyr::bind_rows(
  bio_transf,
  biotransf_append
)

# Potentially filter noisy features with the filt_features() function I made
all_sig_diff <- sort(unique(unlist(upset_comps, use.names = FALSE)))
possible_adducts_signif <- possible_adducts %>%
  dplyr::filter(feature %in% all_sig_diff)