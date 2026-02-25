cli::cli_h1(basename(this.path::this.path()))
# ==============================================================================
# Preparing to search for potential biotransformations -------------------------
# ==============================================================================
cli::cli_h3("Preparing to search for potential biotransformations")
cli::cli_alert_info("Expanding possible adducts")

xchr9_defs <- xcms::featureDefinitions(xchr9) %>%
  tibble::as_tibble(
    x =.,
    rownames = "feature",
    .name_repair = "universal"
  )

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
  cli::cli_alert_success("Adducts expanded")
} else {
  cli::cli_alert_danger("Adducts did not correctly expanded")
}

cli::cli_progress_step("Import biotransformation file")
bio_transf <- import_biotransform_meta(
  file = paste0(opt$data_path, "/", opt$biotransf_file)
)
cli::cli_progress_done()

if (file.exists(opt$rpairs_path)) {
  if (!exists("biotransf_append")) {
    cli::cli_alert_info(
      paste0(
        "'biotransf_append' doesn't exist, generating"
      )
    )
    cli::cli_progress_step("Parsing {.path {basename(opt$rpairs_path)}}")
    data <- readr::read_tsv(
      file = opt$rpairs_path,
      show_col_types = FALSE,
      progress = FALSE
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
        dplyr::across(
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
      cli::cli_progress_done()
  } else {
    cli::cli_alert_info("'biotransf_append' already exists, using it")
  }

  bio_transf2 <- dplyr::bind_rows(
    bio_transf,
    biotransf_append
  )
} else {
  cli::cli_alert_warning(
    paste0(
      "{.val {opt$rpairs_path} doesn't exist, skipping}"
    )
  )
  bio_transf2 <- bio_transf
}

all_sig_diff <- full_limma %>%
  dplyr::filter(adj.P.Val < opt$qvalue) %>%
  pull(feature) %>%
  unique() %>%
  sort()

possible_adducts_signif <- possible_adducts %>%
  dplyr::filter(feature %in% all_sig_diff)

filt_sig_diff <- xchr9_filt$feature[xchr9_filt$feature %in% all_sig_diff]
possible_adducts_filt <- possible_adducts %>%
  dplyr::filter(feature %in% filt_sig_diff)
