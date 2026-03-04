library(janitor)

exp_path <- "V:/aglycone_release_100um_24h"

exp_inds <- readr::read_csv(
  file = file.path(exp_path, "data", "070825_index.csv"),
  show_col_types = FALSE,
  progress = TRUE,
  name_repair = "universal"
) %>%
  dplyr::rename("path" = "Data.File") %>%
  dplyr::mutate(path = basename(path)) %>%
  dplyr::rename(
    "unclean_condition" = "condition",
    "unclean_strain" = "strain"
  ) %>%
  dplyr::mutate(
    condition = janitor::make_clean_names(
      string = unclean_condition,
      case = "snake",
      allow_dupes = TRUE
    ),
    strain = janitor::make_clean_names(
      string = unclean_strain,
      case = "snake",
      allow_dupes = TRUE
    )
  ) %>%
  dplyr::select(
    c(
      "path",
      "strain",
      "condition",
      "unclean_strain",
      "unclean_condition"
    )
  ) %>%
  dplyr::group_by(strain, condition)

all_mzml <- list.files(
  path = file.path(
    exp_path,
    "data",
    "mzml_files"
  ),
  pattern = ".mzML"
)

# Can left_join() this to get the original names
name_lookup <- exp_inds %>%
  dplyr::distinct(
    strain,
    condition,
    unclean_strain,
    unclean_condition
  ) %>%
  dplyr::ungroup() %>%
  tidyr::unite(
    col = "name",
    strain,
    condition,
    remove = FALSE
  ) %>%
  tidyr::unite(
    col = "unclean_name",
    unclean_strain,
    unclean_condition,
    sep = " — ",
    remove = FALSE
  )

if (!all(exp_inds$path %in% all_mzml)) {
  cli::cli_abort("All mzml files don't exist")
}

exp_groups <- list(
  c("b_uniformis_atcc_8492", "bu_gsh_d_ggh_c_gsh_g"),
  c("b_ovatus_atcc_8483", "b_ovatus_atcc_8483_d_operon"),
  "b_thetaiotaomicron_vpi_5482",
  "p_copri_i_ak263"
)
conditions <- unique(exp_inds$condition[exp_inds$condition != "ycfa_glucose"])

combos <- tidyr::expand_grid(
  strains = exp_groups,
  condition = conditions
)

experiments <- purrr::map2(
  .x = combos$strains,
  .y = combos$condition,
  .f = ~ {
    exp_inds %>%
      dplyr::filter(
        strain %in% .x,
        condition %in% c(.y, "ycfa_glucose")
      )
  }
) %>%
  dplyr::bind_rows(.id = "experiment") %>%
  dplyr::arrange(experiment, strain)
