library(janitor)
library(yaml)
library(writexl)

exp_path <- "/Volumes/bluecub/aglycone_release_100um_24h"
run_path <- "/Users/wilhelm/Proton/01_juniper/01_arbete/01_projekt/03_psm/"

# ==============================================================================
# Create experimental metadata from Longs experiments --------------------------
# ==============================================================================
exp_inds <- readr::read_csv(
  file = file.path(
    exp_path,
    "data",
    "experiment",
    "info",
    "fixed_070825_index.csv"
  ),
  show_col_types = FALSE,
  progress = TRUE,
  name_repair = "universal"
) %>%
  dplyr::rename("sample" = "Data.File") %>%
  dplyr::mutate(sample = basename(sample)) %>%
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
      "sample",
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
    "experiment",
    "mzml"
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

if (!all(exp_inds$sample %in% all_mzml)) {
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
) %>%
  dplyr::mutate(
    experiment_id = dplyr::row_number(),
    name = paste0(
      condition, "_", purrr::map_chr(strains, paste, collapse = "_and_")
    )
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
  setNames(combos$name) %>%
  dplyr::bind_rows(.id = "experiment") %>%
  dplyr::left_join(
    dplyr::select(combos, experiment_id, name),
    by = c("experiment" = "name")
  ) %>%
  dplyr::arrange(experiment_id, strain) %>%
  dplyr::mutate(group = paste0(strain, "_", condition)) %>%
  dplyr::select(
    c(
      "sample",
      "group",
      "experiment_id",
      "experiment",
      "strain",
      "condition",
      "unclean_strain",
      "unclean_condition"
    )
  )

# ==============================================================================
# Get information on molecules in experiments ----------------------------------
# ==============================================================================
glycone_list <- readxl::read_xlsx(
  path = file.path(
    run_path,
    "molekyler",
    "information",
    "aglycone_release_100um_24h",
    "glycone_list.xlsx"
  ),
  sheet = "molecules",
  na = c("", "NA")
)

inchis <- glycone_list$InChI[!is.na(glycone_list$InChI)]

pubchem_info_path <- file.path(
  run_path,
  "molekyler",
  "information",
  "aglycone_release_100um_24h",
  "pubchem_info.csv"
)

if (!file.exists(pubchem_info_path)) {
  # httr2's req_perform() errors on a non-2xx status, so a PubChem error page
  # can never reach read_csv(). httr::POST returned it as an ordinary response,
  # and the HTML got parsed as data - that is how a 502 silently turned into
  # eight all-NA rows with an "<html><head>" column, and how a molecule lost
  # its SMILES.
  pubchem_req <- httr2::request(
    paste0(
      "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/",
      "property/Title,MolecularFormula,InChIKey,InChI",
      ",ExactMass,XLogP,TPSA,IsomericSMILES/CSV"
    )
  ) %>%
    # PubChem PUG REST usage policy: no more than 5 requests/second
    httr2::req_throttle(capacity = 5, fill_time_s = 1) %>%
    httr2::req_retry(
      max_tries = 5,
      is_transient = \(resp) httr2::resp_status(resp) %in% c(
        429,
        502,
        503,
        504
      ),
      backoff = \(i) 2^i
    )

  pubchem_info <- inchis %>%
    purrr::map(
      .f = ~ {
        resp <- pubchem_req %>%
          httr2::req_body_form(inchi = .x) %>%
          httr2::req_perform()

        readr::read_csv(
          I(httr2::resp_body_string(resp)),
          show_col_types = FALSE,
          progress = FALSE
        )
      }
    ) %>%
    dplyr::bind_rows()

  readr::write_csv(
    x = pubchem_info,
    file = pubchem_info_path,
  )
} else {
  pubchem_info <- readr::read_csv(
    file = pubchem_info_path,
    show_col_types = FALSE,
    progress = FALSE
  )
}

glycone_info <- glycone_list %>%
  dplyr::left_join(
    x = .,
    y = pubchem_info,
    by = "InChI"
)

# ==============================================================================
# Get chemical formula and SMILES for glycosides & aglycones -------------------
# ==============================================================================
glycone_pairs <- readxl::read_xlsx(
  path = file.path(
    run_path,
    "molekyler",
    "information",
    "aglycone_release_100um_24h",
    "glycone_list.xlsx"
  ),
  sheet = "glycone_pairs",
  na = c("", "NA")
)

# Make pairs for metadata
glycone_pairs_metadata <- glycone_pairs %>%
  dplyr::rename(
    "glycoside_unclean" = "glycoside",
    "aglycone_unclean" = "aglycone"
  ) %>%
  dplyr::mutate(
    dplyr::across(
      .cols = c("glycoside_unclean", "aglycone_unclean"),
      .fns = ~ janitor::make_clean_names(
        string = .,
        case = "snake",
        allow_dupes = TRUE
      ),
      .names = "{stringr::str_remove_all(.col, '_unclean')}"
    )
  ) %>%
  dplyr::left_join(
    x = .,
    y = glycone_info %>%
      dplyr::select(
        c(
          "cayman_id",
          "InChI",
          "SMILES",
          "MolecularFormula"
        )
      ),
    by = c("glycoside_cayman_id" = "cayman_id")
  ) %>%
  dplyr::rename(
    "glycoside_InChI" = "InChI",
    "glycoside_SMILES" = "SMILES",
    "glycoside_form" = "MolecularFormula"
  ) %>%
  dplyr::left_join(
    x = .,
    y = glycone_info %>%
      dplyr::select(
        c(
          "cayman_id",
          "InChI",
          "SMILES",
          "MolecularFormula"
        )
      ),
    by = c("aglycone_cayman_id" = "cayman_id")
  ) %>%
  dplyr::rename(
    "aglycone_InChI" = "InChI",
    "aglycone_SMILES" = "SMILES",
    "aglycone_form" = "MolecularFormula"
  ) %>%
  dplyr::select(
    c(
      "glycoside",
      "aglycone",
      "glycoside_unclean",
      "aglycone_unclean",
      "glycoside_form",
      "aglycone_form",
      "glycoside_SMILES",
      "aglycone_SMILES"
    )
  )

# Create a list of experiments to inspect later
experiment_matrix_path <- file.path(
  run_path,
  "molekyler",
  "information",
  "aglycone_release_100um_24h",
  "experiment_matrix.xlsx"
)

if (!file.exists(experiment_matrix_path)) {
  writexl::write_xlsx(
    x = experiments,
    path = experiment_matrix_path
  )
}

# Create a list of glycones to inspect later
glycones <- readxl::read_xlsx(
  path = file.path(
    run_path,
    "molekyler",
    "information",
    "aglycone_release_100um_24h",
    "glycone_list.xlsx"
  ),
  sheet = "glycone_pairs"
) %>%
  dplyr::select(glycoside, aglycone)

cln_glycone_nms <- glycones %>%
  dplyr::mutate(
    dplyr::across(
      .cols = c("glycoside", "aglycone"),
      .fns = ~ janitor::make_clean_names(
        string = .,
        case = "snake",
        allow_dupes = TRUE
      )
    )
  )

cln_glycone_nms_path <- file.path(
  run_path,
  "molekyler",
  "information",
  "aglycone_release_100um_24h",
  "clean_names_glycone_list.xlsx"
)
if (!file.exists(cln_glycone_nms_path)) {
  writexl::write_xlsx(
    x = cln_glycone_nms,
    path = cln_glycone_nms_path
  )
}


clean_glycone_metadata <- glycone_info %>%
  dplyr::mutate(
    molecule = janitor::make_clean_names(
      string = molecule,
      case = "snake",
      allow_dupes = TRUE
    )
  )

clean_glycone_metadata_path <- file.path(
  run_path,
  "molekyler",
  "information",
  "aglycone_release_100um_24h",
  "clean_glycone_metadata.xlsx"
)

if (!file.exists(clean_glycone_metadata_path)) {
  writexl::write_xlsx(
    x = clean_glycone_metadata,
    path = clean_glycone_metadata_path
  )
}
