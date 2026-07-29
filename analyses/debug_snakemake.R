library(tidyverse)
options(width = 200)

exp_path <- paste0(
  "/Volumes/bluecub/aglycone_release_100um_24h/test"
)

snake_path <- file.path(exp_path, "snakemake_objects")
object_path <- file.path(exp_path, "tables")

list.files(snake_path, full.names = TRUE) %>%
  purrr::walk(
    .x = .,
    .f = ~ read_snakes(.x)
  )

source("scripts/18_biotransformer.R")

afz <- read_csv(file.path(object_path, "prediction_afzelin.csv"))
kaem <- read_csv(file.path(object_path, "prediction_kaempferol.csv"))

biot_dedup <- bind_rows(
  ... = list(
    dplyr::mutate(afz, input = "afzelin", .before = "InChI"),
    dplyr::mutate(kaem, input = "kaempferol", .before = "InChI")
  )
) %>%
  dplyr::group_by(InChIKey, input) %>%
  dplyr::summarize(
    dplyr::across(
      .cols = setdiff(colnames(.), c("InChIKey", "input")),
      .fns = ~ paste(unique(.x), collapse = ", ")
    ),
    .groups = "keep"
  )

biot_mass <- biot_dedup %>%
  dplyr::mutate(
    mass = MetaboCoreUtils::calculateMass(`Molecular formula`)
  ) %>%
  dplyr::relocate(mass, .before = "InChI")

biot_mets <- biot_mass$mass
names(biot_mets) <- biot_mass$InChIKey

# mass2mz() returns one row per input mass, in the same order - bind its
# columns directly onto biot_mass by position instead of round-tripping
# through a named vector, so no key (InChIKey, input, or otherwise) is
# needed to put the pieces back together afterward.
biot_mz <- MetaboCoreUtils::mass2mz(
  x = biot_mass$mass,
  adduct = MetaboCoreUtils::adducts(polarity = snakemake@params$polarity)
) %>%
  tibble::as_tibble()

biot_final <- dplyr::bind_cols(biot_mass, biot_mz) %>%
  tidyr::pivot_longer(
    cols = colnames(biot_mz),
    names_to = "adduct",
    values_to = "mz"
  ) %>%
  dplyr::arrange(InChIKey)

biot_mass_len <- length(biot_mass$InChIKey) *
  nrow(MetaboCoreUtils::adducts(polarity = snakemake@params$polarity))

if (biot_mass_len != nrow(biot_final)) {
  cli::cli_alert_danger(
    "The transformation prediction dataframes are not the same length."
  )
}

def_tib <- xcms::featureDefinitions(xchr9) %>%
  tibble::as_tibble(., rownames = "feature")
# biot_final = mass to mzs - > match the m/zs to the m/zs in the data
predicted_feats <- biot_final %>%
  dplyr::inner_join(
    x = .,
    y = def_tib %>%
      dplyr::mutate(
        tol = MsCoreUtils::ppm(mzmed, snakemake@params$ppm_match),
        mz_lo = mzmed - tol,
        mz_hi = mzmed + tol
      ) %>%
      dplyr::relocate(c("mz_lo", "mz_hi"), .after = "feature"),
    by = dplyr::join_by(dplyr::between(mz, mz_lo, mz_hi))
  )

predicted_feats
