library(pacman)
p_load(tidyverse)
p_load(rcdk)
p_load(httr2)
p_load(readxl)
source("analyses/analyses_functions.R")

feature_map <- readxl::read_xlsx(
  path = paste0(
    "/Users/wilhelm/Proton/01_juniper/01_arbete/01_projekt/03_psm/mapping/",
    "mapping_file.xlsx"
  )
) %>%
  dplyr::mutate(
    dplyr::across(
      .cols = dplyr::everything(),
      .fns = ~ stringr::str_remove_all(., ".pdf")
    )
  )

colnames(feature_map) %>%
  purrr::walk(
    .x = .,
    .f = ~ {
      output_path <- paste0(
        "/Volumes/bluecub/aglycone_release_100um_24h",
        "/output/experiment_analyses/",
        .x
      )

      output_files <- list.files(
        path = file.path(output_path, "objects"),
        full.names = TRUE,
        pattern = ".rds"
      )

      expe <- output_files %>%
        set_names(., nm = stringr::str_remove_all(basename(.), ".rds")) %>%
        purrr::map(
          .x = .,
          .f = ~ readRDS(file = .x)
        )

      good_peaks <- as.vector(na.omit(feature_map[[.x]]))
      message(
        "experiment: ", .x,
        " | peaks: ", length(good_peaks),
        " | expe keys: ", paste(names(expe), collapse = ", ")
      )

      good_peaks %>%
        purrr::walk(
          .x = .,
          .f = ~ {
            bio_sims_smiles <- NULL
            extract_feats_smiles <- NULL

            if (!is.null(expe[["extract_bio_sims"]])) {
              bio_sims_smiles <- expe[["extract_bio_sims"]] %>%
                dplyr::filter(feature == .x) %>%
                dplyr::mutate(sim = round(sim, 2)) %>%
                dplyr::mutate(target_name = NA_character_) %>%
                dplyr::select(
                  feature = feature,
                  sim = sim,
                  target_name = target_name,
                  target_smiles = SMILES
                )
            }
            if (!is.null(expe[["extract_feats"]])) {
              extract_feats_smiles <- expe[["extract_feats"]] %>%
                dplyr::filter(peak_id == .x) %>%
                dplyr::mutate(sim = round(sim, 2)) %>%
                dplyr::select(
                  feature = peak_id,
                  sim = sim,
                  target_name = target_name,
                  target_smiles = target_smiles
                )
            }

            smiles_tib <- dplyr::bind_rows(
              bio_sims_smiles,
              extract_feats_smiles
            )

            message(
              "Peak: ", .x,
              " | bio_sims: ", nrow(bio_sims_smiles),
              " | feats: ", nrow(extract_feats_smiles))

            if (nrow(smiles_tib) == 0) {
              return(invisible(NULL))
            }

            smiles_full <- smiles_tib %>%
              dplyr::mutate(
                row_no = dplyr::row_number(),
                title = dplyr::if_else(
                  is.na(target_name),
                  paste0("Similarity: ", sim),
                  paste0(
                    target_name,
                    ", Sim = ",
                    sim
                  )
                ),
                file = paste0(
                  output_path, "/graphs/", .x, "/", .x, "_", row_no, ".svg"
                )
              )

            dir.create(
              file.path(output_path, "graphs", .x),
              recursive = TRUE,
              showWarnings = FALSE
            )

            purrr::pmap(
              .l = list(
                smiles = smiles_full$target_smiles,
                file   = smiles_full$file,
                title  = smiles_full$title
              ),
              .f = smiles_to_svg
            )
          }
        )
    }
  )


# int_mols %>%
#   purrr::walk(
#     .x = .,
#     .f = ~ {
#       bio_sims_smiles <- experiment$extract_bio_sims %>%
#         dplyr::filter(feature == .x) %>%
#         dplyr::mutate(sim = round(sim, 2)) %>%
#         dplyr::mutate(target_name = NA) %>%
#         dplyr::select(
#           feature = feature,
#           sim = sim,
#           target_name = target_name,
#           target_smiles = SMILES
#         )

#       extract_feats_smiles <- experiment$extract_feats %>%
#         dplyr::filter(peak_id == .x) %>%
#         dplyr::mutate(sim = round(sim, 2)) %>%
#         dplyr::select(
#           feature = peak_id,
#           sim = sim,
#           target_name = target_name,
#           target_smiles = target_smiles
#         )

#       dir.create(
#         file.path(output_path, "/graphs/", .x),
#         recursive = TRUE,
#         showWarnings = FALSE
#       )
#       smiles_tib <- dplyr::bind_rows(bio_sims_smiles, extract_feats_smiles) %>%
#         dplyr::mutate(
#           row_no = dplyr::row_number(),
#           title = dplyr::if_else(
#             is.na(target_name),
#             paste0("Similarity: ", sim),
#             paste0(
#               target_name,
#               ", Sim = ",
#               sim
#             )
#           ),
#           file = paste0(
#             output_path, "/graphs/", .x, "/", .x, "_", row_no, ".svg"
#           )
#         )

#       purrr::pmap(
#         .l = list(
#           smiles = smiles_tib$target_smiles,
#           file   = smiles_tib$file,
#           title  = smiles_tib$title
#         ),
#         .f = smiles_to_svg
#       )
#     }
#   )
