source("analyses/chem_functions.R")

# ==============================================================================
# Plotting glycosides & aglycones ----------------------------------------------
# ==============================================================================
draw_out_path <- paste0(
  "/Users/wilhelm/Library/CloudStorage",
  "/ProtonDrive-wilhelm.sjoland@proton.me-folder",
  "/01_juniper/01_arbete/01_projekt/03_psm",
  "/molekyler/vektorbilder/glykosider_aglykoner"
)

glycone_list <- readxl::read_xlsx(
  path = paste0(
    "/Users/wilhelm/Proton/01_juniper/01_arbete/01_projekt/03_psm/",
    "molekyler/information/glycone_list.xlsx"
  ),
  sheet = "molecules",
  na = c("", "NA")
)

glycone_list_append <- tibble::tibble(
  molecule = c("Saligenin", "Flavan"),
  class = c("AG", "AG"),
  cayman_id = c(NA, NA),
  InChI = c(
    "InChI=1S/C7H8O2/c8-5-6-3-1-2-4-7(6)9/h1-4,8-9H,5H2",
    "InChI=1S/C15H14O/c1-2-6-12(7-3-1)15-11-10-13-8-4-5-9-14(13)16-15/h1-9,15H,10-11H2"
  ),
  note = c(NA, NA)
)

glycone_list <- glycone_list %>%
  dplyr::bind_rows(glycone_list_append)

for (i in seq_len(nrow(glycone_list))) {
  draw_from_inchi(
    inchi = glycone_list$InChI[i],
    name = glycone_list$molecule[i],
    output_path = file.path(
      draw_out_path,
      paste0(glycone_list$molecule[i], ".svg")
    )
  )
}


# ==============================================================================
# Prep plotting ----------------------------------------------------------------
# ==============================================================================
base_dir <- paste0(
  "/Volumes/bluecub/aglycone_release_100um_24h/output/experiment_analyses"
)
dirs <- c(
  "luteolin_7_o_glucuronide_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon",
  "luteolin_7_o_glucuronide_b_thetaiotaomicron_vpi_5482",
  "isorhamnetin_3_o_neohesperidoside_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon",
  "isorhoifolin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon"
)

for (dir in dirs) {
  dir.create(file.path(base_dir, dir, "predictions", "anno_sims"), FALSE, TRUE)
  dir.create(file.path(base_dir, dir, "predictions", "bio_sims"), FALSE, TRUE)

  # ============================================================================
  # Plotting temporary bio_sims matches ----------------------------------------
  # ============================================================================
  bio_sims_path <- file.path(
    base_dir,
    dir,
    "objects/extract_bio_sims.rds"
  )

  bio_sims <- readRDS(bio_sims_path)

  bio_sims_fix <- bio_sims %>%
    dplyr::mutate(
      name = paste0(round(sim, 2), "_", feature, "_", met_id),
      name = gsub("[^A-Za-z0-9._-]", "_", name)
    )
  bio_sims_out <- file.path(
    file.path(base_dir, dir, "predictions", "bio_sims")
  )

  for (i in seq_len(nrow(bio_sims_fix))) {
    draw_from_inchi(
      inchi = bio_sims_fix$InChI[i],
      name = bio_sims_fix$name[i],
      output_path = file.path(
        bio_sims_out,
        paste0(bio_sims_fix$name[i], ".svg")
      )
    )
  }

  # ============================================================================
  # Plotting temporary anno_sims matches ---------------------------------------
  # ============================================================================
  anno_sims_path <- file.path(
    base_dir,
    dir,
    "objects/extract_feats.rds"
  )
  anno_sims <- readRDS(anno_sims_path)
  anno_sims_fix <- anno_sims %>%
    dplyr::mutate(
      name = paste0(round(sim, 2), "_", peak_id, "_", target_name),
      name = gsub("[^A-Za-z0-9._-]", "_", name)
    )
  anno_sims_out <- file.path(
    file.path(base_dir, dir, "predictions", "anno_sims")
  )

  for (i in seq_len(nrow(anno_sims_fix))) {
    draw_from_inchi(
      inchi = anno_sims_fix$target_inchi[i],
      name = anno_sims_fix$name[i],
      output_path = file.path(
        anno_sims_out,
        paste0(anno_sims_fix$name[i], ".svg")
      )
    )
  }
}
