mol_similarity <- function(
  query_smiles = NULL,
  target_smiles = NULL,
  kekulise = TRUE,
  omit_nulls = TRUE,
  fingerprint = "circular",
  circular_type = "ECFP6",
  method = "tanimoto"
) {

  query_mol <- rcdk::parse.smiles(
    smiles = query_smiles,
    kekulise = kekulise,
    omit.nulls = omit_nulls
  )[[1]]

  target_mols <- rcdk::parse.smiles(
    smiles = target_smiles,
    kekulise = kekulise,
    omit.nulls = omit_nulls
  )

  query_fp <- rcdk::get.fingerprint(
    molecule = query_mol,
    type = fingerprint,
    circular.type = circular_type
  )

  target_fps <- lapply(
    X = target_mols,
    FUN = function(x) {
      rcdk::get.fingerprint(
        molecule = x,
        type = fingerprint,
        circular.type = circular_type
      )
    }
  )

  sims <- target_fps %>%
    purrr::map2(
      .x = .,
      .y = names(.),
      .f = ~ tibble::tibble(
        feature = .y,
        sim = fingerprint::distance(
          fp1 = .x,
          fp2 = query_fp,
          method = "tanimoto"
        )
      )
    ) %>%
    purrr::list_rbind() %>%
    dplyr::arrange(dplyr::desc(sim))

  return(sims)
}

# TODO
# check later -> hclust of similarity
# fps <- lapply(mols, get.fingerprint, type='circular')
# fp.sim <- fingerprint::fp.sim.matrix(fps, method='tanimoto')
# fp.dist <- 1 - fp.sim
# cls <- hclust(as.dist(fp.dist))
# plot(cls, labels=FALSE)

run_biotransformer <- function(
  bt_dir = "biotransformer3.0jar",
  smiles = "C1=CC(=CC=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)O",
  b_type = "superbio",
  k_task = "pred", # pred for prediction, or cid for compound identification
  output_file = "apigenin"
) {
  old_wd <- getwd()
  new_wd <- file.path(old_wd, bt_dir)
  biot_output_loc <- file.path(old_wd, res_folder, "tables")
  clean_nm <- gsub("\\..*$", "", output_file)

  cmd <- paste(
    "java -jar BioTransformer3.0_20230525.jar",
    "-k", k_task,
    "-b", b_type,
    # "-isdf", biot_output_loc,
    "-ismi", smiles,
    "-ocsv", paste0(biot_output_loc, "/", clean_nm, ".csv"),
    "-a"
  )

  setwd(new_wd)
  biot_output <- system(cmd, intern = TRUE)
  setwd(old_wd)

  log_file <- file.path(biot_output_loc, paste0(clean_nm, ".log"))
  writeLines(biot_output, log_file)
}