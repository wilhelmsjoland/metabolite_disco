
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