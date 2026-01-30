findSimilarity <- function(
        smiles.obj = NULL,
        query.smiles = NULL,
        anno.obj = NULL,
        type = "circular",
        method = "tanimoto"
        ) {
    
    target.mols <- parse.smiles(smiles = smiles.obj, omit.nulls = TRUE)
    query.mol <- parse.smiles(query.smiles)[[1]]
    query.fp <- get.fingerprint(query.mol, type = type)
    target.fps <- lapply(target.mols, get.fingerprint, type = type)
    
    sims <- data.frame(
        sim = do.call(rbind, 
                      lapply(
                          target.fps,
                          fingerprint::distance,
                          fp2 = query.fp, 
                          method = method
                      )
        ))
    
    sims <- tibble::as_tibble(sims, rownames = "smiles") %>%
        dplyr::arrange(desc(sim)) %>%
        dplyr::left_join(
            x = .,
            y = anno.obj,
            by = c("smiles" = "feature")
        )
    
    return(sims)
}