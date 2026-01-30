query.mol <- parse.smiles("C1[C@@]([C@H]([C@@H](O1)O[C@@H]2[C@H]([C@@H]([C@H](O[C@H]2OC3=CC(=C4C(=C3)OC(=CC4=O)C5=CC=C(C=C5)O)O)CO)O)O)O)(CO)O")[[1]]
target.mols <- mols
query.fp <- get.fingerprint(query.mol, type = 'circular')
target.fps <- lapply(target.mols, get.fingerprint, type = 'circular')
sims <- data.frame(
    sim = do.call(rbind, 
                  lapply(
                      target.fps,
                      fingerprint::distance,
                      fp2 = query.fp, 
                      method = 'tanimoto'
                  )
    ))

tani <- tibble::as_tibble(sims, rownames = "smiles") %>%
    dplyr::arrange(desc(sim))

tani.comb <- tani %>%
    dplyr::left_join(
        x = .,
        y = anno,
        by = c("smiles" = "feature")
    ) %>%
    dplyr::filter(sim > 0.20)

tani.comb %>%
    dplyr::select(target_name, sim, target_compound_id, smiles, peak_id, adduct)

hmp.p <- tani.comb %>%
    dplyr::mutate(apigenin = "Apiin") %>%
    dplyr::distinct(target_inchi, .keep_all = TRUE) %>%
    dplyr::mutate(target_name = forcats::fct_reorder(
        .f = target_name,
        .x = sim,
        .fun = "mean",
        .desc = FALSE
    )) %>%
    ggplot(.,
           aes(
               x = apigenin,
               y = target_name,
               fill = sim
           )) +
    geom_tile() +
    scale_y_discrete(expand = expansion(c(0, 0))) +
    scale_x_discrete(expand = expansion(c(0, 0))) +
    scale_fill_gradient() +
    theme_cowplot() +
    theme(
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.7)
    ) +
    labs(
        fill = "Similarity",
        title = "Structural similarity to Apiin"
    )
hmp.p
