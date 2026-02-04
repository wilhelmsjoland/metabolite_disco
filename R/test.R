try <- matched.diffs %>%
    dplyr::filter(
        dplyr::if_any(
            .cols = c("feat1", "feat2"),
            .fns = ~ .x %in% c("FT02089", "FT04480")
        )
    ) %>%
    dplyr::filter(feat1 != feat2) %>%
    dplyr::mutate(obs_ppm = num_to_ppm(abs(delta_mass - obs_delta_mass))) %>%
    dplyr::arrange(obs_ppm) %>%
    dplyr::mutate(
        pair = purrr::map2(feat1, feat2, ~ c(.x, .y)),
        mz1_forms = purrr::map(mz1, ~ Rdisop::getFormula(Rdisop::decomposeMass(.x, ppm = 2))), # or ppm global
        mz2_forms = purrr::map(mz2, ~ Rdisop::getFormula(Rdisop::decomposeMass(.x, ppm = 2))) # or ppm global
    ) # %>%
    # dplyr::select(name, chem_change, feat1, feat2, mz1, mz2, adduct1, adduct2, mass1, mass2, obs_ppm) %>%
    # dplyr::filter(obs_ppm <= 2)

test <- xcms::featureChromatograms(
    object = xchr9,
    expandRt = 0,
    expandMz = 0,
    aggregationFun = "sum",
    filled = TRUE,
    features = unique(c(try$feat1, try$feat2)),
    missing = 0,
    return.type = "XChromatograms"
)

plotFeatPairs(
    feature_chrom = test,
    filt.match.row = try[3,],
    method = "sum",
    value = "into",
    filled = TRUE,
    missing = 0,
    msLevel = 1,
    save_pairs_loc = NULL,
    device = NULL
)

collect.chrs <- xcms::featureChromatograms(
    object = xchr9,
    expandRt = 0,
    expandMz = 0,
    aggregationFun = "sum",
    filled = TRUE,
    features = unique(gly.agly$feature),
    missing = 0,
    return.type = "XChromatograms"
)

for (i in unique(collect.tib$feature)) {
    tmp.p <- plotFeatChrInt(
        feature_chrom = collect.chrs,
        feature = i,
        method = "sum",
        value = "into",
        filled = TRUE,
        missing = "rowmin_half",
        msLevel = 1,
        save_loc = NULL,
        device = NULL,
        feat_pairs = FALSE
    )
    print(tmp.p$combined)
    
    readline("Enter for next: ")
}